#include "openfhe_dtm/PackedElevation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace openfhe_dtm {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEarthRadiusM = 6371008.8;
constexpr int kStitchW = static_cast<int>(kDemStitchWidth);  // 128

double toRad(double d) { return d * kPi / 180.0; }

std::pair<double, double> localVectorM(const GeoPoint& from, const GeoPoint& to) {
    return {toRad(to.lon - from.lon) * kEarthRadiusM * std::cos(toRad(from.lat)),
            toRad(to.lat - from.lat) * kEarthRadiusM};
}

double distanceM(const GeoPoint& a, const GeoPoint& b) {
    const double lat1 = toRad(a.lat), lat2 = toRad(b.lat);
    const double dlat = toRad(b.lat - a.lat), dlon = toRad(b.lon - a.lon);
    const double h = std::sin(dlat / 2) * std::sin(dlat / 2) +
                     std::cos(lat1) * std::cos(lat2) * std::sin(dlon / 2) * std::sin(dlon / 2);
    return 2.0 * kEarthRadiusM * std::asin(std::min(1.0, std::sqrt(h)));
}

GeoPoint offsetPoint(const GeoPoint& o, double east, double north) {
    return {o.lon + east / (kEarthRadiusM * std::cos(toRad(o.lat))) * 180.0 / kPi,
            o.lat + north / kEarthRadiusM * 180.0 / kPi};
}

// Sample line from base to max_distance in the query direction (base point first).
Route sightLine(const SightSurfaceQuery& q) {
    const auto [e, n] = localVectorM(q.base, q.direction);
    const double len = std::hypot(e, n);
    if (len == 0.0) throw std::invalid_argument("direction vector must not be zero length");
    const std::size_t steps = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(q.max_distance_m / std::max(1.0, q.sample_interval_m))));
    Route r;
    for (std::size_t i = 0; i <= steps; ++i) {
        const double d = q.max_distance_m * static_cast<double>(i) / static_cast<double>(steps);
        r.push_back(offsetPoint(q.base, e / len * d, n / len * d));
    }
    return r;
}

// A whole weight vector below double epsilon is replaced by exact zeros
// (OpenFHE rejects a non-zero plaintext below its resolution).
bool allTiny(const std::vector<double>& w) {
    return std::all_of(w.begin(), w.end(), [](double v) {
        return std::abs(v) < std::numeric_limits<double>::epsilon();
    });
}

} // namespace

PackedElevation::PackedElevation(const RasterDtmProvider& dtm, const IHeBackend& he)
    : dtm_(dtm), he_(he) {
    if (he_.slotCount() != kDemPackSlots)
        throw std::invalid_argument("PackedElevation requires exactly 32768 slots");
}

PackedElevation::Cell PackedElevation::locate(const GeoPoint& p) const {
    const auto [east, north] = localVectorM(dtm_.origin(), p);
    double x = east / dtm_.cellSizeMeters(), y = north / dtm_.cellSizeMeters();
    const double max_x = static_cast<double>(dtm_.widthCells() - 1);
    const double max_y = static_cast<double>(dtm_.heightCells() - 1);
    constexpr double eps = 1e-6;
    if (x < -eps || y < -eps || x > max_x + eps || y > max_y + eps)
        throw std::out_of_range("DTM query point is outside raster bounds");
    x = std::clamp(x, 0.0, max_x);
    y = std::clamp(y, 0.0, max_y);
    const int x0 = static_cast<int>(std::floor(x)), y0 = static_cast<int>(std::floor(y));
    Cell c;
    c.tx = x0 / static_cast<int>(kDemTileWidth);
    c.ty = y0 / static_cast<int>(kDemTileWidth);
    const int lx = x0 - c.tx * static_cast<int>(kDemTileWidth);  // 0..63
    const int ly = y0 - c.ty * static_cast<int>(kDemTileWidth);
    c.slot = static_cast<std::size_t>(packSlotOffset(c.tx) + ly * kStitchW + lx);
    c.dx = x - x0;
    c.dy = y - y0;
    return c;
}

std::vector<double> PackedElevation::stitchValues(int tx, int ty) const {
    std::vector<double> out(kDemStitchSlots, 0.0);
    const double nodata = dtm_.metadata().nodata;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i) {
            const DtmTile t = dtm_.loadTile(tx + i, ty + j);
            for (std::size_t y = 0; y < kDemTileWidth; ++y)
                for (std::size_t x = 0; x < kDemTileWidth; ++x) {
                    const double z = t.elevations_m[y * t.width + x];
                    if (z != nodata)
                        out[(j * kDemTileWidth + y) * kDemStitchWidth + i * kDemTileWidth + x] = z / kDemScale;
                }
        }
    return out;
}

void PackedElevation::encryptFor(const Route& route) {
    const auto t0 = std::chrono::steady_clock::now();
    std::set<PackKey> need;
    for (const auto& p : route) need.insert(keyOf(locate(p)));
    for (const auto& key : need) {
        if (packs_.count(key)) continue;
        std::vector<double> pack(kDemPackSlots, 0.0);
        const auto lower = stitchValues(key.first, key.second);      // tile x = packLeftX0
        const auto upper = stitchValues(key.first + 2, key.second);  // tile x + 2
        std::copy(lower.begin(), lower.end(), pack.begin());
        std::copy(upper.begin(), upper.end(), pack.begin() + kDemStitchSlots);
        packs_.emplace(key, he_.encrypt(he_.encode(pack)));
        ++stats_.packs;
    }
    stats_.provider_encrypt_s +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

IndexedCipher PackedElevation::interpolate(const Route& route) const {
    struct Batch { std::vector<std::size_t> points, slots; std::vector<Cell> cells; };
    std::map<PackKey, std::vector<Batch>> groups;
    for (std::size_t i = 0; i < route.size(); ++i) {
        const Cell c = locate(route[i]);
        auto& batches = groups[keyOf(c)];
        Batch* target = nullptr;  // two points in the same cell need different batches
        for (auto& b : batches)
            if (std::find(b.slots.begin(), b.slots.end(), c.slot) == b.slots.end()) { target = &b; break; }
        if (!target) { batches.emplace_back(); target = &batches.back(); }
        target->points.push_back(i);
        target->slots.push_back(c.slot);
        target->cells.push_back(c);
    }

    IndexedCipher out;
    out.point_count = route.size();
    const std::size_t S = he_.slotCount();
    for (auto& [key, batches] : groups) {
        const auto it = packs_.find(key);
        if (it == packs_.end()) throw std::logic_error("packed stitch was not encrypted by the provider");
        const CipherVector& ct = it->second;
        // Neighbour alignment, once per pack and shared by all its batches.
        const CipherVector r10 = he_.rotate(ct, 1);
        const CipherVector r01 = he_.rotate(ct, kStitchW);
        const CipherVector r11 = he_.rotate(ct, kStitchW + 1);
        stats_.rotations += 3;
        for (auto& b : batches) {
            std::vector<double> w00(S, 0.0), w10(S, 0.0), w01(S, 0.0), w11(S, 0.0);
            for (std::size_t k = 0; k < b.cells.size(); ++k) {
                const double dx = b.cells[k].dx, dy = b.cells[k].dy;
                const std::size_t s = b.slots[k];
                w00[s] = (1 - dx) * (1 - dy) * kDemScale;
                w10[s] = dx * (1 - dy) * kDemScale;
                w01[s] = (1 - dx) * dy * kDemScale;
                w11[s] = dx * dy * kDemScale;
            }
            CipherVector acc;
            bool first = true;
            const std::pair<const CipherVector*, std::vector<double>*> terms[4] = {
                {&ct, &w00}, {&r10, &w10}, {&r01, &w01}, {&r11, &w11}};
            for (const auto& [c, w] : terms) {
                if (allTiny(*w)) continue;
                CipherVector t = he_.mulPlain(*c, he_.encode(*w));
                acc = first ? t : he_.add(acc, t);
                first = false;
            }
            if (first) acc = he_.mulPlain(ct, he_.encode(std::vector<double>(S, 0.0)));
            out.batches.push_back({acc, b.slots, b.points});
            ++stats_.batches;
        }
    }
    return out;
}

IndexedCipher PackedElevation::mergeDisjoint(const IndexedCipher& in) const {
    IndexedCipher out;
    out.point_count = in.point_count;
    std::vector<std::set<std::size_t>> used;
    for (const auto& b : in.batches) {
        std::size_t target = out.batches.size();
        for (std::size_t m = 0; m < out.batches.size(); ++m)
            if (std::none_of(b.slots.begin(), b.slots.end(), [&](std::size_t s) { return used[m].count(s); })) {
                target = m; break;
            }
        if (target == out.batches.size()) { out.batches.push_back(b); used.emplace_back(b.slots.begin(), b.slots.end()); continue; }
        auto& t = out.batches[target];
        t.ciphertext = he_.add(t.ciphertext, b.ciphertext);
        t.slots.insert(t.slots.end(), b.slots.begin(), b.slots.end());
        t.points.insert(t.points.end(), b.points.begin(), b.points.end());
        used[target].insert(b.slots.begin(), b.slots.end());
    }
    return out;
}

std::vector<double> PackedElevation::decrypt(const IndexedCipher& c, std::size_t* decryptions) const {
    std::vector<double> out(c.point_count, std::numeric_limits<double>::quiet_NaN());
    for (const auto& b : c.batches) {
        const PlainVector v = he_.decrypt(b.ciphertext);
        for (std::size_t k = 0; k < b.points.size(); ++k) out[b.points[k]] = v.values.at(b.slots[k]);
        if (decryptions) ++*decryptions;
    }
    return out;
}

PackedSightResult sightSurfacePacked(const SightSurfaceQuery& query, PackedElevation& elevation) {
    if (!(query.angle_degrees > -90.0 && query.angle_degrees < 90.0))
        throw std::invalid_argument("angle must be strictly between -90 and 90");
    const IHeBackend& he = elevation.backend();
    const Route points = sightLine(query);
    PackedSightResult result;
    result.distances_m.assign(points.size(), 0.0);
    for (std::size_t i = 1; i < points.size(); ++i)
        result.distances_m[i] = result.distances_m[i - 1] + distanceM(points[i - 1], points[i]);

    elevation.encryptFor(points);                       // provider
    // Batches are zero outside their own slots, so disjoint ones can be added first.
    const IndexedCipher ground = elevation.mergeDisjoint(elevation.interpolate(points));  // G(d), point 0 = base

    // z_B: mask its slot, then a slot sum copies it into every slot.
    CipherVector zb;
    bool found = false;
    const std::size_t S = he.slotCount();
    for (const auto& b : ground.batches)
        for (std::size_t k = 0; k < b.points.size() && !found; ++k)
            if (b.points[k] == 0) {
                std::vector<double> mask(S, 0.0);
                mask[b.slots[k]] = 1.0;
                zb = he.sumSlots(he.mulPlain(b.ciphertext, he.encode(mask)));
                found = true;
            }
    if (!found) throw std::logic_error("base point missing from interpolation");

    const double tan_theta = std::tan(toRad(query.angle_degrees));
    IndexedCipher h;
    h.point_count = points.size();
    for (const auto& b : ground.batches) {
        std::vector<double> neg(S, 0.0);  // -(d tan theta) at each point's slot
        for (std::size_t k = 0; k < b.points.size(); ++k)
            neg[b.slots[k]] = -result.distances_m[b.points[k]] * tan_theta;
        const CipherVector plane = he.subPlain(zb, he.encode(neg));  // z_B + d tan theta
        h.batches.push_back({he.sub(plane, b.ciphertext), b.slots, b.points});
    }
    const std::vector<double> diff = elevation.decrypt(h, &result.decryptions);
    result.buildable_height_m.resize(diff.size());
    for (std::size_t i = 0; i < diff.size(); ++i) result.buildable_height_m[i] = std::max(0.0, diff[i]);
    return result;
}

} // namespace openfhe_dtm

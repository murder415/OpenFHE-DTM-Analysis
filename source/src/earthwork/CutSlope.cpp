#include "openfhe_dtm/CutSlope.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace cutslope {
namespace {
constexpr double kEarthRadiusM = 6371008.8;
constexpr double kPi = 3.14159265358979323846;
double toRad(double d) { return d * kPi / 180.0; }
using XY = std::pair<double, double>;
XY local(const GeoPoint& o, const GeoPoint& p) {
    return {toRad(p.lon - o.lon) * kEarthRadiusM * std::cos(toRad(o.lat)), toRad(p.lat - o.lat) * kEarthRadiusM};
}
GeoPoint geo(const GeoPoint& o, XY p) {
    return {o.lon + p.first / (kEarthRadiusM * std::cos(toRad(o.lat))) * 180.0 / kPi,
            o.lat + p.second / kEarthRadiusM * 180.0 / kPi};
}
bool inside(const std::vector<XY>& ring, XY p) {
    int w = 0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const XY a = ring[i], b = ring[(i + 1) % ring.size()];
        const double c = (b.first - a.first) * (p.second - a.second) - (p.first - a.first) * (b.second - a.second);
        if (a.second <= p.second && b.second > p.second && c > 0) ++w;
        if (a.second > p.second && b.second <= p.second && c < 0) --w;
    }
    return w != 0;
}
double boundaryDistance(const std::vector<XY>& ring, XY p) {
    double best = 1e300;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const XY a = ring[i], b = ring[(i + 1) % ring.size()];
        const double vx = b.first - a.first, vy = b.second - a.second;
        const double L2 = vx * vx + vy * vy;
        double t = L2 > 0 ? ((p.first - a.first) * vx + (p.second - a.second) * vy) / L2 : 0;
        t = std::clamp(t, 0.0, 1.0);
        best = std::min(best, std::hypot(p.first - (a.first + t * vx), p.second - (a.second + t * vy)));
    }
    return best;
}
} // namespace

CutSlopeGeometry buildGeometry(const CutSlopeQuery& q) {
    if (q.polygon.size() < 3) throw std::invalid_argument("polygon needs 3 points");
    if (!(q.cut_angle_degrees > 0 && q.cut_angle_degrees < 90)) throw std::invalid_argument("0 < theta < 90");
    if (!(q.grid_interval_m > 0) || !(q.slope_width_m >= 0) || !std::isfinite(q.design_height_m))
        throw std::invalid_argument("bad interval/width/H");
    const GeoPoint o = q.polygon.front();
    std::vector<XY> ring;
    for (const auto& p : q.polygon) ring.push_back(local(o, p));
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
    for (auto [x, y] : ring) { x0 = std::min(x0, x); x1 = std::max(x1, x); y0 = std::min(y0, y); y1 = std::max(y1, y); }
    x0 -= q.slope_width_m; x1 += q.slope_width_m; y0 -= q.slope_width_m; y1 += q.slope_width_m;
    const std::size_t cols = static_cast<std::size_t>(std::ceil((x1 - x0) / q.grid_interval_m));
    const std::size_t rows = static_cast<std::size_t>(std::ceil((y1 - y0) / q.grid_interval_m));
    CutSlopeGeometry g;
    g.dx = (x1 - x0) / cols; g.dy = (y1 - y0) / rows;
    const double quarter = g.dx * g.dy / 4.0, tanv = std::tan(toRad(q.cut_angle_degrees));
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> index;
    auto vertex = [&](std::size_t ix, std::size_t iy) {
        auto it = index.find({ix, iy});
        if (it != index.end()) return it->second;
        const XY p{x0 + ix * g.dx, y0 + iy * g.dy};
        const bool in = inside(ring, p);
        const double d = in ? 0.0 : boundaryDistance(ring, p);
        const std::size_t k = g.points.size();
        index[{ix, iy}] = k;
        g.points.push_back(geo(o, p));
        g.inside_w.push_back(0); g.outside_w.push_back(0);
        g.boundary_distance.push_back(d);
        g.slope_height.push_back(q.design_height_m + d * tanv);
        g.outer_edge.push_back(false);
        return k;
    };
    for (std::size_t iy = 0; iy < rows; ++iy)
        for (std::size_t ix = 0; ix < cols; ++ix) {
            const XY c{x0 + (ix + 0.5) * g.dx, y0 + (iy + 0.5) * g.dy};
            const bool in = inside(ring, c);
            const double d = in ? 0.0 : boundaryDistance(ring, c);
            if (!in && d > q.slope_width_m) continue;
            const std::array<std::size_t, 4> v{vertex(ix, iy), vertex(ix + 1, iy), vertex(ix, iy + 1), vertex(ix + 1, iy + 1)};
            for (auto k : v) (in ? g.inside_w : g.outside_w)[k] += quarter;
            if (in) { ++g.inside_cells; g.inside_area += 4 * quarter; }
            else    { ++g.outside_cells; g.outside_area += 4 * quarter; }
        }
    for (std::size_t k = 0; k < g.points.size(); ++k)
        g.outer_edge[k] = g.boundary_distance[k] > q.slope_width_m - std::hypot(g.dx, g.dy);
    return g;
}

CutSlopeResult cutSlopeEarthworkPacked(const CutSlopeQuery& q, PackedElevation& elevation) {
    const IHeBackend& he = elevation.backend();
    CutSlopeResult r;
    r.geometry = buildGeometry(q);
    const auto& g = r.geometry;
    const std::size_t S = he.slotCount();
    elevation.encryptFor(g.points);                                          // provider
    const IndexedCipher cz = elevation.mergeDisjoint(elevation.interpolate(g.points));  // ground z
    const CipherVector cH = he.encrypt(he.encode(std::vector<double>(S, q.design_height_m)));  // H
    CipherVector cin;
    bool first = true;
    IndexedCipher cout;
    cout.point_count = g.points.size();
    for (const auto& b : cz.batches) {
        std::vector<double> Ain(S, 0.0), Aout(S, 0.0), s(S, 0.0);
        for (std::size_t k = 0; k < b.points.size(); ++k) {
            Ain[b.slots[k]] = g.inside_w[b.points[k]];
            Aout[b.slots[k]] = g.outside_w[b.points[k]];
            s[b.slots[k]] = g.slope_height[b.points[k]];
        }
        const CipherVector part = he.sumSlots(he.mulPlain(he.sub(cH, b.ciphertext), he.encode(Ain)));
        cin = first ? part : he.add(cin, part);                             // sum (H - z) A_in
        first = false;
        cout.batches.push_back({he.mulPlain(he.subPlain(b.ciphertext, he.encode(s)), he.encode(Aout)),
                                b.slots, b.points});                         // (z - s) A_out
    }
    // Owner side: one scalar decryption plus one per merged batch.
    const long double in_sum = he.decrypt(cin).values.at(0);
    const std::vector<double> vout = elevation.decrypt(cout, &r.decryptions);
    r.decryptions += 1;
    r.batches = cz.batches.size();
    long double cut = 0;
    for (std::size_t i = 0; i < vout.size(); ++i) {
        if (vout[i] > 0) cut += vout[i];
        if (g.outer_edge[i] && g.outside_w[i] > 0 && vout[i] > 1e-9 * g.outside_w[i]) ++r.unfinished_slope_points;
    }
    r.inside_fill_minus_cut_m3 = static_cast<double>(in_sum);
    r.slope_cut_m3 = static_cast<double>(cut);
    r.net_fill_minus_cut_m3 = static_cast<double>(in_sum - cut);
    return r;
}
} // namespace cutslope

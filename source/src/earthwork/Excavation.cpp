#include "openfhe_dtm/Analysis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <stdexcept>

namespace openfhe_dtm {
namespace {

using XY = std::pair<double, double>;
constexpr double kPi = 3.14159265358979323846;
constexpr double kEarthRadiusM = 6371008.8;

double toRad(double degrees) { return degrees * kPi / 180.0; }

bool samePoint(const GeoPoint& a, const GeoPoint& b) {
    return a.lon == b.lon && a.lat == b.lat;
}

void validatePoint(const GeoPoint& p) {
    if (!std::isfinite(p.lon) || !std::isfinite(p.lat) ||
        p.lon < -180.0 || p.lon > 180.0 || p.lat <= -90.0 || p.lat >= 90.0) {
        throw std::invalid_argument("excavation polygon requires finite geographic coordinates");
    }
}

XY localPoint(const GeoPoint& origin, const GeoPoint& point) {
    return {toRad(point.lon - origin.lon) * kEarthRadiusM * std::cos(toRad(origin.lat)),
            toRad(point.lat - origin.lat) * kEarthRadiusM};
}

GeoPoint geographicPoint(const GeoPoint& origin, const XY& p) {
    return {origin.lon + (p.first / (kEarthRadiusM * std::cos(toRad(origin.lat)))) * 180.0 / kPi,
            origin.lat + (p.second / kEarthRadiusM) * 180.0 / kPi};
}

long double cross(const XY& a, const XY& b, const XY& c) {
    return (static_cast<long double>(b.first) - a.first) *
               (static_cast<long double>(c.second) - a.second) -
           (static_cast<long double>(b.second) - a.second) *
               (static_cast<long double>(c.first) - a.first);
}

bool onSegment(const XY& a, const XY& b, const XY& p) {
    return cross(a, b, p) == 0.0L &&
           p.first >= std::min(a.first, b.first) && p.first <= std::max(a.first, b.first) &&
           p.second >= std::min(a.second, b.second) && p.second <= std::max(a.second, b.second);
}

bool segmentsIntersect(const XY& a, const XY& b, const XY& c, const XY& d) {
    const auto ab_c = cross(a, b, c);
    const auto ab_d = cross(a, b, d);
    const auto cd_a = cross(c, d, a);
    const auto cd_b = cross(c, d, b);
    const bool opposite_ab = (ab_c < 0.0L && ab_d > 0.0L) || (ab_c > 0.0L && ab_d < 0.0L);
    const bool opposite_cd = (cd_a < 0.0L && cd_b > 0.0L) || (cd_a > 0.0L && cd_b < 0.0L);
    return (opposite_ab && opposite_cd) || onSegment(a, b, c) || onSegment(a, b, d) ||
           onSegment(c, d, a) || onSegment(c, d, b);
}

void validateRing(const std::vector<XY>& p) {
    long double twice_area = 0.0L;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const std::size_t next = (i + 1) % p.size();
        if (p[i] == p[next]) {
            throw std::invalid_argument("excavation polygon has a zero-length edge");
        }
        twice_area += static_cast<long double>(p[i].first) * p[next].second -
                      static_cast<long double>(p[next].first) * p[i].second;
        const auto& prev = p[(i + p.size() - 1) % p.size()];
        if (cross(prev, p[i], p[next]) == 0.0L) {
            const long double dot =
                (static_cast<long double>(prev.first) - p[i].first) * (p[next].first - p[i].first) +
                (static_cast<long double>(prev.second) - p[i].second) * (p[next].second - p[i].second);
            if (dot > 0.0L) {
                throw std::invalid_argument("excavation polygon has overlapping adjacent edges");
            }
        }
        for (std::size_t j = i + 1; j < p.size(); ++j) {
            const std::size_t jnext = (j + 1) % p.size();
            if (j == next || jnext == i) continue;
            if (segmentsIntersect(p[i], p[next], p[j], p[jnext])) {
                throw std::invalid_argument("excavation polygon must be a simple ring");
            }
        }
    }
    if (twice_area == 0.0L) {
        throw std::invalid_argument("excavation polygon has zero area");
    }
}

bool pointInPolygon(const std::vector<XY>& polygon, double x, double y) {
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const auto [xi, yi] = polygon[i];
        const auto [xj, yj] = polygon[j];
        const bool crosses = ((yi > y) != (yj > y)) &&
                             (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (crosses) inside = !inside;
    }
    return inside;
}


struct PreparedGeometry {
    EarthworkSlopeGeometry public_data;
    std::vector<XY> volume_xy;
    std::vector<XY> boundary_xy;
};

std::size_t subdivisions(double length, double interval, std::size_t limit) {
    const long double count = std::ceil(static_cast<long double>(length) / interval);
    if (!std::isfinite(count) || count < 1.0L || count > static_cast<long double>(limit)) {
        throw std::length_error("earthwork slope subdivision count exceeds its limit");
    }
    return static_cast<std::size_t>(count);
}

void validateQuery(const EarthworkSlopeQuery& query, std::size_t slots) {
    if (slots == 0) throw std::invalid_argument("earthwork slope requires nonzero slotCount");
    if (!std::isfinite(query.design_height_m)) {
        throw std::invalid_argument("earthwork slope requires an explicit finite design height");
    }
    if (!std::isfinite(query.slope_angle_degrees) ||
        query.slope_angle_degrees <= 0.0 || query.slope_angle_degrees >= 90.0) {
        throw std::invalid_argument("earthwork slope angle must be finite and strictly between 0 and 90 degrees");
    }
    if (!std::isfinite(query.sample_interval_m) || query.sample_interval_m <= 0.0 ||
        !std::isfinite(query.boundary_interval_m) || query.boundary_interval_m < 0.0) {
        throw std::invalid_argument("earthwork slope sampling intervals are invalid");
    }
    if (!std::isfinite(query.area_threshold_m) || query.area_threshold_m < 0.0 ||
        !std::isfinite(query.envelope_tolerance_m) || query.envelope_tolerance_m < 0.0) {
        throw std::invalid_argument("earthwork slope height tolerances must be finite and nonnegative");
    }
    if (query.max_grid_cells == 0 || query.max_boundary_points == 0 ||
        query.max_candidate_ciphertexts == 0) {
        throw std::invalid_argument("earthwork slope resource limits must be positive");
    }
}

PreparedGeometry prepareGeometry(const EarthworkSlopeQuery& query) {
    Route ring = query.polygon;
    if (ring.size() > 1 && samePoint(ring.front(), ring.back())) ring.pop_back();
    if (ring.size() < 3) throw std::invalid_argument("earthwork slope polygon needs at least three vertices");
    for (const auto& p : ring) validatePoint(p);
    if (std::abs(std::cos(toRad(ring.front().lat))) < 1e-8) {
        throw std::invalid_argument("earthwork slope local projection is undefined near the poles");
    }

    PreparedGeometry prepared;
    auto& geometry = prepared.public_data;
    geometry.origin = ring.front();
    std::vector<XY> polygon_xy;
    for (const auto& p : ring) {
        if (std::abs(p.lon - geometry.origin.lon) >= 180.0) {
            throw std::invalid_argument("earthwork slope local projection does not support longitude wrap");
        }
        polygon_xy.push_back(localPoint(geometry.origin, p));
    }
    validateRing(polygon_xy);
    double min_x = polygon_xy.front().first, max_x = min_x;
    double min_y = polygon_xy.front().second, max_y = min_y;
    for (const auto& p : polygon_xy) {
        min_x = std::min(min_x, p.first); max_x = std::max(max_x, p.first);
        min_y = std::min(min_y, p.second); max_y = std::max(max_y, p.second);
    }
    geometry.grid_columns = subdivisions(max_x - min_x, query.sample_interval_m, query.max_grid_cells);
    geometry.grid_rows = subdivisions(max_y - min_y, query.sample_interval_m, query.max_grid_cells);
    if (geometry.grid_columns > query.max_grid_cells / geometry.grid_rows) {
        throw std::length_error("earthwork slope bounding-box grid exceeds max_grid_cells");
    }
    geometry.grid_min_east_m = min_x;
    geometry.grid_min_north_m = min_y;
    geometry.grid_dx_m = (max_x - min_x) / static_cast<double>(geometry.grid_columns);
    geometry.grid_dy_m = (max_y - min_y) / static_cast<double>(geometry.grid_rows);
    const double cell_area = geometry.grid_dx_m * geometry.grid_dy_m;
    const double corner_area = cell_area / 4.0;
    if (!std::isfinite(corner_area) || corner_area <= 0.0) {
        throw std::invalid_argument("earthwork slope quadrature weight is not representable");
    }

    std::map<std::pair<std::size_t, std::size_t>, std::size_t> node_indices;
    const auto add_node = [&](std::size_t ix, std::size_t iy) {
        const auto key = std::make_pair(ix, iy);
        const auto found = node_indices.find(key);
        if (found != node_indices.end()) {
            geometry.point_weights_m2[found->second] += corner_area;
            return found->second;
        }
        const XY xy{
            ix == geometry.grid_columns ? max_x : min_x + static_cast<double>(ix) * geometry.grid_dx_m,
            iy == geometry.grid_rows ? max_y : min_y + static_cast<double>(iy) * geometry.grid_dy_m};
        const std::size_t index = geometry.points.size();
        node_indices.emplace(key, index);
        prepared.volume_xy.push_back(xy);
        geometry.points.push_back(geographicPoint(geometry.origin, xy));
        geometry.point_weights_m2.push_back(corner_area);
        return index;
    };
    for (std::size_t iy = 0; iy < geometry.grid_rows; ++iy) {
        for (std::size_t ix = 0; ix < geometry.grid_columns; ++ix) {
            const double x = min_x + (static_cast<double>(ix) + 0.5) * geometry.grid_dx_m;
            const double y = min_y + (static_cast<double>(iy) + 0.5) * geometry.grid_dy_m;
            if (!pointInPolygon(polygon_xy, x, y)) continue;
            geometry.cells.push_back({add_node(ix, iy), add_node(ix + 1, iy),
                                      add_node(ix + 1, iy + 1), add_node(ix, iy + 1)});
        }
    }
    if (geometry.cells.empty()) {
        throw std::invalid_argument("earthwork slope polygon produced no center-inside grid cells");
    }
    geometry.integration_area_m2 = cell_area * static_cast<double>(geometry.cells.size());
    if (!std::isfinite(geometry.integration_area_m2)) {
        throw std::invalid_argument("earthwork slope integration area is not representable");
    }

    const double boundary_interval = query.boundary_interval_m == 0.0 ?
                                     query.sample_interval_m : query.boundary_interval_m;
    for (std::size_t edge = 0; edge < polygon_xy.size(); ++edge) {
        const auto& a = polygon_xy[edge];
        const auto& b = polygon_xy[(edge + 1) % polygon_xy.size()];
        const double length = std::hypot(b.first - a.first, b.second - a.second);
        const std::size_t steps = subdivisions(length, boundary_interval, query.max_boundary_points);
        if (steps > query.max_boundary_points - geometry.boundary_points.size()) {
            throw std::length_error("earthwork slope boundary sampling exceeds max_boundary_points");
        }
        for (std::size_t j = 0; j < steps; ++j) {
            const double t = static_cast<double>(j) / static_cast<double>(steps);
            const XY xy{a.first + t * (b.first - a.first), a.second + t * (b.second - a.second)};
            prepared.boundary_xy.push_back(xy);
            geometry.boundary_points.push_back(geographicPoint(geometry.origin, xy));
        }
    }
    return prepared;
}

std::vector<double> finalValues(const IHeBackend& he, const CipherVector& cipher,
                                std::size_t valid_count) {
    auto values = he.decrypt(cipher).values;
    if (values.size() < valid_count) {
        throw std::runtime_error("earthwork slope decryption returned too few valid slots");
    }
    values.resize(valid_count);
    for (double v : values) {
        if (!std::isfinite(v)) throw std::runtime_error("earthwork slope decrypted a nonfinite candidate");
    }
    return values;
}

void validateEncryptedShape(const EarthworkSlopeEncryptedResult& encrypted, std::size_t slots) {
    const auto& g = encrypted.geometry;
    if (slots == 0 || g.points.empty() || g.boundary_points.empty() ||
        g.point_weights_m2.size() != g.points.size() || g.cells.empty() ||
        !std::isfinite(g.integration_area_m2) || g.integration_area_m2 <= 0.0 ||
        !std::isfinite(encrypted.area_threshold_m) || encrypted.area_threshold_m < 0.0 ||
        !std::isfinite(encrypted.envelope_tolerance_m) || encrypted.envelope_tolerance_m < 0.0) {
        throw std::invalid_argument("earthwork slope encrypted result has invalid metadata");
    }
    for (double weight : g.point_weights_m2) {
        if (!std::isfinite(weight) || weight <= 0.0) {
            throw std::invalid_argument("earthwork slope quadrature weights must be finite and positive");
        }
    }
    for (const auto& cell : g.cells) {
        for (std::size_t point : cell) {
            if (point >= g.points.size()) throw std::invalid_argument("earthwork slope cell index is invalid");
        }
    }
    if (g.points.size() > std::numeric_limits<std::size_t>::max() - g.boundary_points.size()) {
        throw std::invalid_argument("earthwork slope target count overflow");
    }
    const std::size_t total = g.points.size() + g.boundary_points.size();
    std::size_t expected_begin = 0;
    for (const auto& batch : encrypted.batches) {
        if (batch.begin != expected_begin || batch.sample_count == 0 ||
            batch.sample_count > slots || batch.sample_count > total - expected_begin ||
            batch.upper_candidates.size() != g.boundary_points.size() ||
            batch.lower_candidates.size() != g.boundary_points.size()) {
            throw std::invalid_argument("earthwork slope ciphertext batch layout is invalid");
        }
        expected_begin += batch.sample_count;
    }
    if (expected_begin != total) {
        throw std::invalid_argument("earthwork slope ciphertext batches do not cover every target");
    }
}

} // namespace

EarthworkSlopeEncryptedResult TerrainAnalysis::earthworkSlopeEncrypted(const EarthworkSlopeQuery& query) const {
    const std::size_t slots = he_.slotCount();
    validateQuery(query, slots);
    auto prepared = prepareGeometry(query);
    auto& geometry = prepared.public_data;
    Route targets = geometry.points;
    targets.insert(targets.end(), geometry.boundary_points.begin(), geometry.boundary_points.end());
    std::vector<XY> target_xy = prepared.volume_xy;
    target_xy.insert(target_xy.end(), prepared.boundary_xy.begin(), prepared.boundary_xy.end());
    const std::size_t batch_count = targets.size() / slots + (targets.size() % slots != 0);
    const std::size_t boundary_count = geometry.boundary_points.size();
    if (boundary_count > (std::numeric_limits<std::size_t>::max() - 1) / 2 ||
        batch_count > query.max_candidate_ciphertexts / (2 * boundary_count + 1)) {
        throw std::length_error("earthwork slope candidate ciphertext count exceeds its limit");
    }

    EarthworkSlopeEncryptedResult encrypted;
    encrypted.geometry = geometry;
    encrypted.area_threshold_m = query.area_threshold_m;
    encrypted.envelope_tolerance_m = query.envelope_tolerance_m;
    const double slope = std::tan(toRad(query.slope_angle_degrees));
    if (!std::isfinite(slope)) throw std::invalid_argument("earthwork slope tangent is not finite");

    // All target batches and all boundary candidates are completed here. There
    // are deliberately no decrypt, elevation(), or dtm.elevationAt() calls.
    for (std::size_t begin = 0; begin < targets.size();) {
        const std::size_t count = std::min(slots, targets.size() - begin);
        const Route target_batch(targets.begin() + begin, targets.begin() + begin + count);
        // Keep diagnostic slots on the same magnitude as volume slots. A
        // weight of 1 beside large cell areas amplifies CKKS absolute noise
        // when diagnostic candidates are converted back to height units.
        const double diagnostic_weight = geometry.grid_dx_m * geometry.grid_dy_m / 4.0;
        std::vector<double> weights(count, diagnostic_weight);
        for (std::size_t i = 0; i < count; ++i) {
            if (begin + i < geometry.points.size()) weights[i] = geometry.point_weights_m2[begin + i];
        }
        const CipherVector ground = elevationEncrypted(target_batch);
        const CipherVector height = he_.encrypt(he_.encode(std::vector<double>(count, query.design_height_m)));
        const PlainVector area = he_.encode(weights);
        EarthworkSlopeCipherBatch batch;
        batch.begin = begin;
        batch.sample_count = count;
        batch.target_minus_ground = he_.mulPlain(he_.sub(height, ground), area);
        batch.upper_candidates.reserve(boundary_count);
        batch.lower_candidates.reserve(boundary_count);
        for (std::size_t b = 0; b < boundary_count; ++b) {
            // Repeated public locations reuse the same bilinear HE primitive.
            // This avoids decrypting a boundary elevation for broadcasting.
            const CipherVector boundary = elevationEncrypted(Route(count, geometry.boundary_points[b]));
            const CipherVector difference = he_.sub(boundary, ground);
            std::vector<double> correction(count);
            std::vector<double> negative_correction(count);
            for (std::size_t i = 0; i < count; ++i) {
                const auto& xy = target_xy[begin + i];
                const auto& bxy = prepared.boundary_xy[b];
                correction[i] = slope * std::hypot(xy.first - bxy.first, xy.second - bxy.second);
                if (!std::isfinite(correction[i])) {
                    throw std::invalid_argument("earthwork slope cone correction is not finite");
                }
                negative_correction[i] = -correction[i];
            }
            batch.upper_candidates.push_back(
                he_.mulPlain(he_.subPlain(difference, he_.encode(negative_correction)), area));
            batch.lower_candidates.push_back(
                he_.mulPlain(he_.subPlain(difference, he_.encode(correction)), area));
        }
        encrypted.batches.push_back(std::move(batch));
        begin += count;
    }
    return encrypted;
}

EarthworkSlopeResult TerrainAnalysis::finalizeEarthworkSlope(const EarthworkSlopeEncryptedResult& encrypted) const {
    validateEncryptedShape(encrypted, he_.slotCount());
    EarthworkSlopeResult result;
    result.geometry = encrypted.geometry;
    const auto& geometry = result.geometry;
    const std::size_t volume_count = geometry.points.size();
    std::vector<double> volumes(volume_count);
    std::vector<bool> infeasible(volume_count, false);
    result.point_height_changes_m.resize(volume_count);

    // This is the final owner-side comparison stage. A positive common weight
    // and subtracting the same ground height commute with min/max. Thus:
    // max_b(G_b), min_b(F_b), U reproduce weight*(max(g,min(H,f))-ground).
    for (const auto& batch : encrypted.batches) {
        const auto target = finalValues(he_, batch.target_minus_ground, batch.sample_count);
        std::vector<double> upper(batch.sample_count, std::numeric_limits<double>::infinity());
        std::vector<double> lower(batch.sample_count, -std::numeric_limits<double>::infinity());
        for (std::size_t b = 0; b < geometry.boundary_points.size(); ++b) {
            const auto candidate_upper = finalValues(he_, batch.upper_candidates[b], batch.sample_count);
            const auto candidate_lower = finalValues(he_, batch.lower_candidates[b], batch.sample_count);
            for (std::size_t i = 0; i < batch.sample_count; ++i) {
                upper[i] = std::min(upper[i], candidate_upper[i]);
                lower[i] = std::max(lower[i], candidate_lower[i]);
            }
        }
        for (std::size_t i = 0; i < batch.sample_count; ++i) {
            const std::size_t index = batch.begin + i;
            const bool volume_target = index < volume_count;
            const double weight = volume_target ? geometry.point_weights_m2[index]
                : geometry.grid_dx_m * geometry.grid_dy_m / 4.0;
            const double gap = (lower[i] - upper[i]) / weight;
            if (!std::isfinite(gap)) throw std::runtime_error("earthwork slope envelope gap overflow");
            result.max_envelope_gap_m = std::max(result.max_envelope_gap_m, gap);
            if (gap > encrypted.envelope_tolerance_m) {
                if (volume_target) {
                    infeasible[index] = true;
                    ++result.infeasible_point_count;
                } else {
                    ++result.infeasible_boundary_point_count;
                }
            }
            if (volume_target) {
                const double volume = std::max(lower[i], std::min(target[i], upper[i]));
                const double change = volume / weight;
                if (!std::isfinite(volume) || !std::isfinite(change)) {
                    throw std::runtime_error("earthwork slope final height/volume is not finite");
                }
                volumes[index] = volume;
                result.point_height_changes_m[index] = change;
            }
        }
    }
    for (const auto& cell : geometry.cells) {
        if (std::any_of(cell.begin(), cell.end(), [&](std::size_t p) { return infeasible[p]; })) {
            ++result.infeasible_cell_count;
        }
    }
    result.envelope_valid = result.infeasible_point_count == 0 && result.infeasible_boundary_point_count == 0;
    if (!result.envelope_valid) {
        result.status = "invalid_boundary_envelope";
        // No engineering volume is reported for contradictory constraints.
        result.point_height_changes_m.clear();
        return result;
    }

    // The tolerance only decides whether a tiny positive gap is treated as
    // numerical inconsistency. It never alters the candidates or clips volume.
    long double cut = 0.0L, fill = 0.0L, cut_area = 0.0L, fill_area = 0.0L, unchanged_area = 0.0L;
    for (std::size_t i = 0; i < volume_count; ++i) {
        const double volume = volumes[i];
        const double delta = result.point_height_changes_m[i];
        if (volume > 0.0) fill += volume;
        else cut -= volume;
        if (delta > encrypted.area_threshold_m) fill_area += geometry.point_weights_m2[i];
        else if (delta < -encrypted.area_threshold_m) cut_area += geometry.point_weights_m2[i];
        else unchanged_area += geometry.point_weights_m2[i];
    }
    result.cut_volume_m3 = static_cast<double>(cut);
    result.fill_volume_m3 = static_cast<double>(fill);
    result.net_volume_m3 = static_cast<double>(fill - cut);
    result.absolute_moved_volume_m3 = static_cast<double>(fill + cut);
    result.cut_area_m2 = static_cast<double>(cut_area);
    result.fill_area_m2 = static_cast<double>(fill_area);
    result.unchanged_area_m2 = static_cast<double>(unchanged_area);
    if (!std::isfinite(result.cut_volume_m3) || !std::isfinite(result.fill_volume_m3) ||
        !std::isfinite(result.net_volume_m3) || !std::isfinite(result.absolute_moved_volume_m3) ||
        !std::isfinite(result.cut_area_m2) || !std::isfinite(result.fill_area_m2) ||
        !std::isfinite(result.unchanged_area_m2)) {
        throw std::runtime_error("earthwork slope aggregate overflow");
    }
    result.status = "valid_sampled_boundary_envelope";
    return result;
}

EarthworkSlopeResult TerrainAnalysis::earthworkSlope(const EarthworkSlopeQuery& query) const {
    const auto encrypted = earthworkSlopeEncrypted(query);
    return finalizeEarthworkSlope(encrypted);
}

} // namespace openfhe_dtm

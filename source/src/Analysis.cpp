#include "openfhe_dtm/Analysis.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace openfhe_dtm {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEarthRadiusM = 6371008.8;

double toRad(double degrees) {
    return degrees * kPi / 180.0;
}

double distanceM(const GeoPoint& a, const GeoPoint& b) {
    const double lat1 = toRad(a.lat);
    const double lat2 = toRad(b.lat);
    const double dlat = toRad(b.lat - a.lat);
    const double dlon = toRad(b.lon - a.lon);
    const double h = std::sin(dlat / 2.0) * std::sin(dlat / 2.0) +
                     std::cos(lat1) * std::cos(lat2) *
                         std::sin(dlon / 2.0) * std::sin(dlon / 2.0);
    return 2.0 * kEarthRadiusM * std::asin(std::min(1.0, std::sqrt(h)));
}

GeoPoint offsetPoint(const GeoPoint& origin, double east_m, double north_m) {
    const double lat_rad = toRad(origin.lat);
    return GeoPoint{
        origin.lon + (east_m / (kEarthRadiusM * std::cos(lat_rad))) * 180.0 / kPi,
        origin.lat + (north_m / kEarthRadiusM) * 180.0 / kPi,
    };
}

std::pair<double, double> localVectorM(const GeoPoint& from, const GeoPoint& to) {
    const double lat_rad = toRad(from.lat);
    return {
        toRad(to.lon - from.lon) * kEarthRadiusM * std::cos(lat_rad),
        toRad(to.lat - from.lat) * kEarthRadiusM,
    };
}

std::vector<double> cumulativeDistances(const Route& route) {
    std::vector<double> distances(route.size(), 0.0);
    for (std::size_t i = 1; i < route.size(); ++i) {
        distances[i] = distances[i - 1] + distanceM(route[i - 1], route[i]);
    }
    return distances;
}

double routeLengthM(const Route& route) {
    const std::vector<double> distances = cumulativeDistances(route);
    return distances.empty() ? 0.0 : distances.back();
}

GeoPoint pointAtDistance(const Route& route, double target_m) {
    if (route.empty()) {
        throw std::invalid_argument("route must not be empty");
    }
    if (route.size() == 1 || target_m <= 0.0) {
        return route.front();
    }

    double covered = 0.0;
    for (std::size_t i = 1; i < route.size(); ++i) {
        const double seg_len = distanceM(route[i - 1], route[i]);
        if (covered + seg_len >= target_m) {
            const double t = seg_len == 0.0 ? 0.0 : (target_m - covered) / seg_len;
            const auto [east, north] = localVectorM(route[i - 1], route[i]);
            return offsetPoint(route[i - 1], east * t, north * t);
        }
        covered += seg_len;
    }
    return route.back();
}

std::pair<double, double> unitDirection(const GeoPoint& from, const GeoPoint& to) {
    const auto [east, north] = localVectorM(from, to);
    const double len = std::hypot(east, north);
    if (len == 0.0) {
        throw std::invalid_argument("direction vector must not be zero length");
    }
    return {east / len, north / len};
}

std::pair<double, double> routeTangentAt(const Route& route, double target_m) {
    const double len = routeLengthM(route);
    const double a = std::max(0.0, target_m - 0.5);
    const double b = std::min(len, target_m + 0.5);
    return unitDirection(pointAtDistance(route, a), pointAtDistance(route, b));
}

bool pointInPolygon(const std::vector<std::pair<double, double>>& polygon,
                    double x,
                    double y) {
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const auto [xi, yi] = polygon[i];
        const auto [xj, yj] = polygon[j];
        const bool crosses = ((yi > y) != (yj > y)) &&
                             (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (crosses) {
            inside = !inside;
        }
    }
    return inside;
}

Route lineBetween(const GeoPoint& a, const GeoPoint& b, double interval_m) {
    const double len = distanceM(a, b);
    const std::size_t steps = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(len / std::max(1.0, interval_m))));
    const auto [east, north] = localVectorM(a, b);
    Route route;
    route.reserve(steps + 1);
    for (std::size_t i = 0; i <= steps; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(steps);
        route.push_back(offsetPoint(a, east * t, north * t));
    }
    return route;
}

Route sightLine(const SightSurfaceQuery& query) {
    const auto [east, north] = unitDirection(query.base, query.direction);
    const std::size_t steps = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(query.max_distance_m /
                                             std::max(1.0, query.sample_interval_m))));
    Route route;
    route.reserve(steps + 1);
    for (std::size_t i = 0; i <= steps; ++i) {
        const double d = query.max_distance_m * static_cast<double>(i) /
                         static_cast<double>(steps);
        route.push_back(offsetPoint(query.base, east * d, north * d));
    }
    return route;
}

} // namespace

TerrainAnalysis::TerrainAnalysis(const IDtmProvider& dtm, const IHeBackend& he)
    : dtm_(dtm), he_(he) {}

SlopeResult TerrainAnalysis::slope(const Route& route) const {
    const ElevationResult elev = elevation(route);
    SlopeResult result;
    result.distances_m = elev.distances_m;
    result.points = elev.points;
    result.slope_degrees.resize(elev.elevations_m.size(), 0.0);

    for (std::size_t i = 1; i < elev.elevations_m.size(); ++i) {
        const double rise = elev.elevations_m[i] - elev.elevations_m[i - 1];
        const double run = elev.distances_m[i] - elev.distances_m[i - 1];
        result.slope_degrees[i] = run == 0.0 ? 0.0 : std::atan(rise / run) * 180.0 / kPi;
    }
    return result;
}

SurfaceSlopeResult TerrainAnalysis::surfaceSlope(const Route& points, double radius_m) const {
    if (radius_m <= 0.0) {
        throw std::invalid_argument("surfaceSlope radius must be positive");
    }

    Route neighbor_points;
    neighbor_points.reserve(points.size() * 4);
    for (const GeoPoint& point : points) {
        neighbor_points.push_back(offsetPoint(point, radius_m, 0.0));
        neighbor_points.push_back(offsetPoint(point, -radius_m, 0.0));
        neighbor_points.push_back(offsetPoint(point, 0.0, radius_m));
        neighbor_points.push_back(offsetPoint(point, 0.0, -radius_m));
    }

    const ElevationResult neighbor_elevation = elevation(neighbor_points);
    SurfaceSlopeResult result;
    result.points = points;
    result.slope_degrees.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        const double east = neighbor_elevation.elevations_m[i * 4 + 0];
        const double west = neighbor_elevation.elevations_m[i * 4 + 1];
        const double north = neighbor_elevation.elevations_m[i * 4 + 2];
        const double south = neighbor_elevation.elevations_m[i * 4 + 3];
        const double dz_dx = (east - west) / (2.0 * radius_m);
        const double dz_dy = (north - south) / (2.0 * radius_m);
        result.slope_degrees.push_back(std::atan(std::hypot(dz_dx, dz_dy)) * 180.0 / kPi);
    }
    return result;
}

TerrainSectionResult TerrainAnalysis::terrainSection(const Route& route,
                                                     const TerrainSectionConfig& config) const {
    if (route.size() < 2) {
        throw std::invalid_argument("terrainSection route needs at least two points");
    }
    if (config.section_length_m <= 0.0 || config.interval_m <= 0.0) {
        throw std::invalid_argument("terrainSection length and interval must be positive");
    }

    TerrainSectionResult result;
    result.longitudinal = elevation(route);

    const double total_len = routeLengthM(route);
    for (double center_dist = config.interval_m; center_dist < total_len; center_dist += config.interval_m) {
        const GeoPoint center = pointAtDistance(route, center_dist);
        const auto [east, north] = routeTangentAt(route, center_dist);
        const double perp_east = -north;
        const double perp_north = east;
        const double half = config.section_length_m / 2.0;
        const GeoPoint left = offsetPoint(center, -perp_east * half, -perp_north * half);
        const GeoPoint right = offsetPoint(center, perp_east * half, perp_north * half);
        result.transverse.push_back(elevation(lineBetween(left, right, config.interval_m / 2.0)));
        result.transverse_center_distances_m.push_back(center_dist);
    }
    return result;
}

SightSurfaceEncryptedResult TerrainAnalysis::sightSurfaceEncrypted(const SightSurfaceQuery& query) const {
    if (!std::isfinite(query.max_distance_m) || query.max_distance_m <= 0.0) {
        throw std::invalid_argument("SightSurfaceQuery::max_distance_m must be finite and positive");
    }
    if (!std::isfinite(query.sample_interval_m) || query.sample_interval_m <= 0.0) {
        throw std::invalid_argument("SightSurfaceQuery::sample_interval_m must be finite and positive");
    }
    if (!std::isfinite(query.angle_degrees) || query.angle_degrees <= -90.0 ||
        query.angle_degrees >= 90.0) {
        throw std::invalid_argument("SightSurfaceQuery::angle_degrees must be strictly between -90 and 90");
    }
    if (!std::isfinite(query.base.lon) || !std::isfinite(query.base.lat) ||
        !std::isfinite(query.direction.lon) || !std::isfinite(query.direction.lat)) {
        throw std::invalid_argument("SightSurfaceQuery coordinates must be finite");
    }
    // sightLine uses this same clamped interval and adds the base endpoint.
    // Check the floating-point step count before its conversion to size_t.
    const double steps = std::ceil(query.max_distance_m / std::max(1.0, query.sample_interval_m));
    if (!std::isfinite(steps) || steps >= static_cast<double>(Route{}.max_size())) {
        throw std::length_error("SightSurfaceQuery requests too many samples");
    }
    const std::size_t slots = he_.slotCount();
    if (slots == 0) {
        throw std::invalid_argument("sightSurface requires nonzero slotCount");
    }

    SightSurfaceEncryptedResult result;
    result.points = sightLine(query);
    for (const auto& point : result.points) {
        if (!std::isfinite(point.lon) || !std::isfinite(point.lat)) {
            throw std::invalid_argument("SightSurfaceQuery produced non-finite sample coordinates");
        }
    }
    result.distances_m = cumulativeDistances(result.points);
    const double tan_theta = std::tan(toRad(query.angle_degrees));

    for (std::size_t begin = 0; begin < result.points.size();) {
        const std::size_t count = std::min(slots, result.points.size() - begin);
        const Route route_batch(result.points.begin() + begin,
                                result.points.begin() + begin + count);
        const Route base_batch(count, query.base);
        std::vector<double> negative_correction(count);
        for (std::size_t i = 0; i < count; ++i) {
            negative_correction[i] = -result.distances_m[begin + i] * tan_theta;
            if (!std::isfinite(negative_correction[i])) {
                throw std::invalid_argument("SightSurfaceQuery produced a non-finite plane correction");
            }
        }
        const CipherVector ground_ct = elevationEncrypted(route_batch);
        const CipherVector base_ct = elevationEncrypted(base_batch);
        // base + d*tan(theta), expressed using the backend's existing subPlain.
        const CipherVector plane_ct = he_.subPlain(base_ct, he_.encode(negative_correction));
        result.batches.push_back({begin, count, he_.sub(plane_ct, ground_ct)});
        begin += count;
    }
    return result;
}

SightSurfaceResult TerrainAnalysis::sightSurface(const SightSurfaceQuery& query) const {
    SightSurfaceEncryptedResult encrypted = sightSurfaceEncrypted(query);
    SightSurfaceResult result;
    result.points = std::move(encrypted.points);
    result.distances_m = std::move(encrypted.distances_m);
    result.buildable_height_m.resize(result.points.size());
    for (const auto& batch : encrypted.batches) {
        const PlainVector h = he_.decrypt(batch.height_difference);
        if (h.values.size() < batch.sample_count) {
            throw std::runtime_error("sightSurface decryption returned too few slots");
        }
        for (std::size_t i = 0; i < batch.sample_count; ++i) {
            result.buildable_height_m[batch.begin + i] = std::max(0.0, h.values[i]);
        }
    }
    return result;
}

EarthworkResult TerrainAnalysis::earthwork(const EarthworkQuery& query) const {
    if (query.polygon.size() < 3) {
        throw std::invalid_argument("earthwork polygon needs at least three points");
    }
    if (query.sample_interval_m <= 0.0) {
        throw std::invalid_argument("earthwork sample interval must be positive");
    }

    const GeoPoint origin = query.polygon.front();
    std::vector<std::pair<double, double>> polygon_m;
    polygon_m.reserve(query.polygon.size());
    double min_x = 0.0;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
    for (std::size_t i = 0; i < query.polygon.size(); ++i) {
        const auto xy = localVectorM(origin, query.polygon[i]);
        polygon_m.push_back(xy);
        if (i == 0) {
            min_x = max_x = xy.first;
            min_y = max_y = xy.second;
        } else {
            min_x = std::min(min_x, xy.first);
            max_x = std::max(max_x, xy.first);
            min_y = std::min(min_y, xy.second);
            max_y = std::max(max_y, xy.second);
        }
    }

    Route sample_points;
    std::vector<std::pair<double, double>> sample_offsets_m;
    for (double y = min_y + query.sample_interval_m / 2.0; y <= max_y;
         y += query.sample_interval_m) {
        for (double x = min_x + query.sample_interval_m / 2.0; x <= max_x;
             x += query.sample_interval_m) {
            if (pointInPolygon(polygon_m, x, y)) {
                sample_offsets_m.emplace_back(x, y);
                sample_points.push_back(offsetPoint(origin, x, y));
            }
        }
    }
    if (sample_points.empty()) {
        throw std::runtime_error("earthwork polygon produced no DTM samples");
    }

    const ElevationResult elevation_profile = elevation(sample_points);
    const std::vector<double>& ground = elevation_profile.elevations_m;

    EarthworkResult result;
    result.sample_count = ground.size();
    result.sample_area_m2 = query.sample_interval_m * query.sample_interval_m;
    result.total_area_m2 = result.sample_area_m2 * static_cast<double>(result.sample_count);
    result.min_ground_m = *std::min_element(ground.begin(), ground.end());
    result.max_ground_m = *std::max_element(ground.begin(), ground.end());
    double ground_sum = 0.0;
    for (double z : ground) {
        ground_sum += z;
    }
    result.average_ground_m = ground_sum / static_cast<double>(ground.size());

    std::vector<double> design_m(ground.size(), 0.0);
    for (std::size_t i = 0; i < ground.size(); ++i) {
        const auto [x, y] = sample_offsets_m[i];
        design_m[i] = query.design_height_m +
                      query.design_slope_east * x +
                      query.design_slope_north * y;
        const double height_delta_m = ground[i] - design_m[i];
        const double volume_m3 = height_delta_m * result.sample_area_m2;
        if (height_delta_m > 0.0) {
            result.cut_volume_m3 += volume_m3;
        } else {
            result.fill_volume_m3 += -volume_m3;
        }
    }
    result.net_volume_m3 = result.cut_volume_m3 - result.fill_volume_m3;

    // --- 순 토공량을 암호문 상태로 계산 (절토/성토 분리 없음) ---
    // 회로: 모서리 4점 암호화 -> 평문 가중치 곱(1레벨) -> 4항 덧셈 -> 설계고 뺄셈
    //       -> 셀 넓이 곱(1레벨, 패딩 slot은 0) -> 전체 slot 합(EvalSum) -> 복호화.
    // 비교가 없으므로 얕은 회로이고 부트스트랩이 필요 없다.
    if (ground.size() >= 1 && ground.size() <= he_.slotCount()) {
        const std::size_t n = ground.size();
        std::vector<double> z00(n), z10(n), z01(n), z11(n);
        std::vector<double> w00(n), w10(n), w01(n), w11(n);
        std::vector<double> area_mask(he_.slotCount(), 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const InterpolationStencil s = dtm_.interpolationStencil(sample_points[i]);
            z00[i] = s.z00_m; z10[i] = s.z10_m; z01[i] = s.z01_m; z11[i] = s.z11_m;
            w00[i] = (1.0 - s.dx) * (1.0 - s.dy);
            w10[i] = s.dx * (1.0 - s.dy);
            w01[i] = (1.0 - s.dx) * s.dy;
            w11[i] = s.dx * s.dy;
            area_mask[i] = result.sample_area_m2;
        }
        const CipherVector ground_ct = he_.add(
            he_.add(he_.mulPlain(he_.encrypt(he_.encode(z00)), he_.encode(w00)),
                    he_.mulPlain(he_.encrypt(he_.encode(z10)), he_.encode(w10))),
            he_.add(he_.mulPlain(he_.encrypt(he_.encode(z01)), he_.encode(w01)),
                    he_.mulPlain(he_.encrypt(he_.encode(z11)), he_.encode(w11))));
        const CipherVector delta_ct = he_.subPlain(ground_ct, he_.encode(design_m));
        const CipherVector volume_ct = he_.mulPlain(delta_ct, he_.encode(area_mask));
        const CipherVector net_ct = he_.sumSlots(volume_ct);
        result.net_volume_encrypted_m3 = he_.decrypt(net_ct).values.at(0);
        result.net_volume_abs_error_m3 =
            std::abs(result.net_volume_encrypted_m3 - result.net_volume_m3);
        result.encrypted_net_computed = true;
    }
    return result;
}

} // namespace openfhe_dtm

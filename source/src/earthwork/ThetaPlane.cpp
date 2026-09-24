#include "openfhe_dtm/EarthworkTheta.hpp"
#include "openfhe_dtm/Analysis.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace openfhe_dtm {
namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kEarthRadiusM = 6371008.8;

double toRad(double degrees) { return degrees * kPi / 180.0; }

// The same local coordinates and center-in-polygon quadrature as earthwork().
std::pair<double, double> localVectorM(const GeoPoint& from, const GeoPoint& to) {
    return {toRad(to.lon - from.lon) * kEarthRadiusM * std::cos(toRad(from.lat)),
            toRad(to.lat - from.lat) * kEarthRadiusM};
}

GeoPoint offsetPoint(const GeoPoint& origin, double east, double north) {
    return {origin.lon + east / (kEarthRadiusM * std::cos(toRad(origin.lat))) * 180.0 / kPi,
            origin.lat + north / kEarthRadiusM * 180.0 / kPi};
}

bool pointInPolygon(const std::vector<std::pair<double, double>>& polygon,
                    double x, double y) {
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const auto [xi, yi] = polygon[i];
        const auto [xj, yj] = polygon[j];
        if (((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) inside = !inside;
    }
    return inside;
}
} // namespace

EarthworkThetaEncryptedResult earthworkThetaEncrypted(
    const EarthworkThetaQuery& query, const IDtmProvider& dtm, const IHeBackend& he) {
    const std::size_t slots = he.slotCount();
    if (slots == 0 || query.polygon.size() < 3)
        throw std::invalid_argument("Theta earthwork requires nonzero slots and at least three polygon points");
    if (!std::isfinite(query.design_height_m) || !std::isfinite(query.slope_angle_degrees) ||
        std::abs(query.slope_angle_degrees) >= 90.0)
        throw std::invalid_argument("H must be finite and theta must be strictly between -90 and 90 degrees");
    const double step = query.sample_interval_m;
    const double area = step * step;
    if (!std::isfinite(step) || step <= 0.0 || !std::isfinite(area) || area <= 0.0)
        throw std::invalid_argument("Sample interval and its squared area must be finite and positive");
    for (const auto& p : query.polygon) {
        if (!std::isfinite(p.lon) || !std::isfinite(p.lat) || std::abs(p.lat) >= 90.0)
            throw std::invalid_argument("Polygon coordinates must be finite and away from the poles");
    }

    const GeoPoint origin = query.polygon.front();
    const auto [edge_east, edge_north] = localVectorM(origin, query.polygon[1]);
    const double edge_length = std::hypot(edge_east, edge_north);
    if (!std::isfinite(edge_length) || edge_length <= 0.0)
        throw std::invalid_argument("The first two polygon points must define a nonzero finite direction");
    const double direction_east = edge_east / edge_length;
    const double direction_north = edge_north / edge_length;
    const double slope = std::tan(toRad(query.slope_angle_degrees));
    if (!std::isfinite(slope)) throw std::invalid_argument("Non-finite tangent of theta");

    std::vector<std::pair<double, double>> polygon_m;
    double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0;
    for (const auto& p : query.polygon) {
        const auto xy = localVectorM(origin, p);
        if (!std::isfinite(xy.first) || !std::isfinite(xy.second))
            throw std::invalid_argument("Polygon local coordinates overflow");
        polygon_m.push_back(xy);
        min_x = std::min(min_x, xy.first); max_x = std::max(max_x, xy.first);
        min_y = std::min(min_y, xy.second); max_y = std::max(max_y, xy.second);
    }
    // Bound enumeration even for extremely sparse or tiny-interval polygons.
    const long double nx = std::ceil((static_cast<long double>(max_x) - min_x) / step);
    const long double ny = std::ceil((static_cast<long double>(max_y) - min_y) / step);
    if (nx * ny > 10000000.0L)
        throw std::length_error("Theta earthwork bounding grid exceeds 10000000 cells");

    Route points;
    std::vector<double> design;
    for (double y = min_y + step / 2.0; y <= max_y;) {
        for (double x = min_x + step / 2.0; x <= max_x;) {
            if (pointInPolygon(polygon_m, x, y)) {
                if (points.size() == slots)
                    throw std::length_error("Theta earthwork needs at most slotCount samples; increase the interval or slot count");
                // ax = tan(theta)*direction_east, ay = tan(theta)*direction_north.
                // d is a SIGNED projection, not the radial distance from origin.
                const double distance = x * direction_east + y * direction_north;
                const double height = query.design_height_m + distance * slope;
                if (!std::isfinite(height)) throw std::invalid_argument("Design height overflow");
                points.push_back(offsetPoint(origin, x, y));
                design.push_back(height);
            }
            const double next = x + step;
            if (!(next > x)) throw std::invalid_argument("Sample interval does not advance the x coordinate");
            x = next;
        }
        const double next = y + step;
        if (!(next > y)) throw std::invalid_argument("Sample interval does not advance the y coordinate");
        y = next;
    }
    if (points.empty()) throw std::runtime_error("Theta earthwork polygon produced no samples");

    EarthworkThetaEncryptedResult result;
    result.geometry.sample_count = points.size();
    result.geometry.sample_area_m2 = area;
    result.geometry.total_area_m2 = area * static_cast<double>(points.size());
    result.geometry.origin = origin;
    result.geometry.direction_east = direction_east;
    result.geometry.direction_north = direction_north;
    if (!std::isfinite(result.geometry.total_area_m2))
        throw std::invalid_argument("Total integration area overflow");

    // Fill minus cut directly: (public design - encrypted ground) * area.
    // Plaintext-ciphertext subtraction returns a ciphertext, without decryption.
    // Unused slots MUST retain zero weights before sumSlots.
    std::vector<double> area_mask(slots, 0.0);
    std::fill_n(area_mask.begin(), points.size(), area);
    const TerrainAnalysis analysis(dtm, he);
    const CipherVector ground = analysis.elevationEncrypted(points);
    const CipherVector delta = he.plainSub(he.encode(design), ground);
    const CipherVector volume = he.mulPlain(delta, he.encode(area_mask));
    result.net_volume = he.sumSlots(volume);
    return result;
}

EarthworkThetaResult finalizeEarthworkTheta(
    const EarthworkThetaEncryptedResult& encrypted, const IHeBackend& he) {
    if (encrypted.geometry.sample_count == 0 ||
        encrypted.geometry.sample_count > he.slotCount() ||
        !encrypted.net_volume.native ||
        encrypted.net_volume.slot_count != he.slotCount())
        throw std::invalid_argument("Invalid encrypted theta earthwork result");
    const PlainVector final_value = he.decrypt(encrypted.net_volume);
    if (final_value.values.empty() || !std::isfinite(final_value.values.front()))
        throw std::runtime_error("Final theta earthwork decryption did not produce a finite scalar");
    EarthworkThetaResult result = encrypted.geometry;
    result.net_volume_m3 = final_value.values.front();
    return result;
}

EarthworkThetaResult earthworkTheta(
    const EarthworkThetaQuery& query, const IDtmProvider& dtm, const IHeBackend& he) {
    return finalizeEarthworkTheta(earthworkThetaEncrypted(query, dtm, he), he);
}
} // namespace openfhe_dtm

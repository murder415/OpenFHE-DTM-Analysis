#include "openfhe_dtm/Analysis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace openfhe_dtm {
namespace {

std::vector<double> samplePlainElevation(const IDtmProvider& dtm, const Route& route) {
    std::vector<double> elevations;
    elevations.reserve(route.size());
    for (const auto& point : route) {
        elevations.push_back(dtm.elevationAt(point));
    }
    return elevations;
}

std::vector<double> cumulativeDistances(const Route& route) {
    constexpr double kPi = 3.14159265358979323846;
    constexpr double kEarthRadiusM = 6371008.8;
    const auto to_rad = [=](double degrees) { return degrees * kPi / 180.0; };
    std::vector<double> distances(route.size(), 0.0);
    for (std::size_t i = 1; i < route.size(); ++i) {
        const double lat1 = to_rad(route[i - 1].lat);
        const double lat2 = to_rad(route[i].lat);
        const double dlat = to_rad(route[i].lat - route[i - 1].lat);
        const double dlon = to_rad(route[i].lon - route[i - 1].lon);
        const double h = std::sin(dlat / 2.0) * std::sin(dlat / 2.0) +
                         std::cos(lat1) * std::cos(lat2) *
                             std::sin(dlon / 2.0) * std::sin(dlon / 2.0);
        distances[i] = distances[i - 1] +
                       2.0 * kEarthRadiusM *
                           std::asin(std::min(1.0, std::sqrt(h)));
    }
    return distances;
}

} // namespace

ElevationAccuracy compareElevation(const std::vector<double>& actual,
                                   const std::vector<double>& expected) {
    const std::size_t n = std::min(actual.size(), expected.size());
    ElevationAccuracy accuracy;
    accuracy.compared = n;
    if (n == 0) {
        return accuracy;
    }
    double abs_sum = 0.0;
    double sq_sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double error = actual[i] - expected[i];
        const double abs_error = std::abs(error);
        abs_sum += abs_error;
        sq_sum += error * error;
        accuracy.max_abs_error_m = std::max(accuracy.max_abs_error_m, abs_error);
    }
    accuracy.mae_m = abs_sum / static_cast<double>(n);
    accuracy.rmse_m = std::sqrt(sq_sum / static_cast<double>(n));
    return accuracy;
}

std::vector<double> TerrainAnalysis::plainElevation(const Route& route) const {
    return samplePlainElevation(dtm_, route);
}

CipherVector TerrainAnalysis::elevationEncrypted(const Route& route) const {
    const std::size_t slots = he_.slotCount();
    if (slots == 0 || route.empty() || route.size() > slots) {
        throw std::invalid_argument("elevationEncrypted requires 1..slotCount points and nonzero slotCount");
    }

    std::vector<double> z00, z10, z01, z11;
    std::vector<double> w00, w10, w01, w11;
    for (const auto& point : route) {
        const InterpolationStencil s = dtm_.interpolationStencil(point);
        if (!std::isfinite(s.z00_m) || !std::isfinite(s.z10_m) ||
            !std::isfinite(s.z01_m) || !std::isfinite(s.z11_m) ||
            !std::isfinite(s.dx) || !std::isfinite(s.dy) ||
            s.dx < 0.0 || s.dx > 1.0 || s.dy < 0.0 || s.dy > 1.0) {
            throw std::invalid_argument("interpolation stencil requires finite heights and fractions in [0,1]");
        }
        z00.push_back(s.z00_m);
        z10.push_back(s.z10_m);
        z01.push_back(s.z01_m);
        z11.push_back(s.z11_m);
        w00.push_back((1.0 - s.dx) * (1.0 - s.dy));
        w10.push_back(s.dx * (1.0 - s.dy));
        w01.push_back((1.0 - s.dx) * s.dy);
        w11.push_back(s.dx * s.dy);
    }

    // The provider supplies corner samples at the existing input boundary.
    // Encoding pads each array to the backend's slot count; it is not encryption.
    const auto encodeWeights = [&](const std::vector<double>& weights) {
        // A grid-corner round trip can leave an entire coefficient vector at
        // e.g. 1e-29. OpenFHE rejects this below-resolution nonzero plaintext,
        // whereas an exact zero vector is valid. Discard only a whole vector
        // below double epsilon; the per-height perturbation is bounded by
        // epsilon times the corresponding maximum corner height.
        if (std::all_of(weights.begin(), weights.end(), [](double w) {
                return std::abs(w) < std::numeric_limits<double>::epsilon();
            })) return he_.encode(std::vector<double>(weights.size(), 0.0));
        return he_.encode(weights);
    };
    const CipherVector c00 = he_.mulPlain(he_.encrypt(he_.encode(z00)), encodeWeights(w00));
    const CipherVector c10 = he_.mulPlain(he_.encrypt(he_.encode(z10)), encodeWeights(w10));
    const CipherVector c01 = he_.mulPlain(he_.encrypt(he_.encode(z01)), encodeWeights(w01));
    const CipherVector c11 = he_.mulPlain(he_.encrypt(he_.encode(z11)), encodeWeights(w11));
    return he_.add(he_.add(c00, c10), he_.add(c01, c11));
}

ElevationResult TerrainAnalysis::elevation(const Route& route) const {
    ElevationResult result;
    result.points = route;
    result.distances_m = cumulativeDistances(route);
    result.elevations_m.reserve(route.size());
    const std::size_t slots = he_.slotCount();
    if (slots == 0) {
        throw std::invalid_argument("elevation requires nonzero slotCount");
    }
    for (std::size_t begin = 0; begin < route.size();) {
        const std::size_t count = std::min(slots, route.size() - begin);
        const Route batch(route.begin() + begin, route.begin() + begin + count);
        const PlainVector decrypted = he_.decrypt(elevationEncrypted(batch));
        if (decrypted.values.size() < count) {
            throw std::runtime_error("elevation decryption returned too few slots");
        }
        result.elevations_m.insert(result.elevations_m.end(), decrypted.values.begin(),
                                   decrypted.values.begin() + count);
        begin += count;
    }
    return result;
}

} // namespace openfhe_dtm

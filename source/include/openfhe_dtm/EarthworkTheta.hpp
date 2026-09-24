#pragma once

#include "openfhe_dtm/DtmProvider.hpp"
#include "openfhe_dtm/HeBackend.hpp"

#include <cstddef>
#include <limits>

namespace openfhe_dtm {

// Planar net-volume variant. Positive theta rises from polygon[0] toward
// polygon[1]. H is the plane height at polygon[0]; theta is in degrees.
// This changes the planar ax/ay input convention, not the integration rule.
struct EarthworkThetaQuery {
    Route polygon;
    double design_height_m = 0.0;
    double slope_angle_degrees = 0.0; // finite, -90 < theta < 90
    double sample_interval_m = 5.0;
};

struct EarthworkThetaResult {
    std::size_t sample_count = 0;
    double sample_area_m2 = 0.0; // sample_interval_m squared
    double total_area_m2 = 0.0;
    GeoPoint origin;
    double direction_east = 0.0;  // first edge, normalized in the local plane
    double direction_north = 0.0;
    double net_volume_m3 = std::numeric_limits<double>::quiet_NaN();
    // FILL minus CUT, sum((design - ground) * area).
};

struct EarthworkThetaEncryptedResult {
    EarthworkThetaResult geometry; // net_volume_m3 is not populated yet
    CipherVector net_volume;       // slot sum, not individual ground heights
};

// One batch, 1 <= sample_count <= he.slotCount(), as in the original circuit.
// Uses public H/theta and public geometry, encrypted interpolation, plainSub,
// mulPlain and sumSlots: (public design - encrypted ground) * positive area.
// Unused slots have zero area; the encrypted sum is fill minus cut.
// No decrypt or plain elevationAt call is made here.
EarthworkThetaEncryptedResult earthworkThetaEncrypted(
    const EarthworkThetaQuery& query, const IDtmProvider& dtm, const IHeBackend& he);

// Exactly one final scalar decryption; no encrypted work follows it.
EarthworkThetaResult finalizeEarthworkTheta(
    const EarthworkThetaEncryptedResult& encrypted, const IHeBackend& he);

EarthworkThetaResult earthworkTheta(
    const EarthworkThetaQuery& query, const IDtmProvider& dtm, const IHeBackend& he);

} // namespace openfhe_dtm

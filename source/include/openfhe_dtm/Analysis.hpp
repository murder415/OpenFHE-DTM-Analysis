#pragma once

#include "openfhe_dtm/DtmProvider.hpp"
#include "openfhe_dtm/HeBackend.hpp"

#include <cstddef>
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace openfhe_dtm {

struct ElevationResult {
    std::vector<double> distances_m;
    std::vector<double> elevations_m;
    std::vector<GeoPoint> points;
};

struct ElevationAccuracy {
    double mae_m = 0.0;
    double rmse_m = 0.0;
    double max_abs_error_m = 0.0;
    std::size_t compared = 0;
};

ElevationAccuracy compareElevation(const std::vector<double>& actual,
                                   const std::vector<double>& expected);

struct SlopeResult {
    std::vector<double> distances_m;
    std::vector<double> slope_degrees;
    std::vector<GeoPoint> points;
};

struct SurfaceSlopeResult {
    std::vector<GeoPoint> points;
    std::vector<double> slope_degrees;
};

struct TerrainSectionConfig {
    double section_length_m = 40.0;
    double interval_m = 20.0;
};

struct TerrainSectionResult {
    ElevationResult longitudinal;
    std::vector<ElevationResult> transverse;
    std::vector<double> transverse_center_distances_m;
};

struct SightSurfaceQuery {
    GeoPoint base;
    GeoPoint direction;
    double angle_degrees = 15.0;
    double max_distance_m = 90.0;
    double sample_interval_m = 10.0;
};

struct SightSurfaceResult {
    std::vector<double> distances_m;
    // Intentionally empty: the final-result-only path does not decrypt the
    // intermediate ground elevations or regulation plane.
    std::vector<double> ground_m;
    std::vector<double> regulation_plane_m;
    std::vector<double> buildable_height_m;
    std::vector<GeoPoint> points;
};

struct SightSurfaceCipherBatch {
    std::size_t begin = 0;
    std::size_t sample_count = 0;
    // Signed regulation-plane minus ground height, before final max(0, h).
    CipherVector height_difference;
};

struct SightSurfaceEncryptedResult {
    std::vector<double> distances_m;
    std::vector<GeoPoint> points;
    std::vector<SightSurfaceCipherBatch> batches;
};

struct EarthworkQuery {
    Route polygon;
    double design_height_m = 0.0;
    double design_slope_east = 0.0;
    double design_slope_north = 0.0;
    double sample_interval_m = 5.0;
};

struct EarthworkResult {
    std::size_t sample_count = 0;
    double sample_area_m2 = 0.0;
    double total_area_m2 = 0.0;
    double cut_volume_m3 = 0.0;
    double fill_volume_m3 = 0.0;
    double net_volume_m3 = 0.0;
    double min_ground_m = 0.0;
    double max_ground_m = 0.0;
    double average_ground_m = 0.0;
    // 순 토공량(절토 - 성토)을 암호문 상태로 계산한 값. 셀별 (지반고 - 설계고)를
    // 암호문에서 구해 셀 넓이와 곱하고 전체 slot을 합산한다. 절토/성토 분리는
    // 부호 판정(비교)이 필요해 여기서는 다루지 않는다.
    bool encrypted_net_computed = false;
    double net_volume_encrypted_m3 = 0.0;
    double net_volume_abs_error_m3 = 0.0;  // |encrypted - plaintext|
};

// Boundary-envelope earthwork, separate from the preserved planar earthwork()
// baseline above. Coordinates, the polygon ring, grid and angle are public.
struct EarthworkSlopeQuery {
    Route polygon;
    // Both are explicit caller inputs. NaN defaults reject omitted values.
    double design_height_m = std::numeric_limits<double>::quiet_NaN();
    double slope_angle_degrees = std::numeric_limits<double>::quiet_NaN();
    double sample_interval_m = 5.0;
    // Zero selects sample_interval_m; otherwise a finite positive value.
    double boundary_interval_m = 0.0;
    // Height thresholds in metres: area classification and envelope diagnostics.
    // These defaults are numerical choices, not terrain-accuracy guarantees.
    double area_threshold_m = 1e-6;
    double envelope_tolerance_m = 1e-6;
    std::size_t max_grid_cells = 1000000;
    std::size_t max_boundary_points = 4096;
    std::size_t max_candidate_ciphertexts = 100000;
};

struct EarthworkSlopeGeometry {
    GeoPoint origin;
    double grid_min_east_m = 0.0;
    double grid_min_north_m = 0.0;
    double grid_dx_m = 0.0;
    double grid_dy_m = 0.0;
    std::size_t grid_columns = 0;
    std::size_t grid_rows = 0;
    // Center-inside cells contribute area/4 at each vertex; shared vertices
    // accumulate weights. Boundary cells are not clipped to the polygon.
    Route points;
    std::vector<double> point_weights_m2;
    std::vector<std::array<std::size_t, 4>> cells;
    Route boundary_points;
    double integration_area_m2 = 0.0;
};

struct EarthworkSlopeCipherBatch {
    std::size_t begin = 0;
    std::size_t sample_count = 0;
    // Target order: geometry.points, followed by geometry.boundary_points.
    // Boundary-only diagnostic targets use grid_dx_m*grid_dy_m/4 as a positive
    // numerical scale (removed at final comparison), and are never
    // included in volume/area sums. Padding is zero and ignored on decryption.
    CipherVector target_minus_ground;
    std::vector<CipherVector> upper_candidates;
    std::vector<CipherVector> lower_candidates;
};

struct EarthworkSlopeEncryptedResult {
    EarthworkSlopeGeometry geometry;
    double area_threshold_m = 1e-6;
    double envelope_tolerance_m = 1e-6;
    std::vector<EarthworkSlopeCipherBatch> batches;
};

struct EarthworkSlopeResult {
    EarthworkSlopeGeometry geometry;
    bool envelope_valid = false;
    std::string status;
    std::size_t infeasible_point_count = 0;
    std::size_t infeasible_boundary_point_count = 0;
    std::size_t infeasible_cell_count = 0;
    double max_envelope_gap_m = 0.0;
    // NaN when the boundary envelope is inconsistent beyond its tolerance.
    double cut_volume_m3 = std::numeric_limits<double>::quiet_NaN();
    double fill_volume_m3 = std::numeric_limits<double>::quiet_NaN();
    double net_volume_m3 = std::numeric_limits<double>::quiet_NaN(); // fill - cut
    double absolute_moved_volume_m3 = std::numeric_limits<double>::quiet_NaN();
    double cut_area_m2 = std::numeric_limits<double>::quiet_NaN();
    double fill_area_m2 = std::numeric_limits<double>::quiet_NaN();
    double unchanged_area_m2 = std::numeric_limits<double>::quiet_NaN();
    // Owner-side diagnostics only; derived per-point values are not safe to
    // release merely because the original elevation was not decrypted directly.
    std::vector<double> point_height_changes_m;
};

class TerrainAnalysis {
public:
    TerrainAnalysis(const IDtmProvider& dtm, const IHeBackend& he);

    std::vector<double> plainElevation(const Route& route) const;
    // InterpEnc: one nonempty batch of at most he.slotCount() points.
    // Four encrypted corner arrays are weighted by plaintext bilinear weights;
    // the result stays encrypted. No decryption is performed here.
    CipherVector elevationEncrypted(const Route& route) const;
    // Standalone elevation query: decrypts the interpolated result once per batch.
    ElevationResult elevation(const Route& route) const;
    SlopeResult slope(const Route& route) const;
    SurfaceSlopeResult surfaceSlope(const Route& points, double radius_m) const;
    TerrainSectionResult terrainSection(const Route& route, const TerrainSectionConfig& config) const;
    // All interpolation, plane correction and subtraction stay encrypted.
    SightSurfaceEncryptedResult sightSurfaceEncrypted(const SightSurfaceQuery& query) const;
    // Decrypts only the final signed difference once per batch, then applies max(0,h).
    // ground_m and regulation_plane_m are intentionally left empty.
    SightSurfaceResult sightSurface(const SightSurfaceQuery& query) const;
    EarthworkResult earthwork(const EarthworkQuery& query) const;
    // Encrypt H, interpolate terrain/boundary heights, and construct all upper,
    // lower and target differences times quadrature weights without decryption.
    // The current experimental backend also owns a secret key; this method
    // itself needs no decrypt operation and does not claim key-role separation.
    EarthworkSlopeEncryptedResult earthworkSlopeEncrypted(const EarthworkSlopeQuery& query) const;
    // Owner-side final stage: decrypt candidate vectors, take min/max, diagnose
    // incompatible envelopes, and sum cut/fill. It does not decrypt raw heights,
    // but candidates may reveal heights by inversion when H/weights are known.
    EarthworkSlopeResult finalizeEarthworkSlope(const EarthworkSlopeEncryptedResult& encrypted) const;
    EarthworkSlopeResult earthworkSlope(const EarthworkSlopeQuery& query) const;

private:
    const IDtmProvider& dtm_;
    const IHeBackend& he_;
};

} // namespace openfhe_dtm

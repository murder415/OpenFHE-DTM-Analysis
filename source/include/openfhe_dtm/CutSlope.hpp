#pragma once
// Earthwork with a planned height inside the polygon (spot-height / 점고법 grid)
// and a cut slope outside the polygon at a user excavation angle.
// Written from the public V-World 3D analysis manual description (not the
// boundary-cone min/max method): each outside vertex uses ONE public distance to
// the polygon boundary, s = H + d*tan(theta); only terrain above s is cut.
#include "openfhe_dtm/PackedElevation.hpp"
#include <array>
#include <vector>

namespace cutslope {
using namespace openfhe_dtm;

struct CutSlopeQuery {
    Route polygon;
    double design_height_m = 0;     // H, public
    double cut_angle_degrees = 45;  // theta, 0 < theta < 90
    double grid_interval_m = 5;     // requested cell size (upper bound)
    double slope_width_m = 10;      // public outer ring width examined for the cut slope
};

struct CutSlopeGeometry {
    Route points;                     // grid vertices with any weight
    std::vector<double> inside_w;     // A_c/4 per selected inside cell sharing the vertex
    std::vector<double> outside_w;    // same for outside-ring cells
    std::vector<double> slope_height; // H + d*tan(theta); d = 0 for vertices inside
    std::vector<double> boundary_distance;
    std::vector<bool> outer_edge;     // vertex lies on the outer edge of the examined ring
    double dx = 0, dy = 0;
    std::size_t inside_cells = 0, outside_cells = 0;
    double inside_area = 0, outside_area = 0;
};

struct CutSlopeResult {
    CutSlopeGeometry geometry;
    double inside_fill_minus_cut_m3 = 0;  // sum (H - z) * A_in, one scalar decryption
    double slope_cut_m3 = 0;              // sum max(0, z - s) * A_out, one vector decryption
    double net_fill_minus_cut_m3 = 0;     // inside - slope cut
    std::size_t decryptions = 0, batches = 0;
    std::size_t unfinished_slope_points = 0; // outer-edge vertices where terrain is still above s
};

CutSlopeGeometry buildGeometry(const CutSlopeQuery& q);
// Ground heights come from rotation interpolation on the encrypted packed DEM,
// so the provider never encrypts per-point corner heights.
CutSlopeResult cutSlopeEarthworkPacked(const CutSlopeQuery& q, PackedElevation& elevation);
} // namespace cutslope

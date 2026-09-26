// Rotation interpolation, sight surface and cut-slope earthwork on the encrypted
// packed DEM (32,768 slots), checked against independent plaintext results on a
// synthetic tilted plane, plus rejection of invalid inputs.
#include "openfhe_dtm/CutSlope.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace openfhe_dtm;
using namespace cutslope;

static double bilinear(const IDtmProvider& dtm, const GeoPoint& p) {
    const auto s = dtm.interpolationStencil(p);
    return s.z00_m * (1 - s.dx) * (1 - s.dy) + s.z10_m * s.dx * (1 - s.dy) +
           s.z01_m * (1 - s.dx) * s.dy + s.z11_m * s.dx * s.dy;
}

int main() {
    try {
        RasterGrid grid; grid.origin = {0, 0}; grid.width = 128; grid.height = 128; grid.cell_size_m = 5;
        for (std::size_t y = 0; y < grid.height; ++y)
            for (std::size_t x = 0; x < grid.width; ++x) grid.elevations_m.push_back(0.2 * 5.0 * x + 0.1 * 5.0 * y);
        DtmMetadata m; m.name = "tilted plane"; m.source = "synthetic"; m.crs = "local"; m.resolution_m = 5;
        RasterDtmProvider plane(m, std::move(grid));
        OpenFheBackend he(32768, {}, false);
        PackedElevation elevation(plane, he);

        // 1) interpolation at off-grid points
        Route pts;
        for (int k = 0; k < 20; ++k) pts.push_back(plane.pointFromMeters(13.7 + 29.3 * k, 21.1 + 17.9 * k));
        elevation.encryptFor(pts);
        std::size_t dec = 0;
        const auto z = elevation.decrypt(elevation.mergeDisjoint(elevation.interpolate(pts)), &dec);
        double zerr = 0;
        for (std::size_t i = 0; i < pts.size(); ++i) zerr = std::max(zerr, std::fabs(z[i] - bilinear(plane, pts[i])));
        std::cout << "interpolation max_abs=" << zerr << " decryptions=" << dec << '\n';
        if (!(zerr < 1e-4)) throw std::runtime_error("interpolation differs from plaintext");

        // 2) sight surface
        SightSurfaceQuery sq{plane.pointFromMeters(40, 40), plane.pointFromMeters(400, 300), 12.0, 300.0, 30.0};
        const auto sr = sightSurfacePacked(sq, elevation);
        std::cout << "sight points=" << sr.buildable_height_m.size() << " decryptions=" << sr.decryptions << '\n';
        if (sr.buildable_height_m.empty() || sr.decryptions < 1) throw std::runtime_error("sight surface failed");
        for (double h : sr.buildable_height_m)
            if (!(h >= 0) || !std::isfinite(h)) throw std::runtime_error("sight surface value invalid");

        // 3) cut-slope earthwork
        CutSlopeQuery q;
        q.polygon = {plane.pointFromMeters(150, 150), plane.pointFromMeters(350, 150),
                     plane.pointFromMeters(350, 350), plane.pointFromMeters(150, 350)};
        q.design_height_m = 40; q.cut_angle_degrees = 30; q.grid_interval_m = 10; q.slope_width_m = 80;
        const auto r = cutSlopeEarthworkPacked(q, elevation);
        long double in = 0, cut = 0;
        const auto& g = r.geometry;
        for (std::size_t i = 0; i < g.points.size(); ++i) {
            const long double zi = bilinear(plane, g.points[i]);
            in += (q.design_height_m - zi) * g.inside_w[i];
            const long double v = (zi - g.slope_height[i]) * g.outside_w[i];
            if (v > 0) cut += v;
        }
        const double err = std::fabs(r.net_fill_minus_cut_m3 - static_cast<double>(in - cut));
        std::cout << "earthwork vertices=" << g.points.size() << " decryptions=" << r.decryptions
                  << " net=" << r.net_fill_minus_cut_m3 << " abs_err=" << err << '\n';
        if (!(err < 1e-2)) throw std::runtime_error("net earthwork differs from plaintext");
        if (r.decryptions < 2) throw std::runtime_error("unexpected decryption count");

        // 4) invalid inputs
        int rejected = 0;
        for (double bad : {0.0, 90.0, -5.0}) {
            auto b = q; b.cut_angle_degrees = bad;
            try { buildGeometry(b); } catch (const std::invalid_argument&) { ++rejected; }
        }
        auto b = q; b.polygon.resize(2);
        try { buildGeometry(b); } catch (const std::invalid_argument&) { ++rejected; }
        try { elevation.locate(plane.pointFromMeters(-100, -100)); } catch (const std::out_of_range&) { ++rejected; }
        if (rejected != 5) throw std::runtime_error("invalid input accepted");
        std::cout << "PASS\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; }
}

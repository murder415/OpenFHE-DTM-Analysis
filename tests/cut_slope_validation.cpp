// Cut-slope earthwork check on a synthetic tilted plane against an independent
// plaintext sum over the same vertices and weights (1,024 slots -> 2 batches).
#include "openfhe_dtm/CutSlope.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace openfhe_dtm;
using namespace cutslope;

int main() {
    try {
        RasterGrid grid; grid.origin = {0, 0}; grid.width = 128; grid.height = 128; grid.cell_size_m = 5;
        for (std::size_t y = 0; y < grid.height; ++y)
            for (std::size_t x = 0; x < grid.width; ++x) grid.elevations_m.push_back(0.2 * 5.0 * x);
        DtmMetadata m; m.name = "tilted plane"; m.source = "synthetic"; m.crs = "local"; m.resolution_m = 5;
        RasterDtmProvider plane(m, std::move(grid));
        OpenFheBackend he(1024, {}, false);
        TerrainAnalysis analysis(plane, he);
        CutSlopeQuery q;
        q.polygon = {plane.pointFromMeters(150, 150), plane.pointFromMeters(350, 150),
                     plane.pointFromMeters(350, 350), plane.pointFromMeters(150, 350)};
        q.design_height_m = 40; q.cut_angle_degrees = 30; q.grid_interval_m = 10; q.slope_width_m = 80;
        const auto r = cutSlopeEarthwork(q, analysis, he);
        long double in = 0, cut = 0;
        const auto& g = r.geometry;
        for (std::size_t i = 0; i < g.points.size(); ++i) {
            const auto s = plane.interpolationStencil(g.points[i]);
            const long double dx = s.dx, dy = s.dy;
            const long double z = s.z00_m * (1 - dx) * (1 - dy) + s.z10_m * dx * (1 - dy) + s.z01_m * (1 - dx) * dy + s.z11_m * dx * dy;
            in += (q.design_height_m - z) * g.inside_w[i];
            const long double v = (z - g.slope_height[i]) * g.outside_w[i];
            if (v > 0) cut += v;
        }
        const double err = std::fabs(r.net_fill_minus_cut_m3 - static_cast<double>(in - cut));
        std::cout << "vertices=" << g.points.size() << " batches=" << r.batches << " decryptions=" << r.decryptions
                  << " net=" << r.net_fill_minus_cut_m3 << " abs_err=" << err << '\n';
        if (!(err < 1e-3)) throw std::runtime_error("net earthwork differs from plaintext");
        if (r.batches < 2 || r.decryptions != 2 * r.batches) throw std::runtime_error("unexpected batch/decryption count");
        int rejected = 0;
        for (double bad : {0.0, 90.0, -5.0}) {
            auto b = q; b.cut_angle_degrees = bad;
            try { buildGeometry(b); } catch (const std::invalid_argument&) { ++rejected; }
        }
        auto b = q; b.polygon.resize(2);
        try { buildGeometry(b); } catch (const std::invalid_argument&) { ++rejected; }
        if (rejected != 4) throw std::runtime_error("invalid input accepted");
        std::cout << "PASS\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; }
}

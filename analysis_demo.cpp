// Encrypted terrain analysis on a validation TDB (32,768 slots):
// rotation-based elevation interpolation, sight surface and cut-slope earthwork.
// usage: analysis_demo <input.tdb> [H_m] [cut_angle_deg]
#include "openfhe_dtm/CutSlope.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/TdbDataset.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

using namespace openfhe_dtm;
using namespace cutslope;
using Clock = std::chrono::steady_clock;
static double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

int main(int argc, char** argv) {
    if (argc < 2 || argc > 4) {
        std::cerr << "usage: analysis_demo <input.tdb> [H_m] [cut_angle_deg]\n";
        return 2;
    }
    try {
        const auto dem = loadValidationTdbDataset(argv[1]);
        const double H = argc > 2 ? std::stod(argv[2]) : 14.76;
        const double phi = argc > 3 ? std::stod(argv[3]) : 45.0;
        const double span = std::max(2.0, std::min(dem.widthMeters(), dem.heightMeters()) * 0.55);
        const double x0 = span * 0.12, y0 = span * 0.12;
        auto P = [&](double x, double y) { return dem.pointFromMeters(x, y); };

        auto t = Clock::now();
        OpenFheBackend he(32768, {}, false);
        std::cout << std::fixed << std::setprecision(3) << "keys: " << since(t) << " s\n";
        PackedElevation elevation(dem, he);

        t = Clock::now();
        const auto sight = sightSurfacePacked({P(x0, y0), P(x0 + span, y0 + span * 0.61), 12.0, span,
                                               std::max(2.0, span / 9.0)}, elevation);
        std::cout << "sight surface: " << sight.buildable_height_m.size() << " points, "
                  << sight.decryptions << " decryption(s), " << since(t) << " s\n";
        for (std::size_t i = 0; i < sight.buildable_height_m.size(); ++i)
            std::cout << "  d=" << sight.distances_m[i] << " m  buildable=" << sight.buildable_height_m[i] << " m\n";

        CutSlopeQuery q;
        q.polygon = {P(x0 + span * 0.11, y0 + span * 0.11), P(x0 + span * 0.78, y0 + span * 0.11),
                     P(x0 + span * 0.78, y0 + span * 0.78), P(x0 + span * 0.11, y0 + span * 0.78)};
        q.design_height_m = H; q.cut_angle_degrees = phi;
        q.grid_interval_m = std::max(2.0, span / 18.0); q.slope_width_m = 2 * q.grid_interval_m;
        t = Clock::now();
        const auto r = cutSlopeEarthworkPacked(q, elevation);
        std::cout << std::setprecision(6) << "cut-slope earthwork (H=" << H << " m, angle=" << phi << " deg): "
                  << r.geometry.points.size() << " vertices, " << r.decryptions << " decryption(s), "
                  << std::setprecision(3) << since(t) << " s\n" << std::scientific << std::setprecision(6)
                  << "  inside (fill - cut) = " << r.inside_fill_minus_cut_m3 << " m3\n"
                  << "  slope cut           = " << r.slope_cut_m3 << " m3\n"
                  << "  net (fill - cut)    = " << r.net_fill_minus_cut_m3 << " m3\n";
        const auto& st = elevation.stats();
        std::cout << "packed DEM ciphertexts: " << st.packs << ", rotations: " << st.rotations << '\n';
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}

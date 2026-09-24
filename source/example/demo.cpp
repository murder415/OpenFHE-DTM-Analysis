#include "openfhe_dtm/Analysis.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/RunOptions.hpp"
#include "openfhe_dtm/TdbDataset.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>

namespace {

int runDemo(const openfhe_dtm::IDtmProvider& dtm,
            const std::function<openfhe_dtm::GeoPoint(double, double)>& pointFromMeters,
            double analysis_span_m = 90.0) {
    using namespace openfhe_dtm;

    OpenFheBackend he(16);
    TerrainAnalysis analysis(dtm, he);

    std::cout << "OpenFHE DEM analysis validation\n";
    std::cout << "Validation DEM source: " << dtm.metadata().name << "\n";
    std::cout << "This demo runs DEM samples through real OpenFHE CKKS.\n";

    const auto a = he.encrypt(he.encode({1.0, 2.0, 3.0, 4.0}));
    const auto b = he.encrypt(he.encode({10.0, 20.0, 30.0, 40.0}));
    const auto added = he.decrypt(he.add(a, b));
    const auto rotated = he.decrypt(he.rotate(a, 1));

    std::cout << "OpenFHE add smoke test: "
              << added.values[0] << ", " << added.values[1] << ", "
              << added.values[2] << ", " << added.values[3] << "\n";
    std::cout << "OpenFHE rotate smoke test: "
              << rotated.values[0] << ", " << rotated.values[1] << ", "
              << rotated.values[2] << ", " << rotated.values[3] << "\n";

    const double span = std::max(2.0, analysis_span_m);
    const double x0 = span * 0.12;
    const double y0 = span * 0.12;
    const GeoPoint p0 = pointFromMeters(x0, y0);
    const Route route = {
        p0,
        pointFromMeters(x0 + span * 0.22, y0 + span * 0.11),
        pointFromMeters(x0 + span * 0.50, y0 + span * 0.39),
        pointFromMeters(x0 + span * 0.83, y0 + span * 0.61),
    };

    const auto elevation = analysis.elevation(route);
    const auto elevation_accuracy =
        compareElevation(elevation.elevations_m, analysis.plainElevation(route));
    const auto slope = analysis.slope(route);
    const auto surface_slope = analysis.surfaceSlope(Route{route[1], route[2]},
                                                     std::max(1.0, span / 30.0));
    const double section_interval_m = std::max(1.0, span / 5.0);
    const double section_length_m = std::max(2.0, span / 3.0);
    const auto section =
        analysis.terrainSection(route, TerrainSectionConfig{section_length_m, section_interval_m});
    const auto sight = analysis.sightSurface(
        SightSurfaceQuery{pointFromMeters(x0, y0), pointFromMeters(x0 + span, y0 + span * 0.61),
                          12.0, span, std::max(2.0, span / 9.0)});
    const Route earthwork_area = {
        pointFromMeters(x0 + span * 0.11, y0 + span * 0.11),
        pointFromMeters(x0 + span * 0.78, y0 + span * 0.11),
        pointFromMeters(x0 + span * 0.78, y0 + span * 0.78),
        pointFromMeters(x0 + span * 0.11, y0 + span * 0.78),
    };
    EarthworkQuery earthwork_query;
    earthwork_query.polygon = earthwork_area;
    earthwork_query.design_height_m =
        std::accumulate(elevation.elevations_m.begin(), elevation.elevations_m.end(), 0.0) /
        static_cast<double>(elevation.elevations_m.size());
    earthwork_query.sample_interval_m = std::max(2.0, span / 18.0);
    const auto earthwork = analysis.earthwork(earthwork_query);

    std::cout << "Elevation samples: " << elevation.elevations_m.size()
              << " first=" << elevation.elevations_m.front()
              << " last=" << elevation.elevations_m.back() << "\n";
    std::cout << "Elevation error: mae=" << elevation_accuracy.mae_m
              << " rmse=" << elevation_accuracy.rmse_m
              << " max=" << elevation_accuracy.max_abs_error_m << " m\n";
    std::cout << "Slope samples: " << slope.slope_degrees.size()
              << " last=" << slope.slope_degrees.back() << " deg\n";
    std::cout << "Surface slope samples: " << surface_slope.slope_degrees.size()
              << " first=" << surface_slope.slope_degrees.front() << " deg\n";
    std::cout << "Terrain transverse sections: " << section.transverse.size() << "\n";
    std::cout << "Terrain longitudinal elevations: first="
              << section.longitudinal.elevations_m.front() << " last="
              << section.longitudinal.elevations_m.back() << " m\n";
    if (!section.transverse.empty()) {
        const auto& profile = section.transverse.front().elevations_m;
        const auto range = std::minmax_element(profile.begin(), profile.end());
        std::cout << "First transverse profile: samples=" << profile.size()
                  << " min=" << *range.first << " max=" << *range.second << " m\n";
    }
    std::cout << "Sight-surface samples: " << sight.buildable_height_m.size()
              << " (only final height differences are decrypted)"
              << " buildable[first/last]="
              << sight.buildable_height_m.front() << "/"
              << sight.buildable_height_m.back() << " m\n";
    std::cout << "Earthwork samples: " << earthwork.sample_count
              << " design=" << earthwork_query.design_height_m
              << " ground[min/avg/max]=" << earthwork.min_ground_m << "/"
              << earthwork.average_ground_m << "/" << earthwork.max_ground_m
              << " cut=" << earthwork.cut_volume_m3
              << " fill=" << earthwork.fill_volume_m3
              << " net=" << earthwork.net_volume_m3 << " m3\n";
    if (earthwork.encrypted_net_computed) {
        std::cout << "Earthwork net (ciphertext): "
                  << earthwork.net_volume_encrypted_m3 << " m3 vs plaintext "
                  << earthwork.net_volume_m3 << " m3  abs_error="
                  << earthwork.net_volume_abs_error_m3 << " m3\n";
    }

    return 0;
}

} // namespace

int main(int argc, char** argv) {
    using namespace openfhe_dtm;

    RunOptions options;
    try {
        options = parseRunOptions(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "invalid run options: " << e.what() << "\n";
        return 2;
    }

    try {
        if (originalDataRequested(options)) {
            const TdbDatasetInfo info = inspectTdbDataset(options.tdb_path);
            std::cout << "DEM execution requested: " << options.tdb_path << "\n"
                      << formatTdbDatasetInfo(info);
            RasterDtmProvider tdb = loadValidationTdbDataset(options.tdb_path);
            const double usable_span =
                std::min(tdb.widthMeters(), tdb.heightMeters()) * 0.55;
            return runDemo(tdb,
                           [&](double east_m, double north_m) {
                               return tdb.pointFromMeters(east_m, north_m);
                           },
                           usable_span);
        }

        if (!options.synthetic) {
            std::cerr << "Run with --synthetic for validation, or --tdb <elevation_bil32.tdb>.\n";
            return 2;
        }

        SyntheticRasterProvider synthetic;
        return runDemo(synthetic, [&](double east_m, double north_m) {
            return synthetic.pointFromMeters(east_m, north_m);
        });
    } catch (const std::exception& e) {
        std::cerr << "demo failed: " << e.what() << "\n";
        return 1;
    }
}

#include "openfhe_dtm/Analysis.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/TdbDataset.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    using namespace openfhe_dtm;
    if (argc != 4 && argc != 5) {
        std::cerr << "usage: final_analysis_demo <dem.tdb> <height_m> <slope_angle_degrees> [grid_interval_m]\n"
                  << "The polygon is the paper's reproducible rectangular DEM subregion.\n"
                  << "For an arbitrary public polygon use EarthworkSlopeQuery in the library.\n";
        return 2;
    }
    try {
        auto number = [](const char* text) {
            std::size_t used = 0;
            const std::string input(text);
            const double value = std::stod(input, &used);
            if (used != input.size() || !std::isfinite(value))
                throw std::invalid_argument("Expected a finite numeric argument");
            return value;
        };
        const double height = number(argv[2]);
        const double angle = number(argv[3]);
        if (angle <= 0.0 || angle >= 90.0)
            throw std::invalid_argument("Earthwork slope angle must be between 0 and 90 degrees");
        const RasterDtmProvider terrain = loadValidationTdbDataset(argv[1]);
        const double span = std::max(2.0, std::min(terrain.widthMeters(), terrain.heightMeters()) * 0.55);
        const double interval = argc == 5 ? number(argv[4]) : std::max(2.0, span / 18.0);
        if (interval <= 0.0) throw std::invalid_argument("Grid interval must be positive");
        const double x0 = span * 0.12, y0 = span * 0.12;
        const auto point = [&](double x, double y) { return terrain.pointFromMeters(x, y); };
        // No saved key directory: every process includes fresh context/key generation.
        OpenFheBackend backend(2048);
        TerrainAnalysis analysis(terrain, backend);
        const auto sight_start = std::chrono::steady_clock::now();
        const auto sight = analysis.sightSurface({point(x0,y0),point(x0+span,y0+span*0.61),
                                                 12.0,span,std::max(2.0,span/9.0)});
        const double sight_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-sight_start).count();
        EarthworkSlopeQuery query;
        query.polygon = {point(x0+span*0.11,y0+span*0.11),point(x0+span*0.78,y0+span*0.11),
                         point(x0+span*0.78,y0+span*0.78),point(x0+span*0.11,y0+span*0.78)};
        query.design_height_m = height;
        query.slope_angle_degrees = angle;
        query.sample_interval_m = interval;
        const auto earth_start = std::chrono::steady_clock::now();
        const auto result = analysis.earthworkSlope(query);
        const double earth_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-earth_start).count();
        if (!result.envelope_valid) {
            std::cerr << "Incompatible boundary heights and slope: " << result.status
                      << " gap_m=" << result.max_envelope_gap_m << '\n';
            return 3;
        }
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << "RESULT_JSON {\"height_m\":" << height
                  << ",\"slope_angle_degrees\":" << angle
                  << ",\"requested_grid_interval_m\":" << interval
                  << ",\"sightline_samples\":" << sight.points.size()
                  << ",\"grid_points\":" << result.geometry.points.size()
                  << ",\"grid_cells\":" << result.geometry.cells.size()
                  << ",\"boundary_points\":" << result.geometry.boundary_points.size()
                  << ",\"integration_area_m2\":" << result.geometry.integration_area_m2
                  << ",\"cut_volume_m3\":" << result.cut_volume_m3
                  << ",\"fill_volume_m3\":" << result.fill_volume_m3
                  << ",\"net_fill_minus_cut_m3\":" << result.net_volume_m3
                  << ",\"sightline_seconds_excluding_keys\":" << sight_seconds
                  << ",\"earthwork_seconds_excluding_keys\":" << earth_seconds
                  << ",\"slots\":2048,\"fresh_keys\":true}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

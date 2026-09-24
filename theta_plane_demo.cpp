#include "openfhe_dtm/EarthworkTheta.hpp"
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
        std::cerr << "usage: theta_plane_demo <dem.tdb> <height_m> <theta_degrees> [sample_interval_m]\n"
                  << "Theta tilts a plane toward polygon[1] from polygon[0].\n"
                  << "Net volume is FILL minus CUT; theta=0 is a horizontal plane.\n";
        return 2;
    }
    try {
        const auto number = [](const char* text) {
            std::size_t used = 0;
            const std::string s(text);
            const double value = std::stod(s, &used);
            if (used != s.size() || !std::isfinite(value))
                throw std::invalid_argument("Expected a finite number");
            return value;
        };
        const double height = number(argv[2]);
        const double theta = number(argv[3]);
        if (std::abs(theta) >= 90.0)
            throw std::invalid_argument("Theta must be strictly between -90 and 90 degrees");
        const RasterDtmProvider terrain = loadValidationTdbDataset(argv[1]);
        const double span = std::max(2.0, std::min(terrain.widthMeters(), terrain.heightMeters()) * 0.55);
        const double interval = argc == 5 ? number(argv[4]) : std::max(2.0, span / 18.0);
        const double x0 = span * 0.12, y0 = span * 0.12;
        const auto point = [&](double x, double y) { return terrain.pointFromMeters(x, y); };
        EarthworkThetaQuery query;
        query.polygon = {point(x0+span*0.11,y0+span*0.11),point(x0+span*0.78,y0+span*0.11),
                         point(x0+span*0.78,y0+span*0.78),point(x0+span*0.11,y0+span*0.78)};
        query.design_height_m = height;
        query.slope_angle_degrees = theta;
        query.sample_interval_m = interval;
        OpenFheBackend backend(2048);
        const auto start = std::chrono::steady_clock::now();
        const auto result = earthworkTheta(query, terrain, backend);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << "RESULT_JSON {\"height_m\":" << height
                  << ",\"theta_degrees\":" << theta
                  << ",\"sample_interval_m\":" << interval
                  << ",\"samples\":" << result.sample_count
                  << ",\"sample_area_m2\":" << result.sample_area_m2
                  << ",\"total_area_m2\":" << result.total_area_m2
                  << ",\"direction_east\":" << result.direction_east
                  << ",\"direction_north\":" << result.direction_north
                  << ",\"net_fill_minus_cut_m3\":" << result.net_volume_m3
                  << ",\"seconds_excluding_keys\":" << seconds
                  << ",\"slots\":2048}\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

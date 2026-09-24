#include "openfhe_dtm/TdbDataset.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

double parseDouble(const std::string& value, const std::string& name) {
    std::size_t parsed = 0;
    const double number = std::stod(value, &parsed);
    if (parsed != value.size()) {
        throw std::invalid_argument(name + " must be numeric: " + value);
    }
    return number;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: convert_neh_to_validation_tdb <input-neh.csv> <output.tdb> "
                     "[--origin lon lat] [--bounds minE maxE minN maxN]\n";
        return 2;
    }

    try {
        openfhe_dtm::Twd97NehCsvOptions options;
        for (int i = 3; i < argc; ++i) {
            const std::string arg = argv[i];
            auto needValue = [&](const char* name) -> const char* {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(std::string(name) + " requires a value");
                }
                return argv[++i];
            };

        if (arg == "--origin") {
            options.local_origin.lon = parseDouble(needValue("--origin lon"), "origin lon");
            options.local_origin.lat = parseDouble(needValue("--origin lat"), "origin lat");
        } else if (arg == "--bounds") {
            options.use_bounds = true;
                options.min_east = parseDouble(needValue("--bounds minE"), "bounds minE");
                options.max_east = parseDouble(needValue("--bounds maxE"), "bounds maxE");
            options.min_north = parseDouble(needValue("--bounds minN"), "bounds minN");
            options.max_north = parseDouble(needValue("--bounds maxN"), "bounds maxN");
        } else if (arg == "--allow-incomplete-grid") {
            options.allow_incomplete_grid = true;
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
        }

        if (options.use_bounds &&
            (options.min_east > options.max_east || options.min_north > options.max_north)) {
            throw std::invalid_argument("invalid bounds: min values must not exceed max values");
        }

        openfhe_dtm::writeValidationTdbFromTwd97NehCsv(argv[1], argv[2], options);
        openfhe_dtm::RasterDtmProvider converted =
            openfhe_dtm::loadValidationTdbDataset(argv[2]);
        std::ofstream metadata(std::string(argv[2]) + ".meta.txt");
        if (!metadata) {
            throw std::runtime_error("could not create conversion metadata sidecar");
        }
        metadata << std::setprecision(17);
        metadata << "source_csv=" << argv[1] << "\n";
        metadata << "output_tdb=" << argv[2] << "\n";
        metadata << "source_crs=TWD97 N,E,H\n";
        metadata << "grid_width_cells=" << converted.widthCells() << "\n";
        metadata << "grid_height_cells=" << converted.heightCells() << "\n";
        metadata << "grid_cell_size_m=" << converted.cellSizeMeters() << "\n";
        metadata << "local_origin_lon=" << options.local_origin.lon << "\n";
        metadata << "local_origin_lat=" << options.local_origin.lat << "\n";
        metadata << "use_bounds=" << (options.use_bounds ? "true" : "false") << "\n";
        metadata << "allow_incomplete_grid="
                 << (options.allow_incomplete_grid ? "true" : "false") << "\n";
        if (options.use_bounds) {
            metadata << "min_east=" << options.min_east << "\n";
            metadata << "max_east=" << options.max_east << "\n";
            metadata << "min_north=" << options.min_north << "\n";
            metadata << "max_north=" << options.max_north << "\n";
        }
        std::cout << "converted TWD97 N/E/H CSV to validation .tdb: " << argv[2] << "\n";
        std::cout << "wrote conversion metadata: " << argv[2] << ".meta.txt\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "failed to convert TWD97 N/E/H CSV: " << e.what() << "\n";
        return 2;
    }
}

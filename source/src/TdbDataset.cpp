#include "openfhe_dtm/TdbDataset.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace openfhe_dtm {
namespace {

bool endsWithTdb(std::string path) {
    std::transform(path.begin(), path.end(), path.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return path.size() >= 4 && path.substr(path.size() - 4) == ".tdb";
}

double syntheticFixtureElevation(double east_m, double north_m) {
    const double plane = 100.0 + 0.08 * east_m + 0.04 * north_m;
    const double dx = east_m - 45.0;
    const double dy = north_m - 35.0;
    const double hill = 18.0 * std::exp(-(dx * dx + dy * dy) / 900.0);
    const double trench = -6.0 * std::exp(-std::pow(north_m - 75.0, 2.0) / 160.0);
    return plane + hill + trench;
}

struct NehSample {
    double north = 0.0;
    double east = 0.0;
    double height = 0.0;
};

std::vector<double> uniqueSortedCoordinates(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(),
                             [](double a, double b) { return std::abs(a - b) < 1e-6; }),
                 values.end());
    return values;
}

double inferCellSize(const std::vector<double>& xs, const std::vector<double>& ys) {
    double cell = std::numeric_limits<double>::max();
    auto scan = [&](const std::vector<double>& values) {
        for (std::size_t i = 1; i < values.size(); ++i) {
            const double d = values[i] - values[i - 1];
            if (d > 1e-6) {
                cell = std::min(cell, d);
            }
        }
    };
    scan(xs);
    scan(ys);
    if (!std::isfinite(cell) || cell == std::numeric_limits<double>::max()) {
        throw std::runtime_error("could not infer N/E/H grid cell size");
    }
    return cell;
}

std::size_t coordinateIndex(const std::vector<double>& values, double value) {
    const auto it = std::lower_bound(values.begin(), values.end(), value - 1e-6);
    if (it == values.end() || std::abs(*it - value) > 1e-5) {
        throw std::runtime_error("N/E/H coordinate could not be mapped to grid");
    }
    return static_cast<std::size_t>(std::distance(values.begin(), it));
}

bool sampleInBounds(const NehSample& sample, const Twd97NehCsvOptions& options) {
    if (!options.use_bounds) {
        return true;
    }
    return sample.east >= options.min_east && sample.east <= options.max_east &&
           sample.north >= options.min_north && sample.north <= options.max_north;
}

std::vector<NehSample> readTwd97NehCsv(const std::string& path,
                                       const Twd97NehCsvOptions& options) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("could not open TWD97 N/E/H CSV");
    }

    std::vector<NehSample> samples;
    std::string line;
    while (std::getline(input, line)) {
        for (char& c : line) {
            if (c == ',' || c == ';' || c == '\t') {
                c = ' ';
            }
        }

        std::istringstream row(line);
        NehSample sample;
        if ((row >> sample.north >> sample.east >> sample.height) &&
            sampleInBounds(sample, options)) {
            samples.push_back(sample);
        }
    }
    if (samples.empty()) {
        throw std::runtime_error("TWD97 N/E/H CSV contained no numeric N,E,H rows");
    }
    return samples;
}

RasterGrid gridFromTwd97NehCsv(const std::string& path, const Twd97NehCsvOptions& options) {
    const std::vector<NehSample> samples = readTwd97NehCsv(path, options);

    std::vector<double> easts;
    std::vector<double> norths;
    easts.reserve(samples.size());
    norths.reserve(samples.size());
    for (const NehSample& sample : samples) {
        easts.push_back(sample.east);
        norths.push_back(sample.north);
    }
    easts = uniqueSortedCoordinates(std::move(easts));
    norths = uniqueSortedCoordinates(std::move(norths));
    if (easts.size() < 2 || norths.size() < 2) {
        throw std::runtime_error("TWD97 N/E/H CSV must contain at least a 2 by 2 grid");
    }

    RasterGrid grid;
    grid.origin = options.local_origin;
    grid.width = easts.size();
    grid.height = norths.size();
    grid.cell_size_m = inferCellSize(easts, norths);
    grid.nodata = -9999.0;
    grid.elevations_m.assign(grid.width * grid.height, grid.nodata);

    for (const NehSample& sample : samples) {
        const std::size_t x = coordinateIndex(easts, sample.east);
        const std::size_t y = coordinateIndex(norths, sample.north);
        grid.elevations_m[y * grid.width + x] = sample.height;
    }
    if (!options.allow_incomplete_grid) {
        const auto missing = std::count(grid.elevations_m.begin(), grid.elevations_m.end(),
                                        grid.nodata);
        if (missing != 0) {
            throw std::runtime_error("TWD97 N/E/H CSV does not form a complete grid");
        }
    }
    return grid;
}

void writeValidationTdbGrid(const std::string& path, const RasterGrid& grid) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not create .tdb fixture");
    }

    output << "OPENFHE_DTM_VALIDATION_TDB_V1\n";
    output << std::setprecision(17) << grid.origin.lon << ' ' << grid.origin.lat << ' '
           << grid.width << ' ' << grid.height << ' ' << grid.cell_size_m << ' '
           << grid.nodata << '\n';
    for (double value : grid.elevations_m) {
        output << std::setprecision(17) << value << '\n';
    }
}

} // namespace

TdbDatasetInfo inspectTdbDataset(const std::string& path) {
    TdbDatasetInfo info;
    info.path = path;
    info.extension_ok = endsWithTdb(path);

    try {
        const std::filesystem::path file_path(path);
        info.exists = std::filesystem::exists(file_path);
        if (!info.exists) {
            info.error = "file does not exist";
            return info;
        }
        if (!std::filesystem::is_regular_file(file_path)) {
            info.error = "path is not a regular file";
            return info;
        }

        info.size_bytes = std::filesystem::file_size(file_path);
        std::ifstream input(file_path, std::ios::binary);
        if (!input) {
            info.error = "file could not be opened";
            return info;
        }

        constexpr std::size_t kProbeBytes = 32;
        info.first_bytes.resize(kProbeBytes);
        input.read(reinterpret_cast<char*>(info.first_bytes.data()),
                   static_cast<std::streamsize>(info.first_bytes.size()));
        info.first_bytes.resize(static_cast<std::size_t>(input.gcount()));
    } catch (const std::exception& e) {
        info.error = e.what();
    }

    return info;
}

std::string formatTdbDatasetInfo(const TdbDatasetInfo& info) {
    std::ostringstream out;
    out << "tdb path: " << info.path << '\n';
    out << "extension .tdb: " << (info.extension_ok ? "yes" : "no") << '\n';
    out << "exists: " << (info.exists ? "yes" : "no") << '\n';
    if (info.exists) {
        out << "size bytes: " << info.size_bytes << '\n';
    }
    if (!info.first_bytes.empty()) {
        out << "first bytes:";
        for (unsigned char byte : info.first_bytes) {
            out << ' ' << std::hex << std::setw(2) << std::setfill('0')
                << static_cast<int>(byte);
        }
        out << std::dec << '\n';
    }
    if (!info.error.empty()) {
        out << "preflight note: " << info.error << '\n';
    }
    return out.str();
}

RasterDtmProvider loadValidationTdbDataset(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("could not open .tdb file");
    }

    std::string magic;
    input >> magic;
    if (magic != "OPENFHE_DTM_VALIDATION_TDB_V1") {
        throw std::runtime_error(
            "unsupported .tdb fixture; original DemDatabase format still needs implementation");
    }

    double lon = 0.0;
    double lat = 0.0;
    std::size_t width = 0;
    std::size_t height = 0;
    double cell_size_m = 0.0;
    double nodata = -9999.0;
    input >> lon >> lat >> width >> height >> cell_size_m >> nodata;
    if (!input) {
        throw std::runtime_error("invalid .tdb grid header");
    }
    if (width < 2 || height < 2) {
        throw std::runtime_error("invalid .tdb grid dimensions: width and height must be at least 2");
    }
    if (width > std::numeric_limits<std::size_t>::max() / height) {
        throw std::runtime_error("invalid .tdb grid dimensions: width * height overflows size_t");
    }
    const std::size_t cell_count = width * height;
    constexpr std::size_t kMaxValidationCells = 100000000;
    if (cell_count > kMaxValidationCells) {
        throw std::runtime_error(
            "invalid .tdb grid dimensions: cell count exceeds the 100,000,000-cell safety limit");
    }

    RasterGrid grid;
    grid.origin = GeoPoint{lon, lat};
    grid.width = width;
    grid.height = height;
    grid.cell_size_m = cell_size_m;
    grid.nodata = nodata;
    grid.elevations_m.resize(cell_count);
    for (double& value : grid.elevations_m) {
        input >> value;
        if (!input) {
            throw std::runtime_error("truncated .tdb fixture");
        }
    }

    DtmMetadata metadata{
        "Validation TDB DEM",
        "OpenFHE DEM validation grid",
        "Local tangent plane derived from WGS84 lon/lat",
        "Synthetic metres",
        cell_size_m,
        nodata,
    };
    return RasterDtmProvider(metadata, std::move(grid));
}

void writeValidationTdbDataset(const std::string& path) {
    constexpr double kMarginM = 64.0;
    constexpr std::size_t kWidth = 256;
    constexpr std::size_t kHeight = 256;
    constexpr double kCellSizeM = 1.0;
    constexpr double kNoData = -9999.0;
    const GeoPoint analysis_origin{-1.4700, 53.3800};
    const GeoPoint grid_origin = RasterDtmProvider(
        DtmMetadata{"tmp", "tmp", "tmp", "tmp", 1.0, kNoData},
        RasterGrid{analysis_origin, 2, 2, 1.0, kNoData, {0.0, 0.0, 0.0, 0.0}})
                                     .pointFromMeters(-kMarginM, -kMarginM);

    RasterGrid grid;
    grid.origin = grid_origin;
    grid.width = kWidth;
    grid.height = kHeight;
    grid.cell_size_m = kCellSizeM;
    grid.nodata = kNoData;
    grid.elevations_m.reserve(grid.width * grid.height);
    for (std::size_t y = 0; y < kHeight; ++y) {
        for (std::size_t x = 0; x < kWidth; ++x) {
            grid.elevations_m.push_back(
                syntheticFixtureElevation(static_cast<double>(x) - kMarginM,
                                          static_cast<double>(y) - kMarginM));
        }
    }
    writeValidationTdbGrid(path, grid);
}

RasterDtmProvider loadTwd97NehCsvDataset(const std::string& path,
                                         const GeoPoint& local_origin) {
    Twd97NehCsvOptions options;
    options.local_origin = local_origin;
    return loadTwd97NehCsvDataset(path, options);
}

RasterDtmProvider loadTwd97NehCsvDataset(const std::string& path,
                                         const Twd97NehCsvOptions& options) {
    RasterGrid grid = gridFromTwd97NehCsv(path, options);
    DtmMetadata metadata{
        "Taiwan TWD97 N/E/H DTM",
        "Taiwan 20m DTM N,E,H CSV converted for OpenFHE validation",
        "TWD97 source coordinates normalized into a local tangent plane",
        "Source H metres",
        grid.cell_size_m,
        grid.nodata,
    };
    return RasterDtmProvider(metadata, std::move(grid));
}

void writeValidationTdbFromTwd97NehCsv(const std::string& input_csv_path,
                                       const std::string& output_tdb_path,
                                       const GeoPoint& local_origin) {
    Twd97NehCsvOptions options;
    options.local_origin = local_origin;
    writeValidationTdbFromTwd97NehCsv(input_csv_path, output_tdb_path, options);
}

void writeValidationTdbFromTwd97NehCsv(const std::string& input_csv_path,
                                       const std::string& output_tdb_path,
                                       const Twd97NehCsvOptions& options) {
    const RasterGrid grid = gridFromTwd97NehCsv(input_csv_path, options);
    writeValidationTdbGrid(output_tdb_path, grid);
}

} // namespace openfhe_dtm

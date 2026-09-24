#pragma once

#include "openfhe_dtm/DtmProvider.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace openfhe_dtm {

struct TdbDatasetInfo {
    std::string path;
    bool exists = false;
    bool extension_ok = false;
    std::uintmax_t size_bytes = 0;
    std::vector<unsigned char> first_bytes;
    std::string error;
};

struct Twd97NehCsvOptions {
    GeoPoint local_origin{121.50, 25.05};
    bool use_bounds = false;
    double min_east = 0.0;
    double max_east = 0.0;
    double min_north = 0.0;
    double max_north = 0.0;
    bool allow_incomplete_grid = false;
};

TdbDatasetInfo inspectTdbDataset(const std::string& path);
std::string formatTdbDatasetInfo(const TdbDatasetInfo& info);
RasterDtmProvider loadValidationTdbDataset(const std::string& path);
void writeValidationTdbDataset(const std::string& path);
RasterDtmProvider loadTwd97NehCsvDataset(const std::string& path,
                                         const GeoPoint& local_origin);
RasterDtmProvider loadTwd97NehCsvDataset(const std::string& path,
                                         const Twd97NehCsvOptions& options);
void writeValidationTdbFromTwd97NehCsv(const std::string& input_csv_path,
                                       const std::string& output_tdb_path,
                                       const GeoPoint& local_origin);
void writeValidationTdbFromTwd97NehCsv(const std::string& input_csv_path,
                                       const std::string& output_tdb_path,
                                       const Twd97NehCsvOptions& options);

} // namespace openfhe_dtm

#include "openfhe_dtm/EncryptedDem.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/SlopeGrid.hpp"
#include "openfhe_dtm/TdbDataset.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    using namespace openfhe_dtm;
    if (argc != 3 && argc != 4) {
        std::cerr << "usage: slope_grid_validation <validation.tdb> <cipher-store-dir> "
                     "[--bootstrap]\n";
        return 2;
    }
    try {
        auto dtm = loadValidationTdbDataset(argv[1]);
        const std::filesystem::path store_dir(argv[2]);
        const bool use_bootstrap = argc == 4 && std::string(argv[3]) == "--bootstrap";
        if (argc == 4 && !use_bootstrap)
            throw std::invalid_argument("unknown option; expected --bootstrap");
        OpenFheBackend he(kDemPackSlots,
                          store_dir / (use_bootstrap ? "keys-bootstrap" : "keys"),
                          use_bootstrap);
        SlopeGridAnalysis slope(he);
        const auto expected = SlopeGridAnalysis::plainFull(
            dtm, dtm.cellSizeMeters(), dtm.cellSizeMeters());
        double maximum_tan_squared = 0.0;
        for (double degrees : expected) {
            const double tangent = std::tan(degrees * 3.14159265358979323846 / 180.0);
            const double tan_squared = tangent * tangent;
            maximum_tan_squared = std::max(maximum_tan_squared, tan_squared);
            if (!std::isfinite(tan_squared) || tan_squared > 9.0 + 1e-12)
                throw std::out_of_range("full DEM slope input is outside [0,9]");
        }
        std::cout << "Full DEM: " << dtm.widthCells() << "x" << dtm.heightCells()
                  << ", output slope grid=" << (dtm.widthCells() - 1) << "x"
                  << (dtm.heightCells() - 1) << "\n";
        std::cout << "Slope-domain preflight: tan^2 max=" << maximum_tan_squared
                  << " within [0,9]\n";
        const auto encrypted = slope.computeFullEncrypted(
            dtm, dtm.cellSizeMeters(), dtm.cellSizeMeters(),
            [](std::size_t done, std::size_t total) {
                std::cout << "\rEncrypted slope progress: " << done << "/" << total
                          << std::flush;
                if (done == total) std::cout << "\n";
            });
        const auto actual = slope.decrypt(encrypted);
        const auto original_expected = SlopeGridAnalysis::approximationPlainFull(
            dtm, dtm.cellSizeMeters(), dtm.cellSizeMeters());
        const auto acc = SlopeGridAnalysis::accuracy(actual, expected);
        const auto original_acc = SlopeGridAnalysis::accuracy(actual, original_expected);
        const auto range = std::minmax_element(actual.begin(), actual.end());
        std::size_t block_count = 0;
        for (const auto& pack : encrypted.packs) block_count += pack.blocks.size();
        std::cout << "Encrypted DEM packs: " << encrypted.packs.size()
                  << " ciphertexts, " << block_count
                  << " overlapping 128x128 blocks, slots=" << he.slotCount() << "\n";
        std::cout << "CKKS bootstrap: "
                  << (use_bootstrap ? "enabled and executed before slope evaluation"
                                    : "disabled (leveled memory-constrained mode)")
                  << "\n";
        std::cout << "Encrypted full slope grid: " << encrypted.height << "x"
                  << encrypted.width << "\n";
        std::cout << "Slope range: " << *range.first << " .. " << *range.second
                  << " degrees\n";
        std::cout << "Slope accuracy: mae=" << acc.mae_degrees
                  << " rmse=" << acc.rmse_degrees
                  << " max=" << acc.max_abs_error_degrees << " degrees\n";
        std::cout << "Fixed-polynomial parity: mae=" << original_acc.mae_degrees
                  << " rmse=" << original_acc.rmse_degrees
                  << " max=" << original_acc.max_abs_error_degrees << " degrees\n";
        return original_acc.mae_degrees < 1e-4 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "slope grid validation failed: " << e.what() << "\n";
        return 1;
    }
}

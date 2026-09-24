#pragma once

#include "openfhe_dtm/EncryptedDem.hpp"

#include <functional>
#include <vector>

namespace openfhe_dtm {

struct EncryptedSlopeGrid {
    CipherVector ciphertext;
    std::size_t height = 127;
    std::size_t width = 127;
    std::size_t stride = 128;
    std::size_t slot_offset = 0;
};

struct SlopeGridBlock {
    std::size_t origin_x = 0;
    std::size_t origin_y = 0;
    std::size_t width = 0;
    std::size_t height = 0;
    std::size_t slot_offset = 0;
};

struct EncryptedSlopeGridPack {
    CipherVector ciphertext;
    std::vector<SlopeGridBlock> blocks;
};

struct EncryptedFullSlopeGrid {
    std::size_t width = 0;
    std::size_t height = 0;
    std::vector<EncryptedSlopeGridPack> packs;
};

struct SlopeGridAccuracy {
    double mae_degrees = 0.0;
    double rmse_degrees = 0.0;
    double max_abs_error_degrees = 0.0;
};

class SlopeGridAnalysis {
public:
    explicit SlopeGridAnalysis(const IHeBackend& he);

    EncryptedSlopeGrid computeEncrypted(const EncryptedDemStitch& stitch,
                                        double dx_m, double dy_m) const;
    EncryptedFullSlopeGrid computeFullEncrypted(const RasterDtmProvider& dtm,
                                                double dx_m, double dy_m,
                                                const std::function<void(std::size_t,
                                                                         std::size_t)>& progress = {}) const;
    std::vector<double> decrypt(const EncryptedSlopeGrid& result) const;
    std::vector<double> decrypt(const EncryptedFullSlopeGrid& result) const;
    static std::vector<double> plain(const IDtmProvider& dtm, int tile_x0, int tile_y0,
                                     double dx_m, double dy_m);
    static std::vector<double> originalApproximationPlain(const IDtmProvider& dtm,
                                                          int tile_x0, int tile_y0,
                                                          double dx_m, double dy_m);
    static std::vector<double> plainFull(const RasterDtmProvider& dtm,
                                         double dx_m, double dy_m);
    static std::vector<double> approximationPlainFull(const RasterDtmProvider& dtm,
                                                      double dx_m, double dy_m);
    static double validatePlainDomain(const IDtmProvider& dtm, int tile_x0, int tile_y0,
                                      double dx_m, double dy_m,
                                      double lower = 0.0, double upper = 9.0);
    static SlopeGridAccuracy accuracy(const std::vector<double>& actual,
                                      const std::vector<double>& expected);

private:
    CipherVector computePackedSlope(const CipherVector& packed,
                                    const std::vector<SlopeGridBlock>& blocks,
                                    double dx_m, double dy_m) const;
    const IHeBackend& he_;
};

} // namespace openfhe_dtm

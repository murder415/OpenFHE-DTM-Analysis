#include "openfhe_dtm/SlopeGrid.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace openfhe_dtm {
namespace {
constexpr double kPi = 3.14159265358979323846;
double slopeDegrees(double tan_squared) {
    return std::atan(std::sqrt(std::max(0.0, tan_squared))) * 180.0 / kPi;
}
const std::vector<double> kSlopeChebyshevCoefficients{
    0.984480847858657, 0.389128952267616, -0.192654979045834,
    0.108788772533503, -0.0667328826144381, 0.0436611583613229,
    -0.0301525172043940, 0.0218058843568705, -0.0163995719905173,
    0.0127468714251661, -0.0101841325143747, 0.00832503226359630,
    -0.00693630407845915, 0.00587221321408751, -0.00503875320300484,
    0.00437344600795537, -0.00383360158875261, 0.00338929921774982,
    -0.00301907691735493, 0.00270721480386285, -0.00244198038184240,
    0.00221446826430340, -0.00201781613242581, 0.00184666378562390,
    -0.00169677303078312, 0.00156475561981120, -0.00144787570022889,
    0.00134390392430357, -0.00125100859767283, 0.00116767313159587,
    -0.00109263305474824, 0.00102482713958088, -0.000963359429964294,
    0.000907469165704412, -0.000856507071376238, 0.000809916237674373,
    -0.000767216850213581, 0.000727993689760658, -0.000691886032107898,
    0.000658579274283074, -0.000627798081830240, 0.000599300679067146,
    -0.000572874080039169, 0.000548330124807809, -0.000525502116550742,
    0.000504241990448123, -0.000484417935143491, 0.000465912317352138,
    -0.000448619987745513, 0.000432446716585909, -0.000417308023995179,
    0.000403127910975335, -0.000389838123669875, 0.000377377094781191,
    -0.000365689447362458, 0.000354725210939111, -0.000344439388308244,
    0.000334791432135666, -0.000325744878914955, 0.000317266964425051,
    -0.000309328359041033, 0.000301902915999062, -0.000294967305820652,
    0.00681836815723948};
double originalChebyshevDegrees(double x) {
    const double t = 2.0 * x / 9.0 - 1.0;
    double t0 = 1.0, t1 = t;
    double sum = kSlopeChebyshevCoefficients[0] + kSlopeChebyshevCoefficients[1] * t1;
    for (std::size_t k = 2; k < kSlopeChebyshevCoefficients.size(); ++k) {
        const double next = 2.0 * t * t1 - t0;
        sum += kSlopeChebyshevCoefficients[k] * next;
        t0 = t1; t1 = next;
    }
    return sum * 180.0 / kPi;
}

std::vector<double> rasterCells(const RasterDtmProvider& dtm) {
    const std::size_t width = dtm.widthCells(), height = dtm.heightCells();
    std::vector<double> cells(width * height, dtm.metadata().nodata);
    const std::size_t tiles_x = (width + 63) / 64, tiles_y = (height + 63) / 64;
    for (std::size_t ty = 0; ty < tiles_y; ++ty) {
        for (std::size_t tx = 0; tx < tiles_x; ++tx) {
            const auto tile = dtm.loadTile(static_cast<int>(tx), static_cast<int>(ty));
            for (std::size_t y = 0; y < 64; ++y) {
                const std::size_t gy = ty * 64 + y;
                if (gy >= height) break;
                for (std::size_t x = 0; x < 64; ++x) {
                    const std::size_t gx = tx * 64 + x;
                    if (gx >= width) break;
                    cells[gy * width + gx] = tile.elevations_m[y * tile.width + x];
                }
            }
        }
    }
    return cells;
}
}

SlopeGridAnalysis::SlopeGridAnalysis(const IHeBackend& he) : he_(he) {
    if (he.slotCount() < kDemPackSlots)
        throw std::invalid_argument("packed slope grid requires 32768 CKKS slots");
}

EncryptedSlopeGrid SlopeGridAnalysis::computeEncrypted(const EncryptedDemStitch& stitch,
                                                       double dx_m, double dy_m) const {
    const std::size_t valid_h = std::min<std::size_t>(
        127, stitch.valid_height > 0 ? stitch.valid_height - 1 : 0);
    const std::size_t valid_w = std::min<std::size_t>(
        127, stitch.valid_width > 0 ? stitch.valid_width - 1 : 0);
    const std::vector<SlopeGridBlock> blocks{{0, 0, valid_w, valid_h,
                                               stitch.slot_offset}};
    auto degrees = computePackedSlope(stitch.ciphertext, blocks, dx_m, dy_m);
    return EncryptedSlopeGrid{std::move(degrees), 127, 127, 128,
                              stitch.slot_offset};
}

CipherVector SlopeGridAnalysis::computePackedSlope(
    const CipherVector& packed, const std::vector<SlopeGridBlock>& blocks,
    double dx_m, double dy_m) const {
    if (dx_m <= 0.0 || dy_m <= 0.0)
        throw std::invalid_argument("slope grid spacing must be positive");
    if (blocks.empty() || blocks.size() > 2)
        throw std::invalid_argument("a slope pack must contain one or two blocks");

    // When the backend was created with bootstrap enabled, match the original
    // execution boundary by refreshing the stored DEM ciphertext first.  The
    // leveled fallback remains available for memory-constrained validation.
    const auto refreshed = he_.supportsBootstrap()
        ? he_.bootstrap(packed)
        : packed;
    const auto east = he_.rotate(refreshed, 1);
    const auto south_east = he_.rotate(refreshed, 128);
    auto u = he_.sub(refreshed, east);
    auto v = he_.rotate(he_.sub(refreshed, south_east), 1);

    std::vector<double> u_scale(he_.slotCount(), kDemScale / dx_m);
    std::vector<double> v_scale(he_.slotCount(), kDemScale / dy_m);
    u = he_.mulPlain(u, he_.encode(u_scale));
    v = he_.mulPlain(v, he_.encode(v_scale));
    auto tan_squared = he_.add(he_.mul(u, u), he_.mul(v, v));

    std::vector<double> mask(he_.slotCount(), 0.0);
    for (const auto& block : blocks) {
        if (block.slot_offset + kDemStitchSlots > he_.slotCount() ||
            block.width > 127 || block.height > 127)
            throw std::invalid_argument("invalid packed slope block");
        for (std::size_t y = 0; y < block.height; ++y)
            for (std::size_t x = 0; x < block.width; ++x)
                mask[block.slot_offset + y * 128 + x] = 2.0 / 9.0;
    }
    auto normalized = he_.mulPlain(tan_squared, he_.encode(mask));
    std::vector<double> ones(he_.slotCount(), 1.0);
    normalized = he_.subPlain(normalized, he_.encode(ones));
    auto coefficients = kSlopeChebyshevCoefficients;
    coefficients[0] *= 2.0; // OpenFHE evaluates the constant term as c0 / 2.
    auto degrees = he_.evaluateChebyshevSeries(normalized, coefficients, -1.0, 1.0);
    degrees = he_.mulScalar(degrees, {180.0 / kPi, 0.0});
    for (double& value : mask) value = value == 0.0 ? 0.0 : 1.0;
    degrees = he_.mulPlain(degrees, he_.encode(mask));
    return degrees;
}

EncryptedFullSlopeGrid SlopeGridAnalysis::computeFullEncrypted(
    const RasterDtmProvider& dtm, double dx_m, double dy_m,
    const std::function<void(std::size_t, std::size_t)>& progress) const {
    if (dtm.widthCells() < 2 || dtm.heightCells() < 2)
        throw std::invalid_argument("full slope grid requires at least 2x2 DEM cells");

    const auto cells = rasterCells(dtm);
    const std::size_t source_width = dtm.widthCells(), source_height = dtm.heightCells();
    const std::size_t output_width = source_width - 1, output_height = source_height - 1;
    const std::size_t blocks_x = (output_width + 126) / 127;
    const std::size_t blocks_y = (output_height + 126) / 127;
    const std::size_t total_packs = ((blocks_x + 1) / 2) * blocks_y;
    std::size_t completed_packs = 0;
    EncryptedFullSlopeGrid result{output_width, output_height, {}};
    result.packs.reserve(total_packs);

    for (std::size_t by = 0; by < blocks_y; ++by) {
        for (std::size_t bx = 0; bx < blocks_x; bx += 2) {
            std::vector<double> packed(he_.slotCount(), 0.0);
            std::vector<SlopeGridBlock> blocks;
            for (std::size_t half = 0; half < 2 && bx + half < blocks_x; ++half) {
                const std::size_t origin_x = (bx + half) * 127;
                const std::size_t origin_y = by * 127;
                const std::size_t input_width = std::min<std::size_t>(128, source_width - origin_x);
                const std::size_t input_height = std::min<std::size_t>(128, source_height - origin_y);
                const std::size_t offset = half * kDemStitchSlots;
                for (std::size_t y = 0; y < input_height; ++y) {
                    for (std::size_t x = 0; x < input_width; ++x) {
                        const double value = cells[(origin_y + y) * source_width + origin_x + x];
                        if (value != dtm.metadata().nodata)
                            packed[offset + y * 128 + x] = value / kDemScale;
                    }
                }
                blocks.push_back({origin_x, origin_y, input_width - 1,
                                  input_height - 1, offset});
            }
            auto encrypted_dem = he_.encrypt(he_.encode(packed));
            auto encrypted_slope = computePackedSlope(encrypted_dem, blocks, dx_m, dy_m);
            result.packs.push_back({std::move(encrypted_slope), std::move(blocks)});
            if (progress) progress(++completed_packs, total_packs);
        }
    }
    return result;
}

std::vector<double> SlopeGridAnalysis::decrypt(const EncryptedSlopeGrid& result) const {
    const auto all = he_.decrypt(result.ciphertext).values;
    std::vector<double> out;
    out.reserve(result.height * result.width);
    for (std::size_t y = 0; y < result.height; ++y)
        for (std::size_t x = 0; x < result.width; ++x)
            out.push_back(all[result.slot_offset + y * result.stride + x]);
    return out;
}

std::vector<double> SlopeGridAnalysis::decrypt(const EncryptedFullSlopeGrid& result) const {
    std::vector<double> out(result.width * result.height, 0.0);
    for (const auto& pack : result.packs) {
        const auto values = he_.decrypt(pack.ciphertext).values;
        for (const auto& block : pack.blocks) {
            for (std::size_t y = 0; y < block.height; ++y)
                for (std::size_t x = 0; x < block.width; ++x)
                    out[(block.origin_y + y) * result.width + block.origin_x + x] =
                        values[block.slot_offset + y * 128 + x];
        }
    }
    return out;
}

std::vector<double> SlopeGridAnalysis::plain(const IDtmProvider& dtm, int tile_x0,
                                            int tile_y0, double dx_m, double dy_m) {
    std::vector<double> z(kDemStitchSlots, 0.0);
    for (int ty = 0; ty < 2; ++ty) for (int tx = 0; tx < 2; ++tx) {
        const auto tile = dtm.loadTile(tile_x0 + tx, tile_y0 + ty);
        for (std::size_t y = 0; y < 64; ++y) for (std::size_t x = 0; x < 64; ++x) {
            const double value = tile.elevations_m[y * tile.width + x];
            z[(static_cast<std::size_t>(ty) * 64 + y) * 128 +
              static_cast<std::size_t>(tx) * 64 + x] =
                value == dtm.metadata().nodata ? 0.0 : value;
        }
    }
    std::vector<double> out;
    out.reserve(127 * 127);
    std::size_t valid_w = 127, valid_h = 127;
    if (const auto* raster = dynamic_cast<const RasterDtmProvider*>(&dtm)) {
        const std::size_t sx = static_cast<std::size_t>(std::max(0, tile_x0)) * 64;
        const std::size_t sy = static_cast<std::size_t>(std::max(0, tile_y0)) * 64;
        valid_w = sx < raster->widthCells() ? std::min<std::size_t>(127, raster->widthCells() - sx - 1) : 0;
        valid_h = sy < raster->heightCells() ? std::min<std::size_t>(127, raster->heightCells() - sy - 1) : 0;
    }
    for (std::size_t y = 0; y < 127; ++y) for (std::size_t x = 0; x < 127; ++x) {
        if (x >= valid_w || y >= valid_h) { out.push_back(0.0); continue; }
        const double b = z[y * 128 + x], c = z[y * 128 + x + 1];
        const double d = z[(y + 1) * 128 + x + 1];
        const double ug = (c - b) / dx_m, vg = (d - c) / dy_m;
        out.push_back(slopeDegrees(ug * ug + vg * vg));
    }
    return out;
}

std::vector<double> SlopeGridAnalysis::originalApproximationPlain(
    const IDtmProvider& dtm, int tile_x0, int tile_y0, double dx_m, double dy_m) {
    auto out = plain(dtm, tile_x0, tile_y0, dx_m, dy_m);
    std::size_t valid_w = 127, valid_h = 127;
    if (const auto* raster = dynamic_cast<const RasterDtmProvider*>(&dtm)) {
        const std::size_t sx = static_cast<std::size_t>(std::max(0, tile_x0)) * 64;
        const std::size_t sy = static_cast<std::size_t>(std::max(0, tile_y0)) * 64;
        valid_w = sx < raster->widthCells() ? std::min<std::size_t>(127, raster->widthCells() - sx - 1) : 0;
        valid_h = sy < raster->heightCells() ? std::min<std::size_t>(127, raster->heightCells() - sy - 1) : 0;
    }
    for (std::size_t y = 0; y < 127; ++y) for (std::size_t x = 0; x < 127; ++x) {
        if (x >= valid_w || y >= valid_h) { out[y * 127 + x] = 0.0; continue; }
        const double theta = out[y * 127 + x] * kPi / 180.0;
        const double tangent = std::tan(theta);
        out[y * 127 + x] = originalChebyshevDegrees(tangent * tangent);
    }
    return out;
}

std::vector<double> SlopeGridAnalysis::plainFull(const RasterDtmProvider& dtm,
                                                 double dx_m, double dy_m) {
    if (dx_m <= 0.0 || dy_m <= 0.0 || dtm.widthCells() < 2 || dtm.heightCells() < 2)
        throw std::invalid_argument("invalid full slope grid dimensions or spacing");
    const auto cells = rasterCells(dtm);
    const std::size_t width = dtm.widthCells(), height = dtm.heightCells();
    std::vector<double> out((width - 1) * (height - 1), 0.0);
    for (std::size_t y = 0; y + 1 < height; ++y) {
        for (std::size_t x = 0; x + 1 < width; ++x) {
            const double b = cells[y * width + x], c = cells[y * width + x + 1];
            const double d = cells[(y + 1) * width + x + 1];
            if (b == dtm.metadata().nodata || c == dtm.metadata().nodata ||
                d == dtm.metadata().nodata) continue;
            const double u = (c - b) / dx_m, v = (d - c) / dy_m;
            out[y * (width - 1) + x] = slopeDegrees(u * u + v * v);
        }
    }
    return out;
}

std::vector<double> SlopeGridAnalysis::approximationPlainFull(
    const RasterDtmProvider& dtm, double dx_m, double dy_m) {
    auto out = plainFull(dtm, dx_m, dy_m);
    for (double& degrees : out) {
        const double tangent = std::tan(degrees * kPi / 180.0);
        degrees = originalChebyshevDegrees(tangent * tangent);
    }
    return out;
}

double SlopeGridAnalysis::validatePlainDomain(const IDtmProvider& dtm, int tile_x0,
                                             int tile_y0, double dx_m, double dy_m,
                                             double lower, double upper) {
    if (lower < 0.0 || upper <= lower)
        throw std::invalid_argument("invalid slope approximation domain");
    const auto slopes = plain(dtm, tile_x0, tile_y0, dx_m, dy_m);
    double maximum = 0.0;
    for (double degrees : slopes) {
        const double tangent = std::tan(degrees * kPi / 180.0);
        const double tan_squared = tangent * tangent;
        maximum = std::max(maximum, tan_squared);
        if (!std::isfinite(tan_squared) || tan_squared < lower - 1e-12 ||
            tan_squared > upper + 1e-12) {
            throw std::out_of_range(
                "slope approximation input tan^2(theta) is outside [0,9]; "
                "maximum supported slope is approximately 71.565 degrees");
        }
    }
    return maximum;
}

SlopeGridAccuracy SlopeGridAnalysis::accuracy(const std::vector<double>& actual,
                                              const std::vector<double>& expected) {
    if (actual.size() != expected.size() || actual.empty())
        throw std::invalid_argument("slope grid accuracy shape mismatch");
    double ae = 0.0, se = 0.0, mx = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        const double e = std::abs(actual[i] - expected[i]);
        ae += e; se += e * e; mx = std::max(mx, e);
    }
    const double n = static_cast<double>(actual.size());
    return {ae / n, std::sqrt(se / n), mx};
}

} // namespace openfhe_dtm

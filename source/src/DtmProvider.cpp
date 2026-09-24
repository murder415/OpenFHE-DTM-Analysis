#include "openfhe_dtm/DtmProvider.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace openfhe_dtm {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEarthRadiusM = 6371008.8;

double toRad(double degrees) {
    return degrees * kPi / 180.0;
}

double toDeg(double radians) {
    return radians * 180.0 / kPi;
}

std::pair<double, double> metersFromOrigin(const GeoPoint& origin, const GeoPoint& point) {
    const double lat_rad = toRad(origin.lat);
    const double east = toRad(point.lon - origin.lon) * kEarthRadiusM * std::cos(lat_rad);
    const double north = toRad(point.lat - origin.lat) * kEarthRadiusM;
    return {east, north};
}

GeoPoint pointFromMeters(const GeoPoint& origin, double east_m, double north_m) {
    const double lat_rad = toRad(origin.lat);
    return GeoPoint{
        origin.lon + toDeg(east_m / (kEarthRadiusM * std::cos(lat_rad))),
        origin.lat + toDeg(north_m / kEarthRadiusM),
    };
}

double syntheticBareEarth(double east_m, double north_m) {
    const double plane = 100.0 + 0.08 * east_m + 0.04 * north_m;
    const double dx = east_m - 45.0;
    const double dy = north_m - 35.0;
    const double hill = 18.0 * std::exp(-(dx * dx + dy * dy) / 900.0);
    const double trench = -6.0 * std::exp(-std::pow(north_m - 75.0, 2.0) / 160.0);
    return plane + hill + trench;
}

double sampleGrid(const RasterGrid& grid, std::size_t x, std::size_t y) {
    return grid.elevations_m[y * grid.width + x];
}

RasterGrid makeSyntheticGrid(const GeoPoint& analysis_origin) {
    constexpr double kMarginM = 64.0;
    RasterGrid grid;
    grid.origin = pointFromMeters(analysis_origin, -kMarginM, -kMarginM);
    grid.width = 256;
    grid.height = 256;
    grid.cell_size_m = 1.0;
    grid.nodata = -9999.0;
    grid.elevations_m.reserve(grid.width * grid.height);
    for (std::size_t y = 0; y < grid.height; ++y) {
        for (std::size_t x = 0; x < grid.width; ++x) {
            grid.elevations_m.push_back(
                syntheticBareEarth(static_cast<double>(x) - kMarginM,
                                   static_cast<double>(y) - kMarginM));
        }
    }
    return grid;
}

} // namespace

RasterDtmProvider::RasterDtmProvider(DtmMetadata metadata, RasterGrid grid)
    : metadata_(std::move(metadata)), grid_(std::move(grid)) {
    if (grid_.width < 2 || grid_.height < 2) {
        throw std::invalid_argument("RasterDtmProvider grid must be at least 2 by 2");
    }
    if (grid_.elevations_m.size() != grid_.width * grid_.height) {
        throw std::invalid_argument("RasterDtmProvider grid data size mismatch");
    }
    if (grid_.cell_size_m <= 0.0) {
        throw std::invalid_argument("RasterDtmProvider cell size must be positive");
    }
}

const DtmMetadata& RasterDtmProvider::metadata() const {
    return metadata_;
}

InterpolationStencil RasterDtmProvider::interpolationStencil(const GeoPoint& point) const {
    const auto [east_m, north_m] = metersFromOrigin(grid_.origin, point);
    double x = east_m / grid_.cell_size_m;
    double y = north_m / grid_.cell_size_m;
    constexpr double kBoundsEpsilon = 1e-6;
    const double max_x = static_cast<double>(grid_.width - 1);
    const double max_y = static_cast<double>(grid_.height - 1);
    if (x < -kBoundsEpsilon || y < -kBoundsEpsilon ||
        x > max_x + kBoundsEpsilon || y > max_y + kBoundsEpsilon) {
        throw std::out_of_range("DTM query point is outside raster bounds");
    }
    x = std::clamp(x, 0.0, max_x);
    y = std::clamp(y, 0.0, max_y);

    const std::size_t x0 = static_cast<std::size_t>(std::floor(x));
    const std::size_t y0 = static_cast<std::size_t>(std::floor(y));
    const std::size_t x1 = std::min(x0 + 1, grid_.width - 1);
    const std::size_t y1 = std::min(y0 + 1, grid_.height - 1);
    const double tx = x - static_cast<double>(x0);
    const double ty = y - static_cast<double>(y0);

    InterpolationStencil stencil;
    stencil.z00_m = sampleGrid(grid_, x0, y0);
    stencil.z10_m = sampleGrid(grid_, x1, y0);
    stencil.z01_m = sampleGrid(grid_, x0, y1);
    stencil.z11_m = sampleGrid(grid_, x1, y1);
    stencil.dx = tx;
    stencil.dy = ty;
    if (stencil.z00_m == grid_.nodata || stencil.z10_m == grid_.nodata ||
        stencil.z01_m == grid_.nodata || stencil.z11_m == grid_.nodata) {
        throw std::runtime_error("DTM interpolation touched nodata");
    }
    return stencil;
}

double RasterDtmProvider::elevationAt(const GeoPoint& point) const {
    const InterpolationStencil s = interpolationStencil(point);
    return s.z00_m * (1.0 - s.dx) * (1.0 - s.dy) +
           s.z10_m * s.dx * (1.0 - s.dy) +
           s.z01_m * (1.0 - s.dx) * s.dy +
           s.z11_m * s.dx * s.dy;
}

DtmTile RasterDtmProvider::loadTile(int tile_x, int tile_y) const {
    constexpr std::size_t kTileSize = 64;
    DtmTile tile;
    tile.tile_x = tile_x;
    tile.tile_y = tile_y;
    tile.width = kTileSize;
    tile.height = kTileSize;
    tile.elevations_m.reserve(kTileSize * kTileSize);

    for (std::size_t y = 0; y < kTileSize; ++y) {
        for (std::size_t x = 0; x < kTileSize; ++x) {
            const int source_x = tile_x * static_cast<int>(kTileSize) + static_cast<int>(x);
            const int source_y = tile_y * static_cast<int>(kTileSize) + static_cast<int>(y);
            if (source_x < 0 || source_y < 0 ||
                source_x >= static_cast<int>(grid_.width) ||
                source_y >= static_cast<int>(grid_.height)) {
                tile.elevations_m.push_back(grid_.nodata);
            } else {
                tile.elevations_m.push_back(
                    sampleGrid(grid_, static_cast<std::size_t>(source_x),
                               static_cast<std::size_t>(source_y)));
            }
        }
    }
    return tile;
}

GeoPoint RasterDtmProvider::origin() const {
    return grid_.origin;
}

GeoPoint RasterDtmProvider::pointFromMeters(double east_m, double north_m) const {
    return openfhe_dtm::pointFromMeters(grid_.origin, east_m, north_m);
}

std::size_t RasterDtmProvider::widthCells() const {
    return grid_.width;
}

std::size_t RasterDtmProvider::heightCells() const {
    return grid_.height;
}

double RasterDtmProvider::widthMeters() const {
    return static_cast<double>(grid_.width - 1) * grid_.cell_size_m;
}

double RasterDtmProvider::heightMeters() const {
    return static_cast<double>(grid_.height - 1) * grid_.cell_size_m;
}

double RasterDtmProvider::cellSizeMeters() const {
    return grid_.cell_size_m;
}

const RasterGrid& RasterDtmProvider::grid() const {
    return grid_;
}

SyntheticRasterProvider::SyntheticRasterProvider()
    : origin_{-1.4700, 53.3800},
      raster_(
          DtmMetadata{
              "Synthetic Seoul-like DTM 1m",
              "Deterministic synthetic raster for OpenFHE analysis validation",
              "Local tangent plane derived from WGS84 lon/lat",
              "Synthetic metres",
              1.0,
              -9999.0,
          },
          makeSyntheticGrid(origin_)) {}

const DtmMetadata& SyntheticRasterProvider::metadata() const {
    return raster_.metadata();
}

InterpolationStencil SyntheticRasterProvider::interpolationStencil(
    const GeoPoint& point) const {
    return raster_.interpolationStencil(point);
}

double SyntheticRasterProvider::elevationAt(const GeoPoint& point) const {
    return raster_.elevationAt(point);
}

DtmTile SyntheticRasterProvider::loadTile(int tile_x, int tile_y) const {
    return raster_.loadTile(tile_x, tile_y);
}

GeoPoint SyntheticRasterProvider::origin() const {
    return origin_;
}

GeoPoint SyntheticRasterProvider::pointFromMeters(double east_m, double north_m) const {
    return openfhe_dtm::pointFromMeters(origin_, east_m, north_m);
}

} // namespace openfhe_dtm

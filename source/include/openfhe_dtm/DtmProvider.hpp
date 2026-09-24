#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace openfhe_dtm {

struct GeoPoint {
    double lon = 0.0;
    double lat = 0.0;
};

using Route = std::vector<GeoPoint>;

struct DtmMetadata {
    std::string name;
    std::string source;
    std::string crs;
    std::string vertical_datum;
    double resolution_m = 0.0;
    double nodata = -9999.0;
};

struct DtmTile {
    int tile_x = 0;
    int tile_y = 0;
    std::size_t width = 0;
    std::size_t height = 0;
    std::vector<double> elevations_m;
};

struct RasterGrid {
    GeoPoint origin;
    std::size_t width = 0;
    std::size_t height = 0;
    double cell_size_m = 1.0;
    double nodata = -9999.0;
    std::vector<double> elevations_m;
};

struct InterpolationStencil {
    double z00_m = 0.0;
    double z10_m = 0.0;
    double z01_m = 0.0;
    double z11_m = 0.0;
    double dx = 0.0;
    double dy = 0.0;
};

class IDtmProvider {
public:
    virtual ~IDtmProvider() = default;

    virtual const DtmMetadata& metadata() const = 0;
    virtual InterpolationStencil interpolationStencil(const GeoPoint& point) const = 0;
    virtual double elevationAt(const GeoPoint& point) const = 0;
    virtual DtmTile loadTile(int tile_x, int tile_y) const = 0;
};

class RasterDtmProvider : public IDtmProvider {
public:
    RasterDtmProvider(DtmMetadata metadata, RasterGrid grid);

    const DtmMetadata& metadata() const override;
    InterpolationStencil interpolationStencil(const GeoPoint& point) const override;
    double elevationAt(const GeoPoint& point) const override;
    DtmTile loadTile(int tile_x, int tile_y) const override;

    GeoPoint origin() const;
    GeoPoint pointFromMeters(double east_m, double north_m) const;
    std::size_t widthCells() const;
    std::size_t heightCells() const;
    double widthMeters() const;
    double heightMeters() const;
    double cellSizeMeters() const;

protected:
    const RasterGrid& grid() const;

private:
    DtmMetadata metadata_;
    RasterGrid grid_;
};

class SyntheticRasterProvider final : public IDtmProvider {
public:
    SyntheticRasterProvider();

    const DtmMetadata& metadata() const override;
    InterpolationStencil interpolationStencil(const GeoPoint& point) const override;
    double elevationAt(const GeoPoint& point) const override;
    DtmTile loadTile(int tile_x, int tile_y) const override;

    GeoPoint origin() const;
    GeoPoint pointFromMeters(double east_m, double north_m) const;

private:
    GeoPoint origin_;
    RasterDtmProvider raster_;
};

} // namespace openfhe_dtm

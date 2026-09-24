#include "openfhe_dtm/Analysis.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/RunOptions.hpp"
#include "openfhe_dtm/SlopeGrid.hpp"
#include "openfhe_dtm/TdbDataset.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool ok, const std::string& message) {
    if (!ok) {
        throw std::runtime_error(message);
    }
}

void requireNear(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
                                 " expected=" + std::to_string(expected));
    }
}

} // namespace

int main() {
    using namespace openfhe_dtm;

    try {
        RasterGrid tiny_grid;
        tiny_grid.origin = GeoPoint{0.0, 0.0};
        tiny_grid.width = 3;
        tiny_grid.height = 3;
        tiny_grid.cell_size_m = 1.0;
        tiny_grid.nodata = -9999.0;
        tiny_grid.elevations_m = {
            0.0, 2.0, 4.0,
            2.0, 4.0, 6.0,
            4.0, 6.0, 8.0,
        };
        RasterDtmProvider tiny_dtm(
            DtmMetadata{"Tiny DTM", "unit test", "local tangent", "metres", 1.0, -9999.0},
            tiny_grid);
        requireNear(tiny_dtm.elevationAt(tiny_dtm.pointFromMeters(0.5, 0.5)), 2.0, 0.01,
                    "bilinear interpolation failed");
        const DtmTile tiny_tile = tiny_dtm.loadTile(0, 0);
        require(tiny_tile.elevations_m.size() == 64 * 64, "tile size mismatch");
        requireNear(tiny_tile.elevations_m[0], 0.0, 0.01, "tile origin value mismatch");
        requireNear(tiny_tile.elevations_m[1], 2.0, 0.01, "tile x neighbor value mismatch");
        requireNear(tiny_tile.elevations_m[64], 2.0, 0.01, "tile y neighbor value mismatch");
        requireNear(tiny_tile.elevations_m[3], tiny_grid.nodata, 0.01,
                    "tile padding should use nodata");
        const char* both_args[] = {"synthetic_validation", "--synthetic", "--tdb", "elevation_bil32.tdb"};
        bool rejected_mixed_data_modes = false;
        try {
            (void)parseRunOptions(4, const_cast<char**>(both_args));
        } catch (const std::invalid_argument&) {
            rejected_mixed_data_modes = true;
        }
        require(rejected_mixed_data_modes, "run options should reject mixed synthetic and .tdb modes");
        const char* upper_tdb_args[] = {"synthetic_validation", "--tdb", "ELEVATION_BIL32.TDB"};
        const RunOptions upper_tdb_options = parseRunOptions(3, const_cast<char**>(upper_tdb_args));
        require(upper_tdb_options.tdb_path == "ELEVATION_BIL32.TDB",
                "run options should accept uppercase .TDB extension");
        const TdbDatasetInfo tdb_info = inspectTdbDataset("missing_elevation_bil32.tdb");
        require(tdb_info.extension_ok, ".tdb preflight should recognize extension");
        require(!tdb_info.exists, ".tdb preflight should report missing files");
        const std::filesystem::path fixture_path = "validation_fixture_tmp.tdb";
        writeValidationTdbDataset(fixture_path.string());
        const TdbDatasetInfo fixture_info = inspectTdbDataset(fixture_path.string());
        require(fixture_info.exists, "validation .tdb fixture should exist");
        require(fixture_info.size_bytes > 1000, "validation .tdb fixture should contain raster data");
        RasterDtmProvider fixture_dtm = loadValidationTdbDataset(fixture_path.string());
        requireNear(fixture_dtm.elevationAt(fixture_dtm.pointFromMeters(64.0, 64.0)),
                    100.486, 0.01, "validation .tdb origin sample mismatch");
        std::filesystem::remove(fixture_path);

        const std::filesystem::path oversized_path = "oversized_grid_tmp.tdb";
        {
            std::ofstream oversized(oversized_path);
            oversized << "OPENFHE_DTM_VALIDATION_TDB_V1\n";
            oversized << "0 0 4294967300 2 1 -9999\n";
        }
        bool rejected_oversized_grid = false;
        try {
            (void)loadValidationTdbDataset(oversized_path.string());
        } catch (const std::runtime_error& e) {
            rejected_oversized_grid =
                std::string(e.what()).find("safety limit") != std::string::npos;
        }
        std::filesystem::remove(oversized_path);
        require(rejected_oversized_grid,
                "oversized grid must be rejected before allocation");

        RasterGrid steep_grid = tiny_grid;
        steep_grid.elevations_m = {0, 100, 200, 0, 100, 200, 0, 100, 200};
        RasterDtmProvider steep_dtm(
            DtmMetadata{"Steep DTM", "unit test", "local", "metres", 1, -9999},
            std::move(steep_grid));
        bool rejected_slope_domain = false;
        try {
            (void)SlopeGridAnalysis::validatePlainDomain(steep_dtm, 0, 0, 1, 1);
        } catch (const std::out_of_range&) {
            rejected_slope_domain = true;
        }
        require(rejected_slope_domain, "slope tan-squared values above 9 must be rejected");

        const std::filesystem::path neh_path = "tw_neh_fixture_tmp.csv";
        const std::filesystem::path neh_tdb_path = "tw_neh_fixture_tmp.tdb";
        {
            std::ofstream neh(neh_path);
            neh << "N,E,H\n";
            neh << "2770000,305000,10\n";
            neh << "2770000,305020,12\n";
            neh << "2770020,305000,14\n";
            neh << "2770020,305020,16\n";
        }
        writeValidationTdbFromTwd97NehCsv(neh_path.string(), neh_tdb_path.string(),
                                          GeoPoint{121.50, 25.05});
        RasterDtmProvider neh_dtm = loadValidationTdbDataset(neh_tdb_path.string());
        requireNear(neh_dtm.elevationAt(neh_dtm.pointFromMeters(10.0, 10.0)), 13.0, 0.01,
                    "TWD97 N/E/H validation conversion interpolation mismatch");
        std::filesystem::remove(neh_path);
        std::filesystem::remove(neh_tdb_path);

        {
            std::ofstream neh(neh_path);
            neh << "N,E,H\n";
            neh << "2770000,305000,10\n";
            neh << "2770000,305020,12\n";
            neh << "2770000,305040,14\n";
            neh << "2770020,305000,16\n";
            neh << "2770020,305020,18\n";
            neh << "2770020,305040,20\n";
        }
        Twd97NehCsvOptions neh_options;
        neh_options.local_origin = GeoPoint{121.50, 25.05};
        neh_options.use_bounds = true;
        neh_options.min_east = 305020.0;
        neh_options.max_east = 305040.0;
        neh_options.min_north = 2770000.0;
        neh_options.max_north = 2770020.0;
        writeValidationTdbFromTwd97NehCsv(neh_path.string(), neh_tdb_path.string(),
                                          neh_options);
        RasterDtmProvider clipped_neh_dtm = loadValidationTdbDataset(neh_tdb_path.string());
        requireNear(clipped_neh_dtm.elevationAt(clipped_neh_dtm.pointFromMeters(0.0, 0.0)),
                    12.0, 0.01, "TWD97 N/E/H bounds min corner mismatch");
        requireNear(clipped_neh_dtm.elevationAt(clipped_neh_dtm.pointFromMeters(20.0, 20.0)),
                    20.0, 0.01, "TWD97 N/E/H bounds max corner mismatch");
        std::filesystem::remove(neh_path);
        std::filesystem::remove(neh_tdb_path);

        {
            std::ofstream neh(neh_path);
            neh << "N,E,H\n";
            neh << "2770000,305000,10\n";
            neh << "2770000,305020,12\n";
            neh << "2770020,305000,16\n";
        }
        bool rejected_incomplete_grid = false;
        try {
            writeValidationTdbFromTwd97NehCsv(neh_path.string(), neh_tdb_path.string(),
                                              GeoPoint{121.50, 25.05});
        } catch (const std::runtime_error&) {
            rejected_incomplete_grid = true;
        }
        require(rejected_incomplete_grid, "TWD97 N/E/H conversion should reject incomplete grids");
        std::filesystem::remove(neh_path);
        std::filesystem::remove(neh_tdb_path);

        SyntheticRasterProvider dtm;
        OpenFheBackend he(16);

        TerrainAnalysis tiny_terrain(tiny_dtm, he);
        const GeoPoint tiny_midpoint = tiny_dtm.pointFromMeters(0.5, 0.5);
        const ElevationResult tiny_encrypted = tiny_terrain.elevation(Route{tiny_midpoint});
        requireNear(tiny_encrypted.elevations_m.front(), 2.0, 0.01,
                    "encrypted Method3 bilinear interpolation failed");

        const auto add = he.decrypt(he.add(he.encrypt(he.encode({1.0, 2.0, 3.0, 4.0})),
                                           he.encrypt(he.encode({10.0, 20.0, 30.0, 40.0}))));
        requireNear(add.values[0], 11.0, 0.01, "OpenFHE add slot 0 failed");
        requireNear(add.values[3], 44.0, 0.01, "OpenFHE add slot 3 failed");
        const auto base_ct = he.encrypt(he.encode({1.0, 2.0, 3.0, 4.0}));
        const auto sub = he.decrypt(he.sub(he.encrypt(he.encode({10.0, 20.0, 30.0, 40.0})),
                                           base_ct));
        requireNear(sub.values[0], 9.0, 0.01, "OpenFHE sub slot 0 failed");
        requireNear(sub.values[3], 36.0, 0.01, "OpenFHE sub slot 3 failed");
        const auto mul_plain = he.decrypt(he.mulPlain(base_ct, he.encode({2.0, 3.0, 4.0, 5.0})));
        requireNear(mul_plain.values[0], 2.0, 0.01, "OpenFHE mulPlain slot 0 failed");
        requireNear(mul_plain.values[3], 20.0, 0.01, "OpenFHE mulPlain slot 3 failed");
        const auto rotated = he.decrypt(he.rotate(base_ct, 1));
        requireNear(rotated.values[0], 2.0, 0.01, "OpenFHE rotate direction failed");
        const auto sum_slots = he.decrypt(he.sumSlots(base_ct));
        requireNear(sum_slots.values[0], 10.0, 0.01, "OpenFHE sumSlots slot 0 failed");

        TerrainAnalysis terrain(dtm, he);

        const Route route = {
            dtm.pointFromMeters(0.0, 0.0),
            dtm.pointFromMeters(20.0, 10.0),
            dtm.pointFromMeters(45.0, 35.0),
            dtm.pointFromMeters(75.0, 55.0),
        };

        const ElevationResult elevation = terrain.elevation(route);
        require(elevation.elevations_m.size() == route.size(), "elevation sample count mismatch");
        const ElevationAccuracy elevation_accuracy =
            compareElevation(elevation.elevations_m, terrain.plainElevation(route));
        require(elevation_accuracy.compared == route.size(), "elevation accuracy count mismatch");
        require(elevation_accuracy.max_abs_error_m < 0.05,
                "OpenFHE elevation round trip exceeded 5 cm max error");
        require(elevation.elevations_m.back() > elevation.elevations_m.front(),
                "synthetic elevation should rise along the validation route");

        const SlopeResult slope = terrain.slope(route);
        require(slope.slope_degrees.size() == route.size(), "slope sample count mismatch");
        require(std::abs(slope.slope_degrees.back()) > 1.0, "slope should show measurable relief");
        const SurfaceSlopeResult surface_slope =
            terrain.surfaceSlope(Route{dtm.pointFromMeters(45.0, 35.0)}, 5.0);
        require(surface_slope.slope_degrees.size() == 1,
                "surface slope should return one result per query point");
        require(surface_slope.slope_degrees.front() > 1.0,
                "surface slope should show measurable local relief");

        const TerrainSectionResult section =
            terrain.terrainSection(route, TerrainSectionConfig{40.0, 20.0});
        require(!section.transverse.empty(), "terrain-section should create transverse profiles");
        require(section.transverse.size() == section.transverse_center_distances_m.size(),
                "terrain-section center metadata mismatch");
        require(section.transverse.front().elevations_m.size() >= 3,
                "transverse profile should have sampled points");

        const SightSurfaceResult sight = terrain.sightSurface(
            SightSurfaceQuery{dtm.pointFromMeters(0.0, 0.0), dtm.pointFromMeters(90.0, 55.0), 12.0, 90.0, 10.0});
        require(sight.buildable_height_m.size() == sight.points.size(),
                "sight-surface result vectors should align");
        require(sight.buildable_height_m.size() > 5, "sight-surface should sample the line");
        require(sight.ground_m.empty() && sight.regulation_plane_m.empty(),
                "sight-surface must not expose separately decrypted ground or plane");
        const double sight_base = dtm.elevationAt(sight.points.front());
        for (std::size_t i = 0; i < sight.points.size(); ++i) {
            const double reference = std::max(0.0, sight_base + sight.distances_m[i] *
                std::tan(12.0 * 3.14159265358979323846 / 180.0) - dtm.elevationAt(sight.points[i]));
            require(std::abs(sight.buildable_height_m[i] - reference) < 1e-6,
                    "sight-surface must match an independent terrain reference");
        }

        const Route earthwork_area = {
            dtm.pointFromMeters(10.0, 10.0),
            dtm.pointFromMeters(70.0, 10.0),
            dtm.pointFromMeters(70.0, 70.0),
            dtm.pointFromMeters(10.0, 70.0),
        };
        EarthworkQuery earthwork_query;
        earthwork_query.polygon = earthwork_area;
        earthwork_query.design_height_m = 112.0;
        earthwork_query.sample_interval_m = 5.0;
        const EarthworkResult earthwork = terrain.earthwork(earthwork_query);
        require(earthwork.sample_count > 50, "earthwork should sample the polygon");
        require(earthwork.cut_volume_m3 > 1.0, "earthwork should report cut volume");
        require(earthwork.fill_volume_m3 > 1.0, "earthwork should report fill volume");
        require(std::abs(earthwork.net_volume_m3) > 1.0, "earthwork net volume should be nonzero");

        EarthworkQuery sloped_design_query;
        sloped_design_query.polygon = earthwork_area;
        sloped_design_query.design_height_m = 112.0;
        sloped_design_query.design_slope_east = 0.02;
        sloped_design_query.design_slope_north = -0.01;
        sloped_design_query.sample_interval_m = 5.0;
        const EarthworkResult sloped_earthwork = terrain.earthwork(sloped_design_query);
        require(sloped_earthwork.sample_count == earthwork.sample_count,
                "sloped earthwork should reuse the same sampling grid");
        require(std::abs(sloped_earthwork.net_volume_m3 - earthwork.net_volume_m3) > 1.0,
                "sloped design plane should affect earthwork volume");

        std::cout << "synthetic validation passed\n";
        std::cout << "elevation first/last: " << elevation.elevations_m.front() << " / "
                  << elevation.elevations_m.back() << "\n";
        std::cout << "elevation error mae/rmse/max m: " << elevation_accuracy.mae_m << " / "
                  << elevation_accuracy.rmse_m << " / "
                  << elevation_accuracy.max_abs_error_m << "\n";
        std::cout << "slope last deg: " << slope.slope_degrees.back() << "\n";
        std::cout << "surface slope deg: " << surface_slope.slope_degrees.front() << "\n";
        std::cout << "terrain transverse count: " << section.transverse.size() << "\n";
        std::cout << "sight-surface last buildable m: " << sight.buildable_height_m.back() << "\n";
        std::cout << "earthwork cut/fill/net m3: " << earthwork.cut_volume_m3 << " / "
                  << earthwork.fill_volume_m3 << " / " << earthwork.net_volume_m3 << "\n";
        std::cout << "sloped earthwork net m3: " << sloped_earthwork.net_volume_m3 << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "synthetic validation failed: " << e.what() << "\n";
        return 1;
    }
}

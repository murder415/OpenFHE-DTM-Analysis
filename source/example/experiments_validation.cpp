// 추가 실험용 검증 프로그램 (논문 v1.5 [측정 후 입력] 채우기).
//   실험 1  : 다지점(격자) 고도 보간 오차  -> 표 2 "전 격자 고도 보간 오차"
//   실험 2  : 같은 격자에 대한 평문/암호문 실행시간 비교 -> 표 3 "평문 기준 동일 분석"
//   실험 2b : demo 와 동일한 다각형/설계고의 순 토공량을 암호문 상태로 계산 -> 표 2 "순 토공량(암호문)"
// demo.cpp 의 runDemo 기하(route, earthwork 다각형, sample_interval, 설계고)를 그대로 재현하되
// slot 수만 2048 로 올려 순 토공량 암호문 블록(ground.size() <= slotCount())이 실행되게 한다.
#include "openfhe_dtm/Analysis.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/TdbDataset.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

int main(int argc, char** argv) {
    using namespace openfhe_dtm;
    if (argc != 2) {
        std::cerr << "usage: experiments_validation <dem.tdb>\n";
        return 2;
    }
    try {
        RasterDtmProvider tdb = loadValidationTdbDataset(argv[1]);
        const double span = std::max(
            2.0, std::min(tdb.widthMeters(), tdb.heightMeters()) * 0.55);
        const auto pt = [&](double e, double n) { return tdb.pointFromMeters(e, n); };

        OpenFheBackend he(2048);
        TerrainAnalysis analysis(tdb, he);

        // ---- demo.cpp 와 동일한 route / earthwork 설정 ----
        const double x0 = span * 0.12;
        const double y0 = span * 0.12;
        const Route route = {
            pt(x0, y0),
            pt(x0 + span * 0.22, y0 + span * 0.11),
            pt(x0 + span * 0.50, y0 + span * 0.39),
            pt(x0 + span * 0.83, y0 + span * 0.61),
        };
        const auto route_elev = analysis.elevation(route);
        const double design_height_m =
            std::accumulate(route_elev.elevations_m.begin(),
                            route_elev.elevations_m.end(), 0.0) /
            static_cast<double>(route_elev.elevations_m.size());

        EarthworkQuery q;
        q.polygon = {
            pt(x0 + span * 0.11, y0 + span * 0.11),
            pt(x0 + span * 0.78, y0 + span * 0.11),
            pt(x0 + span * 0.78, y0 + span * 0.78),
            pt(x0 + span * 0.11, y0 + span * 0.78),
        };
        q.design_height_m = design_height_m;
        q.sample_interval_m = std::max(2.0, span / 18.0);
        const auto ew = analysis.earthwork(q);

        std::cout << "EARTHWORK samples=" << ew.sample_count
                  << " design=" << design_height_m
                  << " cut=" << ew.cut_volume_m3
                  << " fill=" << ew.fill_volume_m3
                  << " net_plain=" << ew.net_volume_m3 << " m3\n";
        if (ew.encrypted_net_computed) {
            const double rel = ew.net_volume_abs_error_m3 /
                               std::max(1.0, std::abs(ew.net_volume_m3));
            std::cout << "EARTHWORK net_ciphertext=" << ew.net_volume_encrypted_m3
                      << " m3 abs_error=" << ew.net_volume_abs_error_m3
                      << " m3 rel_error=" << rel << "\n";
        } else {
            std::cout << "EARTHWORK net_ciphertext=SKIPPED (samples "
                      << ew.sample_count << " > slots " << he.slotCount() << ")\n";
        }

        // ---- 실험 1 + 2 : 다지점 격자 고도 보간 오차 & 평문/암호문 시간 ----
        const int side = 23;  // 23 x 23 = 529 표본
        // grid  : 기존(격자와 정렬될 수 있음)
        // grid2 : 셀 대각 방향으로 무리수 오프셋을 줘서 dx,dy 가 확실히 (0,1) 내부
        Route grid, grid2;
        grid.reserve(static_cast<std::size_t>(side) * side);
        grid2.reserve(static_cast<std::size_t>(side) * side);
        for (int gy = 0; gy < side; ++gy) {
            for (int gx = 0; gx < side; ++gx) {
                const double fx = 0.08 + 0.84 * gx / (side - 1);
                const double fy = 0.08 + 0.84 * gy / (side - 1);
                grid.push_back(pt(span * fx, span * fy));
                grid2.push_back(pt(span * fx + 0.5 * tdb.cellSizeMeters() * 1.4142135623730951,
                                   span * fy + 0.5 * tdb.cellSizeMeters() * 0.7071067811865476));
            }
        }
        // 셀 내 위치(dx,dy) 진단
        double dxmin = 1e9, dxmax = -1e9, dxsum = 0.0;
        double dxmin2 = 1e9, dxmax2 = -1e9, dxsum2 = 0.0;
        for (std::size_t i = 0; i < grid.size(); ++i) {
            const auto s  = tdb.interpolationStencil(grid[i]);
            const auto s2 = tdb.interpolationStencil(grid2[i]);
            const double d  = std::abs(s.dx) + std::abs(s.dy);
            const double d2 = std::abs(s2.dx) + std::abs(s2.dy);
            dxmin = std::min(dxmin, d); dxmax = std::max(dxmax, d); dxsum += d;
            dxmin2 = std::min(dxmin2, d2); dxmax2 = std::max(dxmax2, d2); dxsum2 += d2;
        }
        std::cout << "GRIDDIAG |dx|+|dy|  grid : min=" << dxmin << " mean="
                  << dxsum / grid.size() << " max=" << dxmax << "\n";
        std::cout << "GRIDDIAG |dx|+|dy|  grid2: min=" << dxmin2 << " mean="
                  << dxsum2 / grid2.size() << " max=" << dxmax2 << "\n";

        const auto c0 = std::chrono::steady_clock::now();
        const auto grid_enc = analysis.elevation(grid);
        const auto c1 = std::chrono::steady_clock::now();
        const auto grid_ref = analysis.plainElevation(grid);
        const auto c2 = std::chrono::steady_clock::now();
        const auto acc = compareElevation(grid_enc.elevations_m, grid_ref);
        // grid2 (셀 내부 위치) 오차
        const auto acc2 = compareElevation(analysis.elevation(grid2).elevations_m,
                                           analysis.plainElevation(grid2));
        std::cout << "GRID2 (cell-interior) N=" << grid2.size()
                  << " MAE=" << acc2.mae_m << " RMSE=" << acc2.rmse_m
                  << " MAX=" << acc2.max_abs_error_m << " m\n";
        const double enc_ms =
            std::chrono::duration<double, std::milli>(c1 - c0).count();
        const double plain_ms =
            std::chrono::duration<double, std::milli>(c2 - c1).count();
        std::cout << "GRID N=" << grid.size() << " MAE=" << acc.mae_m
                  << " RMSE=" << acc.rmse_m << " MAX=" << acc.max_abs_error_m
                  << " m  enc_ms=" << enc_ms << " plain_ms=" << plain_ms
                  << " ratio=" << (plain_ms > 0.0 ? enc_ms / plain_ms : 0.0)
                  << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "experiments validation failed: " << e.what() << "\n";
        return 1;
    }
}

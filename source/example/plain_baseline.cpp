// 평문 기준선(암호화 없음) — 표 3 "평문 기준 동일 분석" 행 측정용.
// OpenFHE 컨텍스트를 만들지 않고 동일한 격자/다각형에 대해 평문 쌍선형 보간과
// 평문 토공량(절토/성토/순)을 계산한다. 프로세스 시작-종료 전체 시간과 peak WS를
// 외부(Measure-Command / Get-Process)에서 측정한다.
#include "openfhe_dtm/DtmProvider.hpp"
#include "openfhe_dtm/TdbDataset.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
    using namespace openfhe_dtm;
    if (argc != 2) {
        std::cerr << "usage: plain_baseline <dem.tdb>\n";
        return 2;
    }
    try {
        RasterDtmProvider tdb = loadValidationTdbDataset(argv[1]);
        const double span = std::max(
            2.0, std::min(tdb.widthMeters(), tdb.heightMeters()) * 0.55);
        const auto pt = [&](double e, double n) { return tdb.pointFromMeters(e, n); };

        // 529점 격자 평문 보간
        const int side = 23;
        double sum = 0.0;
        for (int gy = 0; gy < side; ++gy) {
            for (int gx = 0; gx < side; ++gx) {
                const double fx = 0.08 + 0.84 * gx / (side - 1);
                const double fy = 0.08 + 0.84 * gy / (side - 1);
                sum += tdb.elevationAt(pt(span * fx, span * fy));
            }
        }

        // demo 와 동일한 다각형/설계고의 평문 토공량
        const double x0 = span * 0.12, y0 = span * 0.12;
        const double dh = tdb.elevationAt(pt(x0, y0));  // 대략치(설계고 자리표시)
        const double si = std::max(2.0, span / 18.0);
        const double minx = x0 + span * 0.11, maxx = x0 + span * 0.78;
        const double miny = y0 + span * 0.11, maxy = y0 + span * 0.78;
        double cut = 0.0, fill = 0.0;
        long n = 0;
        for (double y = miny + si / 2.0; y <= maxy; y += si) {
            for (double x = minx + si / 2.0; x <= maxx; x += si) {
                const double g = tdb.elevationAt(pt(x, y));
                const double d = (g - dh) * si * si;
                if (g - dh > 0.0) cut += d; else fill += -d;
                ++n;
            }
        }
        std::cout << "PLAIN grid_sum=" << sum << " earthwork_samples=" << n
                  << " cut=" << cut << " fill=" << fill
                  << " net=" << (cut - fill) << " m3\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "plain baseline failed: " << e.what() << "\n";
        return 1;
    }
}

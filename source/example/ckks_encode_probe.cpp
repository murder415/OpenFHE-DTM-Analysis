// "초반 근사" = CKKS 인코딩(실수 -> scaling factor 곱 -> 정수 반올림 -> 다항식) 오차만 분리 측정.
//   암호화 X, 동형연산 X.  MakeCKKSPackedPlaintext -> GetRealPackedValue 왕복.
//   현실적인 값: DEM 고도(비정수 소수) + 비정수 쌍선형 가중치.
//   기준선: __float128.

#include "openfhe.h"
#include <quadmath.h>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace lbcrypto;
using q = __float128;

static void sweepN(uint32_t modsize) {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> elev(-6.5, 1114.0);   // 대만 타일 고도 범위
    std::uniform_real_distribution<double> frac(0.0, 1.0);

    std::printf("  scaleModSize = %u\n", modsize);
    std::printf("    N     | MAX 인코딩 abs(m) | MAX rel      | RMS rel\n");
    std::printf("    ------+-------------------+--------------+---------\n");
    for (uint32_t N : {1u, 4u, 16u, 64u, 256u, 529u, 2048u}) {
        CCParams<CryptoContextCKKSRNS> p;
        p.SetMultiplicativeDepth(2);
        p.SetScalingModSize(modsize);
        p.SetBatchSize(4096);
        p.SetCKKSDataType(COMPLEX);
        p.SetSecurityLevel(HEStd_128_classic);
        auto cc = GenCryptoContext(p);
        cc->Enable(PKE); cc->Enable(LEVELEDSHE); cc->Enable(KEYSWITCH);

        // 값 1: 고도, 값 2: 쌍선형 가중치 (dx,dy 무작위)
        std::vector<double> zval(N), wval(N);
        for (uint32_t k = 0; k < N; ++k) {
            zval[k] = elev(rng);
            double dx = frac(rng), dy = frac(rng);
            wval[k] = (1.0 - dx) * (1.0 - dy);       // w00 형태, 0~1
        }
        auto keys = cc->KeyGen();
        auto test = [&](const std::vector<double>& v, double& mabs, double& mrel, double& rrel) {
            // encode -> encrypt -> decrypt  (동형연산 0). '초반 근사(인코딩) + 암호화 잡음' 바닥.
            auto pt = cc->MakeCKKSPackedPlaintext(v);
            auto ct = cc->Encrypt(keys.publicKey, pt);
            Plaintext o; cc->Decrypt(keys.secretKey, ct, &o); o->SetLength(N);
            auto out = o->GetRealPackedValue();
            mabs = mrel = 0.0;
            long double sq = 0.0L;
            for (uint32_t k = 0; k < N; ++k) {
                q e = fabsq((q)out[k] - (q)v[k]);       // v[k] 자체가 double 이므로 이게 인코딩 오차
                double rel = (double)(e / fabsq((q)v[k] + (q)1e-30));
                mabs = std::max(mabs, (double)e);
                mrel = std::max(mrel, rel);
                sq += (long double)rel * (long double)rel;
            }
            rrel = std::sqrt((double)(sq / N));
        };
        double za, zr, zrms;  test(zval, za, zr, zrms);
        std::printf("    %5u | z(고도)  %.3e  | %.3e  | %.3e\n", N, za, zr, zrms);
        double wa, wr, wrms;  test(wval, wa, wr, wrms);
        std::printf("          | w(가중치) %.3e | %.3e  | %.3e\n", wa, wr, wrms);
    }
    std::printf("\n");
}

int main() {
    std::printf("== CKKS 인코딩(초반 근사) 오차만 — 암호화·연산 없음, 기준 __float128 ==\n\n");
    for (uint32_t b : {30u, 40u, 50u, 59u}) sweepN(b);
    std::printf("해석:\n");
    std::printf(" - 인코딩 오차 = 실수를 (값*2^b) 반올림해 정수화할 때의 양자화 오차.\n");
    std::printf(" - 이론 상한: 슬롯당 ~ (0.5 * sqrt(N_ring)) / 2^b.  b 가 클수록, N 이 작을수록 작다.\n");
    std::printf(" - 50비트에서 z(고도, 수백 m)의 상대 인코딩 오차가 이후 암호문 오차의 바닥이 된다.\n");
    return 0;
}

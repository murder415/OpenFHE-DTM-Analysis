// "10^-13 보다 커야 한다" 반론에 대한 엄격 재검증.
//
//  - 프로젝트 elevation() 회로를 그대로:  encrypt(zij) -> EvalMult(ct, pt_wij) -> EvalAdd x3 -> Decrypt
//  - 컨텍스트도 OpenFheBackend 와 동일: depth 3, ScalingModSize 50, COMPLEX, BatchSize 2048
//  - 어려운 입력: 한 셀 안에서 코너 고도가 크게 다르고 dx=dy=0.5 (보간 기여 최대)
//  - 기준선을 __float128 (113-bit) 로 잡아 double 반올림 의혹 제거
//  - OpenFHE 자체 오차 추정치 plaintext->GetLogError() 도 출력 (라이브러리가 스스로 보는 값)
//  - 스케일링 기법 FIXEDMANUAL(최악) 도 비교

#include "openfhe.h"
#include <quadmath.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

using namespace lbcrypto;

static void run(const std::string& tag, ScalingTechnique tech, uint32_t depth,
                uint32_t modsize, uint32_t batch) {
    // 한 셀: 코너가 평지~산지처럼 크게 다름
    const std::vector<double> z00v{ 12.0, 900.0,  30.0, 1113.0,  600.0 };
    const std::vector<double> z10v{ 480.0,  15.0, 780.0,  60.0,  210.0 };
    const std::vector<double> z01v{ 250.0, 640.0, 120.0, 300.0,  990.0 };
    const std::vector<double> z11v{ 1100.0, 70.0, 555.0, 820.0,   45.0 };
    const double dx = 0.5, dy = 0.5;                    // 보간 최대 기여
    const double w00 = (1-dx)*(1-dy), w10 = dx*(1-dy),
                 w01 = (1-dx)*dy,     w11 = dx*dy;      // 각 0.25
    const std::size_t n = z00v.size();

    // 기준선: __float128
    std::vector<__float128> refq(n);
    std::vector<double> refd(n);
    for (std::size_t i = 0; i < n; ++i) {
        __float128 r = (__float128)z00v[i]*(__float128)w00 + (__float128)z10v[i]*(__float128)w10
                     + (__float128)z01v[i]*(__float128)w01 + (__float128)z11v[i]*(__float128)w11;
        refq[i] = r;
        refd[i] = z00v[i]*w00 + z10v[i]*w10 + z01v[i]*w01 + z11v[i]*w11;
    }

    CCParams<CryptoContextCKKSRNS> p;
    p.SetMultiplicativeDepth(depth);
    p.SetScalingModSize(modsize);
    p.SetBatchSize(batch);
    p.SetCKKSDataType(COMPLEX);
    p.SetSecurityLevel(HEStd_128_classic);
    if (tech != INVALID_RS_TECHNIQUE) p.SetScalingTechnique(tech);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);
    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    auto enc = [&](const std::vector<double>& v) {
        return cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(v));
    };
    auto ct00 = enc(z00v), ct10 = enc(z10v), ct01 = enc(z01v), ct11 = enc(z11v);
    auto pw00 = cc->MakeCKKSPackedPlaintext(std::vector<double>(n, w00));
    auto pw10 = cc->MakeCKKSPackedPlaintext(std::vector<double>(n, w10));
    auto pw01 = cc->MakeCKKSPackedPlaintext(std::vector<double>(n, w01));
    auto pw11 = cc->MakeCKKSPackedPlaintext(std::vector<double>(n, w11));
    auto sum = cc->EvalAdd(cc->EvalAdd(cc->EvalMult(ct00, pw00), cc->EvalMult(ct10, pw10)),
                           cc->EvalAdd(cc->EvalMult(ct01, pw01), cc->EvalMult(ct11, pw11)));
    Plaintext out;
    cc->Decrypt(keys.secretKey, sum, &out);
    out->SetLength(n);
    auto got = out->GetRealPackedValue();

    double max_abs = 0, max_rel = 0, max_abs_vs_d = 0;
    for (std::size_t i = 0; i < n; ++i) {
        __float128 eq = fabsq((__float128)got[i] - refq[i]);
        double e = (double)eq;
        double rel = (double)(eq / fabsq(refq[i]));
        max_abs = std::max(max_abs, e);
        max_rel = std::max(max_rel, rel);
        max_abs_vs_d = std::max(max_abs_vs_d, std::fabs(got[i] - refd[i]));
    }
    // double 기준선 자체가 __float128 대비 얼마나 틀어졌나
    double base_roundoff = 0;
    for (std::size_t i = 0; i < n; ++i)
        base_roundoff = std::max(base_roundoff, (double)fabsq((__float128)refd[i] - refq[i]));

    std::printf("%-22s depth=%u mod=%u  ", tag.c_str(), depth, modsize);
    std::printf("MAX abs=%.3e m  rel=%.3e  | vs double=%.3e | double바닥=%.3e | OpenFHE GetLogError=%.1f bit (=2^-x -> %.2e)\n",
                max_abs, max_rel, max_abs_vs_d, base_roundoff,
                out->GetLogError(), std::pow(2.0, out->GetLogError()));
}

int main() {
    std::cout << "== elevation() 회로, 어려운 입력(코너 크게 다름, dx=dy=0.5), 값 ~수백 m ==\n";
    std::cout << "== 기준선: __float128 (113-bit).  '10^-13 보다 커야 한다' 검증 ==\n\n";

    run("project-equivalent", INVALID_RS_TECHNIQUE, 3, 50, 2048);   // 백엔드와 동일(기본 스케일링)
    run("FIXEDMANUAL(최악)",  FIXEDMANUAL,          3, 50, 2048);
    run("FIXEDAUTO",          FIXEDAUTO,            3, 50, 2048);
    run("FLEXIBLEAUTO",       FLEXIBLEAUTO,         3, 50, 2048);
    std::cout << "\n-- scaling modulus 를 낮추면 (이론상 오차 2^-b 로 증가) --\n";
    run("mod 40",             INVALID_RS_TECHNIQUE, 3, 40, 2048);
    run("mod 35",             INVALID_RS_TECHNIQUE, 3, 35, 2048);
    run("mod 30",             INVALID_RS_TECHNIQUE, 3, 30, 2048);

    std::cout << "\n판정: project-equivalent 의 rel 이 10^-13 근처이고, mod 를 낮추면 커지고,\n"
                 "      OpenFHE 자체 GetLogError 도 같은 값을 말하면 -> 측정이 맞다.\n";
    return 0;
}

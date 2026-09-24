// CKKS 정밀도 검증 프로브 — "암호문 보간 오차가 이렇게 작을 수 있나?" 에 답하기 위한 것.
//
// 검증하는 것:
//  (A) 회로가 진짜 암호문 상태인가:  Encrypt -> EvalMult(cipher,plaintext) -> EvalAdd -> Decrypt
//      (프로젝트 elevation() 경로와 동일한 depth-1 쌍선형 보간)
//  (B) 오차가 CKKS scaling factor 에 이론대로 반응하는가:
//      ScalingModSize 를 30~59 비트로 바꾸며 오차를 측정.
//      진짜 CKKS 라면 +10비트마다 오차가 약 2^10(~1000)배 줄어야 한다.
//      평문을 몰래 쓰고 있었다면 오차가 modulus 와 무관하게 일정하다.
//  (C) 값이 정말 암호문 안에 있는가(tamper 검사):
//      복호 전에 암호문에 +K 를 더하면 복호값도 정확히 +K 만큼 변해야 한다.
//  (D) double 기준선의 자체 반올림:  long double 기준과도 비교해
//      측정된 "오차" 가 CKKS 잡음인지 double 반올림인지 가른다.

#include "openfhe.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace lbcrypto;

int main() {
    // 실제 대만 타일에서 뽑은 것과 비슷한 값 범위(평지 ~ 산지 1100 m 대)
    const std::vector<double> corners = {12.531, 803.274, 45.118, 1113.902};
    const std::vector<double> weights = {0.20, 0.30, 0.35, 0.15};  // 합 = 1 (쌍선형 가중치)

    // 기준선 두 개: double / long double
    double ref_d = 0.0;
    long double ref_ld = 0.0L;
    for (int i = 0; i < 4; ++i) {
        ref_d  += corners[i] * weights[i];
        ref_ld += static_cast<long double>(corners[i]) * static_cast<long double>(weights[i]);
    }
    const long double baseline_roundoff = std::fabsl(static_cast<long double>(ref_d) - ref_ld);

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "interpolated value      = " << ref_d << " m\n";
    std::cout << "double vs long double   = " << (double)baseline_roundoff
              << "  (double 기준선 자체 반올림 바닥)\n\n";
    std::cout << "  bits |   abs_err (m) |  rel_err     | ratio vs prev | note\n";
    std::cout << "  -----+---------------+--------------+---------------+------\n";

    double prev_abs = 0.0;
    for (uint32_t bits : {30u, 35u, 40u, 45u, 50u, 59u}) {
        CCParams<CryptoContextCKKSRNS> p;
        p.SetMultiplicativeDepth(2);
        p.SetScalingModSize(bits);
        p.SetBatchSize(8);
        p.SetSecurityLevel(HEStd_128_classic);
        auto cc = GenCryptoContext(p);
        cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);

        // ---- 프로젝트 elevation() 과 동일한 depth-1 회로 ----
        Ciphertext<DCRTPoly> sum;
        for (int i = 0; i < 4; ++i) {
            auto pt_z = cc->MakeCKKSPackedPlaintext(std::vector<double>{corners[i]});
            auto pt_w = cc->MakeCKKSPackedPlaintext(std::vector<double>{weights[i]});
            auto ct   = cc->Encrypt(keys.publicKey, pt_z);         // 암호화
            auto term = cc->EvalMult(ct, pt_w);                    // 암호문 x 평문가중치
            sum = (i == 0) ? term : cc->EvalAdd(sum, term);        // 암호문 덧셈
        }

        Plaintext out;
        cc->Decrypt(keys.secretKey, sum, &out);
        out->SetLength(1);
        const double got = out->GetRealPackedValue()[0];
        const double abs_err = std::fabs(got - ref_d);
        const double rel_err = abs_err / std::fabs(ref_d);
        const double ratio   = (prev_abs > 0.0) ? prev_abs / abs_err : 0.0;

        std::cout << "   " << std::setw(3) << bits << " | " << std::setw(12) << abs_err
                  << "  | " << std::setw(11) << rel_err << "  | "
                  << std::setw(11) << (ratio > 0 ? std::to_string(ratio).substr(0, 9) : "   -   ")
                  << "   |";
        if (bits == 50) std::cout << " <- 논문 파라미터";
        std::cout << "\n";
        prev_abs = abs_err;

        // ---- (C) tamper 검사: 50비트에서만 ----
        if (bits == 50) {
            const double K = 137.0;
            auto pt_k = cc->MakeCKKSPackedPlaintext(std::vector<double>{K});
            auto ct_k = cc->EvalAdd(sum, pt_k);
            Plaintext out_k;
            cc->Decrypt(keys.secretKey, ct_k, &out_k);
            out_k->SetLength(1);
            const double got_k = out_k->GetRealPackedValue()[0];
            std::cout << "        tamper: decrypt(sum + " << K << ") - decrypt(sum) = "
                      << (got_k - got) << "  (기대: 1.3700e+02)\n";
        }
    }

    // ---- (E) 밀집 패킹 재현: 50비트 + 다중 slot + 큰 값 범위 ----
    //      experiments_validation 의 529점(MAX ~1.7e-11 m)이 왜 나오는지 확인.
    std::cout << "\n[E] 50비트, N개 값을 한 암호문에 패킹했을 때 MAX abs_err:\n";
    std::cout << "  N     | 값 범위        | MAX abs_err (m) | MAX rel_err\n";
    std::cout << "  ------+---------------+-----------------+-----------\n";
    for (uint32_t N : {1u, 16u, 64u, 256u, 529u, 2048u}) {
        CCParams<CryptoContextCKKSRNS> p;
        p.SetMultiplicativeDepth(2);
        p.SetScalingModSize(50);
        p.SetBatchSize(4096);
        p.SetSecurityLevel(HEStd_128_classic);
        auto cc = GenCryptoContext(p);
        cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);

        // N개 지점: 고도 5~1114 m, 가중치 0.1~0.4 를 순환
        std::vector<double> z0(N), z1(N), w0(N), w1(N), ref(N);
        for (uint32_t k = 0; k < N; ++k) {
            z0[k] = 5.0 + (1109.0 * k) / std::max(1u, N - 1);
            z1[k] = 1114.0 - (900.0 * k) / std::max(1u, N - 1);
            w0[k] = 0.35; w1[k] = 0.65;
            ref[k] = z0[k] * w0[k] + z1[k] * w1[k];
        }
        auto ct0 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(z0));
        auto ct1 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(z1));
        auto acc = cc->EvalAdd(cc->EvalMult(ct0, cc->MakeCKKSPackedPlaintext(w0)),
                               cc->EvalMult(ct1, cc->MakeCKKSPackedPlaintext(w1)));
        Plaintext out;
        cc->Decrypt(keys.secretKey, acc, &out);
        out->SetLength(N);
        auto got = out->GetRealPackedValue();
        double max_abs = 0.0, max_rel = 0.0;
        for (uint32_t k = 0; k < N; ++k) {
            double e = std::fabs(got[k] - ref[k]);
            max_abs = std::max(max_abs, e);
            max_rel = std::max(max_rel, e / std::fabs(ref[k]));
        }
        std::cout << "  " << std::setw(5) << N << " | 5 ~ 1114 m     | "
                  << std::setw(14) << max_abs << "  | " << std::setw(10) << max_rel << "\n";
    }

    std::cout << "\n판정:\n";
    std::cout << " - +10비트마다 abs_err 가 ~1000배 줄면  => 진짜 CKKS 잡음.\n";
    std::cout << " - modulus 를 바꿔도 abs_err 가 일정하면 => 평문 우회(가짜).\n";
    std::cout << " - 50비트 abs_err 가 'double vs long double' 값의 수 배~수십 배면\n";
    std::cout << "   측정 오차의 상당부분이 double 기준선 반올림. 상대오차로 보고할 것.\n";
    return 0;
}

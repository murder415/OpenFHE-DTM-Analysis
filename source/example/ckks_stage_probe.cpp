// elevation() 회로의 단계별 오차 분해.
//   각 노드에서 복호 -> __float128 정확값과 비교 -> 그 단계까지 누적된 절대/상대오차.
//   목적: "인코딩 / 암호화 / 곱셈(ct*pt) / 덧셈 / 최종 합"  중 어디서 오차가 생기는지.
//
//   회로 (프로젝트 elevation() 와 동일):
//     pt_zij  = Encode(zij)               [E1] 인코딩
//     ct_zij  = Encrypt(pt_zij)           [E2] 암호화(신선 잡음)
//     pt_wij  = Encode(wij)               [E1'] 가중치 인코딩
//     c_ij    = EvalMult(ct_zij, pt_wij)  [E3] ct*pt 곱 + rescale
//     L       = EvalAdd(c00, c10)         [E4] 덧셈
//     R       = EvalAdd(c01, c11)         [E4]
//     sum     = EvalAdd(L, R)             [E5] 덧셈
//     out     = Decrypt(sum)              [E6] 복호

#include "openfhe.h"
#include <quadmath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <iostream>
#include <vector>

using namespace lbcrypto;
using q = __float128;
using Clock = std::chrono::steady_clock;

static double elapsedMs(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

static void printTiming(const char* name, std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    const double median = values.size() % 2
        ? values[values.size() / 2]
        : (values[values.size() / 2 - 1] + values[values.size() / 2]) / 2.0;
    std::printf("%-37s %9.4f  %9.4f  %9.4f  %9.4f\n",
                name, mean, median, values.front(), values.back());
}

static double maxAbsRel(const std::vector<double>& got, const std::vector<q>& ref,
                        double& rel_out) {
    double a = 0, r = 0;
    for (size_t i = 0; i < got.size(); ++i) {
        q e = fabsq((q)got[i] - ref[i]);
        a = std::max(a, (double)e);
        r = std::max(r, (double)(e / fabsq(ref[i])));
    }
    rel_out = r;
    return a;
}

int main() {
    // 한 셀, 코너 크게 다름, dx=dy=0.5
    const std::vector<double> z00{ 12.0, 900.0,  30.0, 1113.0, 600.0 };
    const std::vector<double> z10{ 480.0,  15.0, 780.0,  60.0,  210.0 };
    const std::vector<double> z01{ 250.0, 640.0, 120.0, 300.0,  990.0 };
    const std::vector<double> z11{ 1100.0, 70.0, 555.0, 820.0,   45.0 };
    const double dx = 0.5, dy = 0.5;
    const double w00 = (1-dx)*(1-dy), w10 = dx*(1-dy), w01 = (1-dx)*dy, w11 = dx*dy;
    const size_t n = z00.size();

    // __float128 정확 기준값(단계별)
    std::vector<q> ref_z00(n), ref_c00(n), ref_c10(n), ref_c01(n), ref_c11(n),
                   ref_L(n), ref_R(n), ref_sum(n);
    for (size_t i = 0; i < n; ++i) {
        ref_z00[i] = (q)z00[i];
        ref_c00[i] = (q)z00[i]*(q)w00; ref_c10[i] = (q)z10[i]*(q)w10;
        ref_c01[i] = (q)z01[i]*(q)w01; ref_c11[i] = (q)z11[i]*(q)w11;
        ref_L[i]   = ref_c00[i] + ref_c10[i];
        ref_R[i]   = ref_c01[i] + ref_c11[i];
        ref_sum[i] = ref_L[i] + ref_R[i];
    }

    CCParams<CryptoContextCKKSRNS> p;
    p.SetMultiplicativeDepth(3);
    p.SetScalingModSize(50);
    p.SetBatchSize(2048);
    p.SetCKKSDataType(COMPLEX);
    p.SetSecurityLevel(HEStd_128_classic);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);
    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    auto dec = [&](auto ct) {
        Plaintext o; cc->Decrypt(keys.secretKey, ct, &o); o->SetLength(n);
        return o->GetRealPackedValue();
    };
    double rel;
    std::printf("stage                                    abs(m)        rel\n");
    std::printf("--------------------------------------------------------------\n");

    // [E1] 인코딩만 (encode -> 곧바로 decode, 암호화 없음)
    {
        auto pt = cc->MakeCKKSPackedPlaintext(z00);
        pt->SetLength(n);
        auto v = pt->GetRealPackedValue();
        double a = maxAbsRel(v, ref_z00, rel);
        std::printf("[E1] Encode(z) 만 (decode 왕복)          %.3e   %.3e\n", a, rel);
        // 가중치 인코딩 오차
        std::vector<q> refw{ (q)w00 };
        auto pw = cc->MakeCKKSPackedPlaintext(std::vector<double>{w00});
        pw->SetLength(1);
        double rr; double aw = maxAbsRel(pw->GetRealPackedValue(), refw, rr);
        std::printf("     Encode(w=%.2f) 오차                  %.3e   %.3e\n", w00, aw, rr);
    }

    // [E2] Encrypt->Decrypt 왕복 (연산 없음)
    auto ct00 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(z00));
    auto ct10 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(z10));
    auto ct01 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(z01));
    auto ct11 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(z11));
    {
        double a = maxAbsRel(dec(ct00), ref_z00, rel);
        std::printf("[E2] Encrypt->Decrypt (연산 0)           %.3e   %.3e\n", a, rel);
    }

    auto pw00 = cc->MakeCKKSPackedPlaintext(std::vector<double>(n, w00));
    auto pw10 = cc->MakeCKKSPackedPlaintext(std::vector<double>(n, w10));
    auto pw01 = cc->MakeCKKSPackedPlaintext(std::vector<double>(n, w01));
    auto pw11 = cc->MakeCKKSPackedPlaintext(std::vector<double>(n, w11));

    // [E3] EvalMult(ct, pt) + rescale — 항별
    auto c00 = cc->EvalMult(ct00, pw00);
    auto c10 = cc->EvalMult(ct10, pw10);
    auto c01 = cc->EvalMult(ct01, pw01);
    auto c11 = cc->EvalMult(ct11, pw11);
    { double a=maxAbsRel(dec(c00),ref_c00,rel); std::printf("[E3] c00 = EvalMult(ct00, w00)           %.3e   %.3e\n",a,rel);}
    { double a=maxAbsRel(dec(c10),ref_c10,rel); std::printf("[E3] c10 = EvalMult(ct10, w10)           %.3e   %.3e\n",a,rel);}
    { double a=maxAbsRel(dec(c01),ref_c01,rel); std::printf("[E3] c01 = EvalMult(ct01, w01)           %.3e   %.3e\n",a,rel);}
    { double a=maxAbsRel(dec(c11),ref_c11,rel); std::printf("[E3] c11 = EvalMult(ct11, w11)           %.3e   %.3e\n",a,rel);}

    // 항별 오차 합(코히런트 상한 예측)
    {
        auto d00=dec(c00), d10=dec(c10), d01=dec(c01), d11=dec(c11);
        double sum_of_term_abs = 0;
        for (size_t i=0;i<n;++i){
            double s = (double)(fabsq((q)d00[i]-ref_c00[i]) + fabsq((q)d10[i]-ref_c10[i])
                              + fabsq((q)d01[i]-ref_c01[i]) + fabsq((q)d11[i]-ref_c11[i]));
            sum_of_term_abs = std::max(sum_of_term_abs, s);
        }
        std::printf("     (항별 오차 절대 합, 최악 상한 예측)   %.3e\n", sum_of_term_abs);
    }

    // [E4] 덧셈
    auto L = cc->EvalAdd(c00, c10);
    auto R = cc->EvalAdd(c01, c11);
    { double a=maxAbsRel(dec(L),ref_L,rel); std::printf("[E4] L = EvalAdd(c00, c10)               %.3e   %.3e\n",a,rel);}
    { double a=maxAbsRel(dec(R),ref_R,rel); std::printf("[E4] R = EvalAdd(c01, c11)               %.3e   %.3e\n",a,rel);}

    // [E5] 최종 덧셈
    auto sum = cc->EvalAdd(L, R);
    { double a=maxAbsRel(dec(sum),ref_sum,rel);
      std::printf("[E5] sum = EvalAdd(L, R)  == 최종 결과   %.3e   %.3e\n",a,rel);}

    // [E6] 복호는 위 dec() 에 포함. Decrypt 자체 추가 오차는 무시 수준(같은 키, 반올림).

    std::cout << "\n해석:\n"
              << " - E1(인코딩) + E2(암호화) 가 바닥.  E3(곱+rescale) 에서 조금 늘고,\n"
              << "   E4/E5(덧셈) 는 거의 안 늘림(덧셈은 잡음도 더해질 뿐 증폭 없음).\n"
              << " - 최종 오차 ≈ 항별 오차들의 부분적 합(√n~n).  '합산이라 커진다'는 게 이 부분.\n"
              << " - 전부 scaling factor 2^-50 근처에서 시작하므로 최종도 ~10^-13 상대.\n";

    // 529개 보간값을 한 암호문에 넣었을 때의 실제 회로 단계별 시간.
    // context/key 생성은 제외하고, 5회 준비 실행 뒤 50회 반복한다.
    constexpr size_t timing_n = 529;
    constexpr int warmups = 5;
    constexpr int repeats = 50;
    std::vector<double> tz00(timing_n), tz10(timing_n), tz01(timing_n), tz11(timing_n);
    std::vector<double> tw00(timing_n, w00), tw10(timing_n, w10);
    std::vector<double> tw01(timing_n, w01), tw11(timing_n, w11);
    for (size_t i = 0; i < timing_n; ++i) {
        tz00[i] = z00[i % n]; tz10[i] = z10[i % n];
        tz01[i] = z01[i % n]; tz11[i] = z11[i % n];
    }

    std::vector<double> encode_z_ms, encode_w_ms, encrypt_ms, multiply_ms;
    std::vector<double> add_ms, decrypt_ms, total_ms;
    double checksum = 0.0;
    for (int run = -warmups; run < repeats; ++run) {
        const auto total_begin = Clock::now();

        auto begin = Clock::now();
        auto pz00 = cc->MakeCKKSPackedPlaintext(tz00);
        auto pz10 = cc->MakeCKKSPackedPlaintext(tz10);
        auto pz01 = cc->MakeCKKSPackedPlaintext(tz01);
        auto pz11 = cc->MakeCKKSPackedPlaintext(tz11);
        auto end = Clock::now();
        const double this_encode_z = elapsedMs(begin, end);

        begin = Clock::now();
        auto tct00 = cc->Encrypt(keys.publicKey, pz00);
        auto tct10 = cc->Encrypt(keys.publicKey, pz10);
        auto tct01 = cc->Encrypt(keys.publicKey, pz01);
        auto tct11 = cc->Encrypt(keys.publicKey, pz11);
        end = Clock::now();
        const double this_encrypt = elapsedMs(begin, end);

        begin = Clock::now();
        auto tpw00 = cc->MakeCKKSPackedPlaintext(tw00);
        auto tpw10 = cc->MakeCKKSPackedPlaintext(tw10);
        auto tpw01 = cc->MakeCKKSPackedPlaintext(tw01);
        auto tpw11 = cc->MakeCKKSPackedPlaintext(tw11);
        end = Clock::now();
        const double this_encode_w = elapsedMs(begin, end);

        begin = Clock::now();
        auto tc00 = cc->EvalMult(tct00, tpw00);
        auto tc10 = cc->EvalMult(tct10, tpw10);
        auto tc01 = cc->EvalMult(tct01, tpw01);
        auto tc11 = cc->EvalMult(tct11, tpw11);
        end = Clock::now();
        const double this_multiply = elapsedMs(begin, end);

        begin = Clock::now();
        auto tleft = cc->EvalAdd(tc00, tc10);
        auto tright = cc->EvalAdd(tc01, tc11);
        auto tsum = cc->EvalAdd(tleft, tright);
        end = Clock::now();
        const double this_add = elapsedMs(begin, end);

        begin = Clock::now();
        Plaintext tout;
        cc->Decrypt(keys.secretKey, tsum, &tout);
        tout->SetLength(timing_n);
        const auto decoded = tout->GetRealPackedValue();
        end = Clock::now();
        const double this_decrypt = elapsedMs(begin, end);
        const double this_total = elapsedMs(total_begin, end);
        checksum += decoded.front();

        if (run >= 0) {
            encode_z_ms.push_back(this_encode_z);
            encrypt_ms.push_back(this_encrypt);
            encode_w_ms.push_back(this_encode_w);
            multiply_ms.push_back(this_multiply);
            add_ms.push_back(this_add);
            decrypt_ms.push_back(this_decrypt);
            total_ms.push_back(this_total);
        }
    }

    std::printf("\n529-slot timing (warm-up=%d, repeats=%d, key generation excluded)\n",
                warmups, repeats);
    std::printf("stage                                  mean(ms) median(ms)    min(ms)    max(ms)\n");
    std::printf("-------------------------------------------------------------------------------\n");
    printTiming("Encode elevation vectors x4", encode_z_ms);
    printTiming("Encrypt x4", encrypt_ms);
    printTiming("Encode weight vectors x4", encode_w_ms);
    printTiming("EvalMult(ciphertext, plaintext) x4", multiply_ms);
    printTiming("EvalAdd x3", add_ms);
    printTiming("Decrypt + decode x1", decrypt_ms);
    printTiming("End-to-end circuit", total_ms);
    std::printf("checksum=%.6f\n", checksum);
    return 0;
}

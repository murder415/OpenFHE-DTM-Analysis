#pragma once

#include "openfhe_dtm/HeBackend.hpp"

#ifdef OPENFHE_DTM_ENABLE_OPENFHE
#include "pke/openfhe.h"
#endif

namespace openfhe_dtm {

class OpenFheBackend final : public IHeBackend {
public:
    explicit OpenFheBackend(std::size_t slots = 16,
                            std::filesystem::path state_directory = {},
                            bool enable_bootstrap = false);

    std::size_t slotCount() const override;
    PlainVector encode(const std::vector<double>& values, int level = 0) const override;
    PlainVector encodeComplex(const std::vector<std::complex<double>>& values,
                              int level = 0) const override;
    CipherVector encrypt(const PlainVector& plain) const override;
    PlainVector decrypt(const CipherVector& cipher) const override;
    std::vector<std::complex<double>> decryptComplex(const CipherVector& cipher) const override;

    CipherVector add(const CipherVector& a, const CipherVector& b) const override;
    CipherVector sub(const CipherVector& a, const CipherVector& b) const override;
    CipherVector mul(const CipherVector& a, const CipherVector& b) const override;
    CipherVector mulPlain(const CipherVector& a, const PlainVector& b) const override;
    CipherVector subPlain(const CipherVector& a, const PlainVector& b) const override;
    CipherVector mulScalar(const CipherVector& a, std::complex<double> scalar) const override;
    CipherVector conjugate(const CipherVector& a) const override;
    CipherVector rotate(const CipherVector& a, int offset) const override;
    CipherVector sumSlots(const CipherVector& a) const override;
    bool supportsBootstrap() const override;
    CipherVector bootstrap(const CipherVector& a) const override;
    CipherVector evaluateFunction(const CipherVector& a, double lower, double upper,
                                  std::size_t degree,
                                  double (*function)(double)) const override;
    CipherVector evaluateChebyshevSeries(const CipherVector& a,
                                         const std::vector<double>& coefficients,
                                         double lower, double upper) const override;
    void saveCipher(const CipherVector& cipher, const std::filesystem::path& path) const override;
    CipherVector loadCipher(const std::filesystem::path& path) const override;

private:
    struct Impl;
    std::shared_ptr<Impl> impl_;
};

} // namespace openfhe_dtm

#pragma once

#include <memory>
#include <filesystem>
#include <complex>
#include <string>
#include <vector>

namespace openfhe_dtm {

struct PlainVector {
    std::vector<double> values;
    int level = 0;
    std::shared_ptr<void> native;
};

struct CipherVector {
    std::shared_ptr<void> native;
    std::size_t slot_count = 0;
    int level = 0;
    std::string backend;
};

class IHeBackend {
public:
    virtual ~IHeBackend() = default;

    virtual std::size_t slotCount() const = 0;
    virtual PlainVector encode(const std::vector<double>& values, int level = 0) const = 0;
    virtual PlainVector encodeComplex(const std::vector<std::complex<double>>& values,
                                      int level = 0) const = 0;
    virtual CipherVector encrypt(const PlainVector& plain) const = 0;
    virtual PlainVector decrypt(const CipherVector& cipher) const = 0;
    virtual std::vector<std::complex<double>> decryptComplex(const CipherVector& cipher) const = 0;

    virtual CipherVector add(const CipherVector& a, const CipherVector& b) const = 0;
    virtual CipherVector sub(const CipherVector& a, const CipherVector& b) const = 0;
    virtual CipherVector mul(const CipherVector& a, const CipherVector& b) const = 0;
    virtual CipherVector mulPlain(const CipherVector& a, const PlainVector& b) const = 0;
    virtual CipherVector subPlain(const CipherVector& a, const PlainVector& b) const = 0;
    // Plaintext minus ciphertext, with an encrypted result.
    virtual CipherVector plainSub(const PlainVector& a, const CipherVector& b) const = 0;
    virtual CipherVector mulScalar(const CipherVector& a,
                                   std::complex<double> scalar) const = 0;
    virtual CipherVector conjugate(const CipherVector& a) const = 0;
    virtual CipherVector rotate(const CipherVector& a, int offset) const = 0;
    virtual CipherVector sumSlots(const CipherVector& a) const = 0;
    virtual bool supportsBootstrap() const = 0;
    virtual CipherVector bootstrap(const CipherVector& a) const = 0;
    virtual CipherVector evaluateFunction(const CipherVector& a, double lower, double upper,
                                          std::size_t degree,
                                          double (*function)(double)) const = 0;
    virtual CipherVector evaluateChebyshevSeries(const CipherVector& a,
                                                 const std::vector<double>& coefficients,
                                                 double lower, double upper) const = 0;
    virtual void saveCipher(const CipherVector& cipher,
                            const std::filesystem::path& path) const = 0;
    virtual CipherVector loadCipher(const std::filesystem::path& path) const = 0;
};

} // namespace openfhe_dtm

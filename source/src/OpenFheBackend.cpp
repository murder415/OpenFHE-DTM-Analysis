#include "openfhe_dtm/OpenFheBackend.hpp"

#ifdef OPENFHE_DTM_ENABLE_OPENFHE
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"
#include "scheme/ckksrns/ckksrns-fhe.h"
#include "version.h"
#endif

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace openfhe_dtm {

#ifdef OPENFHE_DTM_ENABLE_OPENFHE
namespace {

using CryptoContextT = lbcrypto::CryptoContext<lbcrypto::DCRTPoly>;
using CiphertextT = lbcrypto::Ciphertext<lbcrypto::DCRTPoly>;
using KeyPairT = lbcrypto::KeyPair<lbcrypto::DCRTPoly>;

std::shared_ptr<lbcrypto::Plaintext> asPlain(const PlainVector& plain) {
    if (!plain.native) {
        throw std::runtime_error("PlainVector does not contain an OpenFHE plaintext");
    }
    return std::static_pointer_cast<lbcrypto::Plaintext>(plain.native);
}

std::shared_ptr<CiphertextT> asCipher(const CipherVector& cipher) {
    if (cipher.backend != "openfhe" || !cipher.native) {
        throw std::runtime_error("CipherVector does not contain an OpenFHE ciphertext");
    }
    return std::static_pointer_cast<CiphertextT>(cipher.native);
}

CipherVector wrapCipher(CiphertextT ciphertext, std::size_t slots, int level) {
    return CipherVector{
        std::make_shared<CiphertextT>(std::move(ciphertext)),
        slots,
        level,
        "openfhe",
    };
}

std::string readTextFile(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("missing OpenFHE key manifest: " + path.string());
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

void requireManifestField(const std::string& manifest, const std::string& field) {
    if (manifest.find(field) == std::string::npos)
        throw std::runtime_error("OpenFHE key manifest mismatch: expected " + field);
}

std::string manifestJson(std::size_t slots, uint32_t depth, bool bootstrap,
                         const std::string& key_tag) {
    std::ostringstream out;
    out << "{\n"
        << "  \"format_version\": 1,\n"
        << "  \"library\": \"OpenFHE\",\n"
        << "  \"openfhe_version\": \"" << GetOPENFHEVersion() << "\",\n"
        << "  \"scheme\": \"CKKS\",\n"
        << "  \"batch_size\": " << slots << ",\n"
        << "  \"multiplicative_depth\": " << depth << ",\n"
        << "  \"scaling_modulus_bits\": 50,\n"
        << "  \"ckks_data_type\": \"COMPLEX\",\n"
        << "  \"bootstrap_enabled\": " << (bootstrap ? "true" : "false") << ",\n"
        << "  \"bootstrap_level_budget\": [3, 3],\n"
        << "  \"key_set_id\": \"" << key_tag << "\",\n"
        << "  \"rotation_indices\": [1, -1, 2, -2, 4, -4, 8, -8";
    if (slots >= 16384) out << ", 16, -16, 32, -32, 64, -64, 128, -128";
    out << "]\n}\n";
    return out.str();
}

void writeManifestAtomically(const std::filesystem::path& path,
                             const std::string& contents) {
    const auto temporary = path.string() + ".tmp";
    {
        std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        if (!output || !(output << contents) || !output.flush())
            throw std::runtime_error("failed to write OpenFHE key manifest");
    }
    std::filesystem::rename(temporary, path);
}

} // namespace
#endif

struct OpenFheBackend::Impl {
    std::size_t slots = 16;
    bool bootstrap_enabled = false;
    std::vector<uint32_t> bootstrap_level_budget{3, 3};

#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    CryptoContextT cc;
    KeyPairT keys;
#endif
};

OpenFheBackend::OpenFheBackend(std::size_t slots, std::filesystem::path state_directory,
                               bool enable_bootstrap)
    : impl_(std::make_shared<Impl>()) {
    impl_->slots = slots;
    impl_->bootstrap_enabled = enable_bootstrap;

#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    const uint32_t levels_after_bootstrap = slots >= 16384 ? 16 : 3;
    const auto secret_key_distribution = lbcrypto::UNIFORM_TERNARY;
    uint32_t multiplicative_depth = levels_after_bootstrap;
    if (impl_->bootstrap_enabled) {
        multiplicative_depth += lbcrypto::FHECKKSRNS::GetBootstrapDepth(
            impl_->bootstrap_level_budget, secret_key_distribution);
    }
    const auto context_path = state_directory / "cryptocontext.bin";
    if (!state_directory.empty() && std::filesystem::exists(context_path)) {
        const auto manifest = readTextFile(state_directory / "manifest.json");
        requireManifestField(manifest, "\"format_version\": 1");
        requireManifestField(manifest, "\"library\": \"OpenFHE\"");
        requireManifestField(manifest, "\"openfhe_version\": \"" +
                                           GetOPENFHEVersion() + "\"");
        requireManifestField(manifest, "\"scheme\": \"CKKS\"");
        requireManifestField(manifest, "\"batch_size\": " +
                                           std::to_string(slots) + ",");
        requireManifestField(manifest, "\"multiplicative_depth\": " +
                                           std::to_string(multiplicative_depth) + ",");
        requireManifestField(manifest, "\"scaling_modulus_bits\": 50");
        requireManifestField(manifest, "\"ckks_data_type\": \"COMPLEX\"");
        requireManifestField(manifest, std::string("\"bootstrap_enabled\": ") +
                                           (impl_->bootstrap_enabled ? "true," : "false,"));
        requireManifestField(manifest, "\"bootstrap_level_budget\": [3, 3]");
        const std::string expected_rotations = slots >= 16384
            ? "\"rotation_indices\": [1, -1, 2, -2, 4, -4, 8, -8, 16, -16, 32, -32, 64, -64, 128, -128]"
            : "\"rotation_indices\": [1, -1, 2, -2, 4, -4, 8, -8]";
        requireManifestField(manifest, expected_rotations);
        if (!lbcrypto::Serial::DeserializeFromFile(context_path.string(), impl_->cc,
                                                   lbcrypto::SerType::BINARY) ||
            !lbcrypto::Serial::DeserializeFromFile((state_directory / "public.key").string(),
                                                   impl_->keys.publicKey,
                                                   lbcrypto::SerType::BINARY) ||
            !lbcrypto::Serial::DeserializeFromFile((state_directory / "secret.key").string(),
                                                   impl_->keys.secretKey,
                                                   lbcrypto::SerType::BINARY)) {
            throw std::runtime_error("failed to restore OpenFHE context/key set");
        }
        requireManifestField(manifest, "\"key_set_id\": \"" +
                                           impl_->keys.publicKey->GetKeyTag() + "\"");
        if (impl_->bootstrap_enabled) {
            impl_->cc->Enable(lbcrypto::FHE);
            impl_->cc->EvalBootstrapSetup(impl_->bootstrap_level_budget, {0, 0},
                                          static_cast<uint32_t>(impl_->slots));
        }
        // OpenFHE keeps evaluation keys in process-wide maps. Remove an older
        // copy of this same key set before restoring it again in one process.
        CryptoContextT::element_type::ClearEvalMultKeys(
            impl_->keys.publicKey->GetKeyTag());
        CryptoContextT::element_type::ClearEvalAutomorphismKeys(
            impl_->keys.publicKey->GetKeyTag());
        std::ifstream mult_keys(state_directory / "eval-mult.keys", std::ios::binary);
        std::ifstream rotate_keys(state_directory / "eval-rotate.keys", std::ios::binary);
        if (!mult_keys || !impl_->cc->DeserializeEvalMultKey(mult_keys, lbcrypto::SerType::BINARY) ||
            !rotate_keys || !impl_->cc->DeserializeEvalAutomorphismKey(
                                rotate_keys, lbcrypto::SerType::BINARY)) {
            throw std::runtime_error("failed to restore OpenFHE evaluation keys");
        }
        return;
    }

    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> parameters;
    if (impl_->bootstrap_enabled) {
        parameters.SetSecretKeyDist(secret_key_distribution);
        parameters.SetScalingTechnique(lbcrypto::FLEXIBLEAUTO);
        parameters.SetKeySwitchTechnique(lbcrypto::HYBRID);
    }
    parameters.SetMultiplicativeDepth(multiplicative_depth);
    parameters.SetScalingModSize(50);
    parameters.SetBatchSize(static_cast<uint32_t>(slots));
    parameters.SetCKKSDataType(lbcrypto::COMPLEX);

    impl_->cc = lbcrypto::GenCryptoContext(parameters);
    impl_->cc->Enable(lbcrypto::PKE);
    impl_->cc->Enable(lbcrypto::KEYSWITCH);
    impl_->cc->Enable(lbcrypto::LEVELEDSHE);
    impl_->cc->Enable(lbcrypto::ADVANCEDSHE);
    if (impl_->bootstrap_enabled) {
        impl_->cc->Enable(lbcrypto::FHE);
        impl_->cc->EvalBootstrapSetup(impl_->bootstrap_level_budget, {0, 0},
                                      static_cast<uint32_t>(impl_->slots));
    }

    impl_->keys = impl_->cc->KeyGen();
    impl_->cc->EvalMultKeyGen(impl_->keys.secretKey);
    if (impl_->bootstrap_enabled) {
        impl_->cc->EvalBootstrapKeyGen(impl_->keys.secretKey,
                                       static_cast<uint32_t>(impl_->slots));
    }
    std::vector<int32_t> rotations{1, -1, 2, -2, 4, -4, 8, -8};
    if (slots >= 16384) {
        rotations.insert(rotations.end(), {16, -16, 32, -32, 64, -64, 128, -128});
    }
    impl_->cc->EvalRotateKeyGen(impl_->keys.secretKey, rotations);
    // 전체 slot 합(EvalSum)에 필요한 자동형태(automorphism) 키. 순 토공량 등 격자 전역
    // 집계를 암호문 상태로 수행할 때 사용한다. batch 크기에 로그 비례하는 소수의 키만 생성된다.
    impl_->cc->EvalSumKeyGen(impl_->keys.secretKey);
    const uint32_t conjugation_index = impl_->cc->GetCyclotomicOrder() - 1;
    auto conjugation_key = lbcrypto::FHECKKSRNS::ConjugateKeyGen(impl_->keys.secretKey);
    lbcrypto::CryptoContextImpl<lbcrypto::DCRTPoly>::GetEvalAutomorphismKeyMap(
        impl_->keys.secretKey->GetKeyTag())[conjugation_index] = std::move(conjugation_key);

    if (!state_directory.empty()) {
        std::filesystem::create_directories(state_directory);
        if (!lbcrypto::Serial::SerializeToFile(context_path.string(), impl_->cc,
                                               lbcrypto::SerType::BINARY) ||
            !lbcrypto::Serial::SerializeToFile((state_directory / "public.key").string(),
                                               impl_->keys.publicKey,
                                               lbcrypto::SerType::BINARY) ||
            !lbcrypto::Serial::SerializeToFile((state_directory / "secret.key").string(),
                                               impl_->keys.secretKey,
                                               lbcrypto::SerType::BINARY)) {
            throw std::runtime_error("failed to persist OpenFHE context/key set");
        }
        std::ofstream mult_keys(state_directory / "eval-mult.keys", std::ios::binary);
        std::ofstream rotate_keys(state_directory / "eval-rotate.keys", std::ios::binary);
        if (!mult_keys || !impl_->cc->SerializeEvalMultKey(mult_keys, lbcrypto::SerType::BINARY) ||
            !rotate_keys || !impl_->cc->SerializeEvalAutomorphismKey(
                                rotate_keys, lbcrypto::SerType::BINARY)) {
            throw std::runtime_error("failed to persist OpenFHE evaluation keys");
        }
        mult_keys.close();
        rotate_keys.close();
        writeManifestAtomically(state_directory / "manifest.json",
                                manifestJson(slots, multiplicative_depth,
                                             impl_->bootstrap_enabled,
                                             impl_->keys.publicKey->GetKeyTag()));
    }
#endif
}

std::size_t OpenFheBackend::slotCount() const { return impl_->slots; }

PlainVector OpenFheBackend::encode(const std::vector<double>& values, int level) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    std::vector<double> padded = values;
    if (padded.size() < impl_->slots) {
        padded.resize(impl_->slots, 0.0);
    }
    auto native = impl_->cc->MakeCKKSPackedPlaintext(padded);
    return PlainVector{std::move(padded), level, std::make_shared<lbcrypto::Plaintext>(native)};
#else
    (void)values;
    (void)level;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

PlainVector OpenFheBackend::encodeComplex(const std::vector<std::complex<double>>& values,
                                          int level) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    std::vector<std::complex<double>> padded = values;
    if (padded.size() < impl_->slots) padded.resize(impl_->slots, {0.0, 0.0});
    auto native = impl_->cc->MakeCKKSPackedPlaintext(padded);
    std::vector<double> real(padded.size());
    std::transform(padded.begin(), padded.end(), real.begin(),
                   [](const auto& z) { return z.real(); });
    return PlainVector{std::move(real), level,
                       std::make_shared<lbcrypto::Plaintext>(native)};
#else
    (void)values; (void)level;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::encrypt(const PlainVector& plain) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    auto ciphertext = impl_->cc->Encrypt(impl_->keys.publicKey, *asPlain(plain));
    return wrapCipher(std::move(ciphertext), impl_->slots, plain.level);
#else
    (void)plain;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

PlainVector OpenFheBackend::decrypt(const CipherVector& cipher) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    lbcrypto::Plaintext result;
    impl_->cc->Decrypt(impl_->keys.secretKey, *asCipher(cipher), &result);
    result->SetLength(impl_->slots);
    std::vector<double> values = result->GetRealPackedValue();
    return PlainVector{std::move(values), cipher.level, std::make_shared<lbcrypto::Plaintext>(result)};
#else
    (void)cipher;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

std::vector<std::complex<double>> OpenFheBackend::decryptComplex(const CipherVector& cipher) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    lbcrypto::Plaintext result;
    impl_->cc->Decrypt(impl_->keys.secretKey, *asCipher(cipher), &result);
    result->SetLength(impl_->slots);
    return result->GetCKKSPackedValue();
#else
    (void)cipher;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::add(const CipherVector& a, const CipherVector& b) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    return wrapCipher(impl_->cc->EvalAdd(*asCipher(a), *asCipher(b)), impl_->slots,
                      std::min(a.level, b.level));
#else
    (void)a;
    (void)b;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::sub(const CipherVector& a, const CipherVector& b) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    return wrapCipher(impl_->cc->EvalSub(*asCipher(a), *asCipher(b)), impl_->slots,
                      std::min(a.level, b.level));
#else
    (void)a;
    (void)b;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::mul(const CipherVector& a, const CipherVector& b) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    return wrapCipher(impl_->cc->EvalMult(*asCipher(a), *asCipher(b)), impl_->slots,
                      std::min(a.level, b.level) - 1);
#else
    (void)a; (void)b;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::mulPlain(const CipherVector& a, const PlainVector& b) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    return wrapCipher(impl_->cc->EvalMult(*asCipher(a), *asPlain(b)), impl_->slots, a.level - 1);
#else
    (void)a;
    (void)b;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::subPlain(const CipherVector& a, const PlainVector& b) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    return wrapCipher(impl_->cc->EvalSub(*asCipher(a), *asPlain(b)), impl_->slots, a.level);
#else
    (void)a; (void)b;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::mulScalar(const CipherVector& a,
                                       std::complex<double> scalar) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    return wrapCipher(impl_->cc->EvalMult(*asCipher(a), scalar), impl_->slots, a.level - 1);
#else
    (void)a; (void)scalar;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::conjugate(const CipherVector& a) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    auto native = *asCipher(a);
    const uint32_t index = impl_->cc->GetCyclotomicOrder() - 1;
    const auto& keys = CryptoContextT::element_type::GetEvalAutomorphismKeyMap(
        native->GetKeyTag());
    (void)index;
    return wrapCipher(lbcrypto::FHECKKSRNS::Conjugate(native, keys), impl_->slots, a.level);
#else
    (void)a;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::rotate(const CipherVector& a, int offset) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    return wrapCipher(impl_->cc->EvalRotate(*asCipher(a), offset), impl_->slots, a.level);
#else
    (void)a;
    (void)offset;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::sumSlots(const CipherVector& a) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    return wrapCipher(impl_->cc->EvalSum(*asCipher(a), static_cast<uint32_t>(impl_->slots)),
                      impl_->slots, a.level);
#else
    (void)a;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

bool OpenFheBackend::supportsBootstrap() const {
    return impl_->bootstrap_enabled;
}

CipherVector OpenFheBackend::bootstrap(const CipherVector& a) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    if (!impl_->bootstrap_enabled) {
        throw std::runtime_error("OpenFHE bootstrap was not enabled for this backend");
    }
    auto native = *asCipher(a);
    auto refreshed = impl_->cc->EvalBootstrap(native);
    return wrapCipher(std::move(refreshed), impl_->slots, 16);
#else
    (void)a;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::evaluateFunction(const CipherVector& a, double lower,
                                              double upper, std::size_t degree,
                                              double (*function)(double)) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    auto native = *asCipher(a);
    auto result = impl_->cc->EvalChebyshevFunction(function, native, lower, upper,
                                                   static_cast<uint32_t>(degree));
    return wrapCipher(std::move(result), impl_->slots, a.level - 5);
#else
    (void)a; (void)lower; (void)upper; (void)degree; (void)function;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::evaluateChebyshevSeries(
    const CipherVector& a, const std::vector<double>& coefficients,
    double lower, double upper) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    auto native = *asCipher(a);
    auto result = impl_->cc->EvalChebyshevSeries(native, coefficients, lower, upper);
    return wrapCipher(std::move(result), impl_->slots, a.level - 5);
#else
    (void)a; (void)coefficients; (void)lower; (void)upper;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

void OpenFheBackend::saveCipher(const CipherVector& cipher,
                                const std::filesystem::path& path) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    std::filesystem::create_directories(path.parent_path());
    if (!lbcrypto::Serial::SerializeToFile(path.string(), *asCipher(cipher),
                                           lbcrypto::SerType::BINARY)) {
        throw std::runtime_error("failed to serialize OpenFHE ciphertext");
    }
#else
    (void)cipher; (void)path;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

CipherVector OpenFheBackend::loadCipher(const std::filesystem::path& path) const {
#ifdef OPENFHE_DTM_ENABLE_OPENFHE
    CiphertextT ciphertext;
    if (!lbcrypto::Serial::DeserializeFromFile(path.string(), ciphertext,
                                               lbcrypto::SerType::BINARY)) {
        throw std::runtime_error("failed to deserialize OpenFHE ciphertext");
    }
    return wrapCipher(std::move(ciphertext), impl_->slots, 0);
#else
    (void)path;
    throw std::runtime_error("OpenFHE support is disabled");
#endif
}

} // namespace openfhe_dtm

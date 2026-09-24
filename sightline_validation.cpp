#include "openfhe_dtm/Analysis.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/TdbDataset.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace openfhe_dtm;
namespace fs = std::filesystem;

namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double earth_radius_m = 6371008.8;

void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

std::string jsonString(const std::string& s) {
    std::string out = "\"";
    for (const unsigned char c : s) {
        if (c == '\\' || c == '"') { out += '\\'; out += static_cast<char>(c); }
        else if (c == '\n') out += "\\n";
        else if (c == '\r') out += "\\r";
        else if (c == '\t') out += "\\t";
        else if (c < 32) throw std::runtime_error("Unexpected control character in JSON");
        else out += static_cast<char>(c);
    }
    return out + '"';
}

// This wrapper observes only calls made by the production analysis. The
// independent baseline reads the original provider directly, after HE finishes.
class GuardedProvider final : public IDtmProvider {
public:
    explicit GuardedProvider(const IDtmProvider& provider) : provider_(provider) {}
    const DtmMetadata& metadata() const override { return provider_.metadata(); }
    InterpolationStencil interpolationStencil(const GeoPoint& p) const override {
        ++stencil_calls;
        return provider_.interpolationStencil(p);
    }
    double elevationAt(const GeoPoint&) const override {
        throw std::runtime_error("Plain elevationAt was used inside the HE sightline path");
    }
    DtmTile loadTile(int x, int y) const override { return provider_.loadTile(x, y); }
    mutable std::size_t stencil_calls = 0;
private:
    const IDtmProvider& provider_;
};

struct Event {
    std::string operation;
    std::size_t output = 0, input_a = 0, input_b = 0;
};

// Enforce the stronger contract of the new two-phase API: all encrypted
// computation finishes before the first decrypt. Only further final-batch
// decrypt calls are permitted afterwards; a re-encryption or even encode fails.
class AuditedBackend final : public IHeBackend {
public:
    explicit AuditedBackend(const IHeBackend& backend) : backend_(backend) {}
    void reset() const {
        events.clear(); decrypted_values.clear(); nodes_.clear(); lifetimes_.clear();
        decrypt_started_ = false; encrypt_calls = decrypt_calls = 0; next_id_ = 1;
    }
    std::size_t slotCount() const override { return backend_.slotCount(); }
    PlainVector encode(const std::vector<double>& v, int l = 0) const override {
        before("encode"); events.push_back({"encode"}); return backend_.encode(v, l);
    }
    PlainVector encodeComplex(const std::vector<std::complex<double>>& v, int l = 0) const override {
        before("encodeComplex"); events.push_back({"encodeComplex"}); return backend_.encodeComplex(v, l);
    }
    CipherVector encrypt(const PlainVector& p) const override {
        before("encrypt"); ++encrypt_calls; return record("encrypt", backend_.encrypt(p));
    }
    PlainVector decrypt(const CipherVector& c) const override {
        const auto input = id(c);
        require(nodes_.at(c.native.get()).second == "sub",
                "Final sightline decrypt did not receive a ciphertext height difference");
        decrypt_started_ = true; ++decrypt_calls;
        events.push_back({"decrypt", 0, input});
        auto p = backend_.decrypt(c); decrypted_values.push_back(p.values); return p;
    }
    std::vector<std::complex<double>> decryptComplex(const CipherVector&) const override {
        throw std::runtime_error("Unexpected complex decryption in the sightline path");
    }
    CipherVector add(const CipherVector& a, const CipherVector& b) const override {
        before("add"); return record("add", backend_.add(a,b), id(a),id(b));
    }
    CipherVector sub(const CipherVector& a, const CipherVector& b) const override {
        before("sub"); return record("sub", backend_.sub(a,b), id(a),id(b));
    }
    CipherVector mul(const CipherVector& a, const CipherVector& b) const override {
        before("mul"); return record("mul", backend_.mul(a,b), id(a),id(b));
    }
    CipherVector mulPlain(const CipherVector& a, const PlainVector& b) const override {
        before("mulPlain"); return record("mulPlain", backend_.mulPlain(a,b), id(a));
    }
    CipherVector subPlain(const CipherVector& a, const PlainVector& b) const override {
        before("subPlain"); return record("subPlain", backend_.subPlain(a,b), id(a));
    }
    CipherVector mulScalar(const CipherVector& a, std::complex<double> b) const override {
        before("mulScalar"); return record("mulScalar", backend_.mulScalar(a,b), id(a));
    }
    CipherVector conjugate(const CipherVector& a) const override {
        before("conjugate"); return record("conjugate", backend_.conjugate(a), id(a));
    }
    CipherVector rotate(const CipherVector& a, int n) const override {
        before("rotate"); return record("rotate", backend_.rotate(a,n), id(a));
    }
    CipherVector sumSlots(const CipherVector& a) const override {
        before("sumSlots"); return record("sumSlots", backend_.sumSlots(a), id(a));
    }
    bool supportsBootstrap() const override { return backend_.supportsBootstrap(); }
    CipherVector bootstrap(const CipherVector& a) const override {
        before("bootstrap"); return record("bootstrap", backend_.bootstrap(a), id(a));
    }
    CipherVector evaluateFunction(const CipherVector& a, double lo, double hi,
                                   std::size_t n, double (*f)(double)) const override {
        before("evaluateFunction");
        return record("evaluateFunction", backend_.evaluateFunction(a,lo,hi,n,f), id(a));
    }
    CipherVector evaluateChebyshevSeries(const CipherVector& a, const std::vector<double>& c,
                                          double lo, double hi) const override {
        before("evaluateChebyshevSeries");
        return record("evaluateChebyshevSeries", backend_.evaluateChebyshevSeries(a,c,lo,hi), id(a));
    }
    void saveCipher(const CipherVector& a, const fs::path& p) const override {
        before("saveCipher"); events.push_back({"saveCipher", 0,id(a)}); backend_.saveCipher(a,p);
    }
    CipherVector loadCipher(const fs::path& p) const override {
        before("loadCipher"); return record("loadCipher", backend_.loadCipher(p));
    }
    mutable std::vector<Event> events;
    mutable std::vector<std::vector<double>> decrypted_values;
    mutable std::size_t encrypt_calls = 0, decrypt_calls = 0;
private:
    void before(const char* operation) const {
        require(!decrypt_started_, std::string("HE operation after decryption: ") + operation);
    }
    std::size_t id(const CipherVector& c) const {
        auto it = nodes_.find(c.native.get());
        require(it != nodes_.end(), "Ciphertext bypassed the audited backend");
        return it->second.first;
    }
    CipherVector record(const char* operation, CipherVector c,
                        std::size_t a = 0, std::size_t b = 0) const {
        require(c.native != nullptr, "Backend returned a null ciphertext");
        const auto output = next_id_++;
        // Keep every native handle alive, so pointer identity cannot be reused.
        lifetimes_.push_back(c.native);
        nodes_[c.native.get()] = {output,operation};
        events.push_back({operation,output,a,b}); return c;
    }
    const IHeBackend& backend_;
    mutable bool decrypt_started_ = false;
    mutable std::size_t next_id_ = 1;
    mutable std::map<const void*,std::pair<std::size_t,std::string>> nodes_;
    mutable std::vector<std::shared_ptr<void>> lifetimes_;
};

// Higher precision independent geodesic distance, rather than reusing the
// production helper. The plaintext baseline never calls analysis.elevation().
double independentDistance(const GeoPoint& a, const GeoPoint& b) {
    const long double r = 6371008.8L;
    const long double rad = std::acos(-1.0L) / 180.0L;
    const long double x = (static_cast<long double>(b.lat)-a.lat)*rad/2;
    const long double y = (static_cast<long double>(b.lon)-a.lon)*rad/2;
    const long double hav = std::sin(x)*std::sin(x) +
        std::cos(a.lat*rad)*std::cos(b.lat*rad)*std::sin(y)*std::sin(y);
    return static_cast<double>(2*r*std::asin(std::sqrt(std::min(1.0L,hav))));
}

struct CaseResult {
    std::string name;
    std::size_t samples = 0, batches = 0, encryptions = 0, decryptions = 0;
    std::size_t stencil_calls = 0, near_zero = 0, near_zero_sign_mismatches = 0;
    double angle = 0, max_error = 0, rmse = 0, mae = 0, max_raw_error = 0;
    double max_distance_error = 0, max_analytic_error = 0, elapsed_seconds = 0;
    double first_height = 0, last_height = 0;
};

CaseResult checkCase(const std::string& name, const IDtmProvider& provider,
                     const SightSurfaceQuery& query, AuditedBackend& he, const fs::path& out,
                     std::size_t expected_samples,
                     const std::function<double(const GeoPoint&)>& analytic = {}) {
    he.reset(); GuardedProvider guarded(provider); TerrainAnalysis analysis(guarded,he);
    const auto start = std::chrono::steady_clock::now();
    const auto result = analysis.sightSurface(query);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    const auto n = result.points.size(), slots = he.slotCount();
    require(n == expected_samples, name + ": unexpected sample count");
    require(result.distances_m.size() == n && result.buildable_height_m.size() == n,
            name + ": output sizes differ");
    require(result.ground_m.empty() && result.regulation_plane_m.empty(),
            name + ": decrypted intermediate ground or plane was returned");
    const auto batches = (n+slots-1)/slots;
    require(he.decrypt_calls == batches, name + ": expected exactly one terminal decrypt per batch");
    require(!he.events.empty() && he.events.back().operation == "decrypt",
            name + ": final backend call is not decrypt");
    // before() already rejects any HE or encoding call after first decrypt.
    bool seen_decrypt = false;
    for (const auto& e : he.events) {
        if (e.operation == "decrypt") seen_decrypt = true;
        else require(!seen_decrypt, name + ": a backend operation followed decrypt");
    }
    std::ofstream trace(out/(name+"_trace.csv"));
    require(static_cast<bool>(trace), "Cannot open trace CSV");
    trace << "event,operation,output_cipher,input_a,input_b\n";
    for (std::size_t i=0;i<he.events.size();++i) {
        const auto& e=he.events[i];
        trace<<i<<','<<e.operation<<','<<e.output<<','<<e.input_a<<','<<e.input_b<<'\n';
    }
    std::ofstream csv(out/(name+"_samples.csv"));
    require(static_cast<bool>(csv), "Cannot open samples CSV");
    csv << std::setprecision(std::numeric_limits<double>::max_digits10);
    csv << "sample,lon,lat,independent_distance_m,returned_distance_m,plain_ground_m,plain_base_m,plain_plane_m,plain_difference_m,decoded_difference_m,plain_height_m,he_height_m,abs_error_m,near_zero\n";
    CaseResult summary; summary.name=name; summary.samples=n; summary.batches=batches;
    summary.encryptions=he.encrypt_calls; summary.decryptions=he.decrypt_calls;
    summary.angle=query.angle_degrees; summary.stencil_calls=guarded.stencil_calls;
    summary.elapsed_seconds=seconds;
    const double base=provider.elevationAt(query.base);
    const double tangent=std::tan(query.angle_degrees*pi/180.0);
    double distance=0, squared=0;
    // 1e-7 m is only a reporting band for numerically delicate sign tests.
    // It is NOT a clipping threshold: production still computes max(0,h).
    constexpr double near_zero_band_m=1e-7;
    for (std::size_t i=0;i<n;++i) {
        if(i) distance+=independentDistance(result.points[i-1],result.points[i]);
        const double z=provider.elevationAt(result.points[i]);
        const double plane=base+distance*tangent;
        const double diff=plane-z;
        const double expected=std::max(0.0,diff);
        const auto& decoded=he.decrypted_values.at(i/slots);
        require(decoded.size()>i%slots, "Decoded batch shorter than active samples");
        const double raw=decoded[i%slots];
        const double actual=result.buildable_height_m[i];
        require(std::isfinite(raw) && std::isfinite(actual), name+": non-finite output");
        require(actual == std::max(0.0,raw), name+": output does not match final max(0,decoded)");
        const double error=std::abs(actual-expected);
        const bool near=std::abs(diff)<=near_zero_band_m;
        if(near) {
            ++summary.near_zero;
            if((diff>0)!=(raw>0)) ++summary.near_zero_sign_mismatches;
        }
        // Away from zero, a sign disagreement is a functional failure.
        if(!near) require((diff>0)==(raw>0), name+": wrong sign away from zero");
        summary.max_error=std::max(summary.max_error,error);
        summary.max_raw_error=std::max(summary.max_raw_error,std::abs(raw-diff));
        summary.max_distance_error=std::max(summary.max_distance_error,std::abs(distance-result.distances_m[i]));
        if(analytic) summary.max_analytic_error=std::max(summary.max_analytic_error,std::abs(z-analytic(result.points[i])));
        summary.mae+=error; squared+=error*error;
        csv<<i<<','<<result.points[i].lon<<','<<result.points[i].lat<<','<<distance<<','
           <<result.distances_m[i]<<','<<z<<','<<base<<','<<plane<<','<<diff<<','<<raw<<','
           <<expected<<','<<actual<<','<<error<<','<<(near?1:0)<<'\n';
    }
    summary.mae/=static_cast<double>(n); summary.rmse=std::sqrt(squared/static_cast<double>(n));
    summary.first_height=result.buildable_height_m.front();
    summary.last_height=result.buildable_height_m.back();
    require(summary.max_error<1e-6 && summary.max_raw_error<1e-6,
            name+": error exceeds predeclared 1e-6 m validation tolerance");
    require(summary.max_distance_error<1e-7, name+": independent distance mismatch");
    if(analytic) require(summary.max_analytic_error<1e-9, name+": plaintext raster is not the analytic plane");
    return summary;
}

RasterDtmProvider planeProvider(double slope_east=0, double slope_north=0) {
    RasterGrid grid; grid.origin={0,0}; grid.width=128; grid.height=128; grid.cell_size_m=5;
    for(std::size_t y=0;y<grid.height;++y) for(std::size_t x=0;x<grid.width;++x)
        grid.elevations_m.push_back(100+slope_east*(5*x)+slope_north*(5*y));
    DtmMetadata metadata; metadata.name="Synthetic analytic plane"; metadata.source="test fixture";
    metadata.crs="local equirectangular over geographic coordinates"; metadata.resolution_m=5;
    return RasterDtmProvider(metadata,std::move(grid));
}

std::function<double(const GeoPoint&)> planeFormula(double sx=0,double sy=0) {
    return [sx,sy](const GeoPoint& p) {
        return 100+sx*p.lon*pi/180.0*earth_radius_m+sy*p.lat*pi/180.0*earth_radius_m;
    };
}

void writeCaseJson(std::ostream& out,const CaseResult& r) {
    out<<"{\"name\":"<<jsonString(r.name)<<",\"passed\":true"
       <<",\"samples\":"<<r.samples<<",\"batches\":"<<r.batches
       <<",\"angle_degrees\":"<<r.angle<<",\"encrypt_calls\":"<<r.encryptions
       <<",\"decrypt_calls\":"<<r.decryptions<<",\"stencil_calls\":"<<r.stencil_calls
       <<",\"all_he_before_first_decrypt\":true,\"intermediate_plain_lookup\":false"
       <<",\"max_abs_error_m\":"<<r.max_error<<",\"mae_m\":"<<r.mae
       <<",\"rmse_m\":"<<r.rmse<<",\"max_raw_difference_error_m\":"<<r.max_raw_error
       <<",\"max_distance_error_m\":"<<r.max_distance_error
       <<",\"max_analytic_fixture_error_m\":"<<r.max_analytic_error
       <<",\"near_zero_samples\":"<<r.near_zero
       <<",\"near_zero_sign_mismatches\":"<<r.near_zero_sign_mismatches
       <<",\"first_height_m\":"<<r.first_height<<",\"last_height_m\":"<<r.last_height
       <<",\"elapsed_seconds\":"<<r.elapsed_seconds<<'}';
}
} // namespace

int main(int argc,char** argv) {
    if(argc<3 || argc>4) {
        std::cerr<<"usage: sightline_validation <taiwan.tdb> <output_dir> [runs=4]\n"; return 2;
    }
    try {
        const fs::path output=argv[2]; fs::create_directories(output);
        const int runs=argc==4?std::stoi(argv[3]):4;
        require(runs>0 && runs<=100,"runs must be in [1,100]");
        const auto taiwan=loadValidationTdbDataset(argv[1]);
        const double span=std::max(2.0,std::min(taiwan.widthMeters(),taiwan.heightMeters())*0.55);
        const double x0=span*0.12,y0=span*0.12;
        const SightSurfaceQuery taiwan_query{taiwan.pointFromMeters(x0,y0),
            taiwan.pointFromMeters(x0+span,y0+span*0.61),12.0,span,std::max(2.0,span/9.0)};
        auto flat=planeProvider(); auto uphill=planeProvider(0.3,0.07); auto downhill=planeProvider(-0.1,0);
        std::cout<<std::setprecision(std::numeric_limits<double>::max_digits10);
        std::vector<std::vector<CaseResult>> all_results;
        for(int run=1;run<=runs;++run) {
            // Empty state directory means a fresh key pair for each repeat.
            OpenFheBackend backend(16); AuditedBackend he(backend);
            const fs::path run_dir=output/("run_"+std::to_string(run)); fs::create_directories(run_dir);
            std::vector<CaseResult> results;
            results.push_back(checkCase("taiwan_demo_10",taiwan,taiwan_query,he,run_dir,10));
            const auto exercise=[&](const std::string& name,const RasterDtmProvider& p,
                                     double angle,double length,double interval,std::size_t count,
                                     double sx=0,double sy=0) {
                const SightSurfaceQuery q{p.pointFromMeters(80,80),p.pointFromMeters(180,80),angle,length,interval};
                results.push_back(checkCase(name,p,q,he,run_dir,count,planeFormula(sx,sy)));
            };
            exercise("flat_positive_angle",flat,12,90,10,10);
            exercise("flat_zero_angle",flat,0,90,10,10);
            exercise("flat_negative_angle",flat,-12,90,10,10);
            exercise("uphill_above_plane",uphill,12,90,10,10,0.3,0.07);
            exercise("downhill_zero_angle",downhill,0,90,10,10,-0.1,0);
            exercise("exact_slot_boundary",flat,12,150,10,16);
            exercise("one_past_slots",flat,12,160,10,17);
            exercise("three_batches_tail_9",flat,12,400,10,41);
            exercise("minimum_two_samples",flat,12,0.1,10,2);
            exercise("near_zero_positive",flat,1e-10,90,10,10);
            exercise("near_zero_negative",flat,-1e-10,90,10,10);
            // Sightline geometry deliberately has at least two samples. Test
            // n=1 independently at the encrypted interpolation API boundary.
            he.reset(); GuardedProvider single_guard(flat); TerrainAnalysis single_analysis(single_guard,he);
            const GeoPoint single_point=flat.pointFromMeters(83.2,87.4);
            const CipherVector single=single_analysis.elevationEncrypted(Route{single_point});
            require(he.decrypt_calls==0 && he.encrypt_calls==4,
                    "Single-point encrypted interpolation performed a decrypt or changed corner packing");
            // This explicit test readout is outside the production API, after
            // all production calls. It is not an intermediate decryption.
            const PlainVector single_readout=backend.decrypt(single);
            require(!single_readout.values.empty(),"Single-point test readout is empty");
            const double single_error=std::abs(single_readout.values.front()-flat.elevationAt(single_point));
            require(single_error<1e-6,"Single-point interpolation error exceeds tolerance");
            std::ofstream json(run_dir/"results.json"); require(static_cast<bool>(json),"Cannot open run JSON");
            json<<std::setprecision(std::numeric_limits<double>::max_digits10);
            json<<"{\"run\":"<<run<<",\"fresh_keys\":true,\"slots\":16,\"cases\":[";
            for(std::size_t i=0;i<results.size();++i) { if(i)json<<',';writeCaseJson(json,results[i]); }
            json<<"],\"single_point_encrypted_interpolation\":{\"samples\":1,\"api_decrypt_calls\":0,"
                <<"\"external_validation_readout_decrypt_calls\":1,\"abs_error_m\":"<<single_error<<"}}\n";
            std::cout<<"{\"run\":"<<run<<",\"taiwan\":";writeCaseJson(std::cout,results.front());
            std::cout<<",\"all_cases_passed\":"<<results.size()<<"}\n";
            all_results.push_back(std::move(results));
        }
        std::ofstream summary(output/"summary.json"); require(static_cast<bool>(summary),"Cannot open summary JSON");
        summary<<std::setprecision(std::numeric_limits<double>::max_digits10);
        summary<<"{\"passed\":true,\"runs\":"<<runs
            <<",\"numeric_tolerance_m\":1e-6,\"near_zero_reporting_band_m\":1e-7"
            <<",\"near_zero_policy\":\"Report sign changes within the band; do not clamp or mask them beyond max(0,decoded).\""
            <<",\"intermediate_decrypts\":0,\"he_operations_after_first_decrypt\":0"
            <<",\"plaintext_baseline\":\"Original DEM bilinear interpolation plus independently computed geodesic cumulative distance and input angle, then max(0,z_base+d*tan(angle)-z).\""
            <<",\"scope\":\"Final-only sightline path; does not validate terrain occlusion, encrypted source storage, distributed key separation, or field survey accuracy.\""
            <<",\"minimum_sightline_samples\":2,\"results\":[";
        bool first=true;
        for(std::size_t r=0;r<all_results.size();++r) for(const auto& c:all_results[r]) {
            if (!first) summary << ',';
            first = false;
            summary<<"{\"run\":"<<(r+1)<<",\"result\":";writeCaseJson(summary,c);summary<<'}';
        }
        summary<<"]}\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"sightline validation failed: "<<e.what()<<'\n';return 1;
    }
}

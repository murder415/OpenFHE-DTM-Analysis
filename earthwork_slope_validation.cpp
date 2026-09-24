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
        throw std::runtime_error("Plain elevationAt was used inside the HE earthwork slope path");
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
        if(live_trace_.is_open()) live_trace_.close();
        live_trace_.clear();case_name_.clear();failure_path_.clear();
        events.clear(); decrypted_values.clear(); decrypted_ids.clear(); encrypted_inputs.clear(); nodes_.clear(); lifetimes_.clear();
        decrypt_started_ = false; encrypt_calls = decrypt_calls = 0; next_id_ = 1;
    }
    void startCase(const std::string& name,const fs::path& output) const {
        case_name_=name;failure_path_=output/(name+"_encode_failure.json");
        live_trace_.open(output/(name+"_live_trace.csv"));
        require(static_cast<bool>(live_trace_),"Cannot create live trace CSV");
        live_trace_<<"event,operation,values,maximum_abs,minimum_nonzero_abs,output_cipher,input_a,input_b\n";
        live_trace_<<std::setprecision(std::numeric_limits<double>::max_digits10);
    }
    std::size_t slotCount() const override { return backend_.slotCount(); }
    PlainVector encode(const std::vector<double>& v, int l = 0) const override {
        before("encode"); events.push_back({"encode"});
        double maximum=0,minimum=std::numeric_limits<double>::infinity();
        for(double value:v) {maximum=std::max(maximum,std::abs(value));if(value!=0)minimum=std::min(minimum,std::abs(value));}
        if(live_trace_.is_open()) {
            live_trace_<<(events.size()-1)<<",encode,"<<v.size()<<','<<maximum<<','<<minimum<<",0,0,0\n";
            live_trace_.flush();
        }
        try {return backend_.encode(v,l);}
        catch(const std::exception& error) {
            if(!failure_path_.empty()) {
                std::ofstream failure(failure_path_);
                failure<<std::setprecision(std::numeric_limits<double>::max_digits10)
                       <<"{\"case\":"<<jsonString(case_name_)<<",\"slots\":"<<slotCount()
                       <<",\"event\":"<<(events.size()-1)<<",\"level_argument\":"<<l
                       <<",\"error\":"<<jsonString(error.what())<<",\"maximum_abs\":"<<maximum
                       <<",\"values\":[";
                for(std::size_t i=0;i<v.size();++i) {
                    if(i) failure<<',';
                    failure<<v[i];
                }
                failure<<"]}\n";
            }
            throw std::runtime_error(case_name_+": encode failed at event "+std::to_string(events.size()-1)+
                                     " for "+std::to_string(v.size())+" values: "+error.what());
        }
    }
    PlainVector encodeComplex(const std::vector<std::complex<double>>& v, int l = 0) const override {
        before("encodeComplex"); events.push_back({"encodeComplex"}); return backend_.encodeComplex(v, l);
    }
    CipherVector encrypt(const PlainVector& p) const override {
        before("encrypt"); ++encrypt_calls;
        auto c = record("encrypt", backend_.encrypt(p));
        encrypted_inputs[id(c)] = p.values; return c;
    }
    PlainVector decrypt(const CipherVector& c) const override {
        const auto input = id(c);
        require(nodes_.at(c.native.get()).second == "mulPlain",
                "Final earthwork slope decrypt did not receive a area-weighted ciphertext difference");
        decrypt_started_ = true; ++decrypt_calls;
        events.push_back({"decrypt", 0, input}); decrypted_ids.push_back(input);
        if(live_trace_.is_open()) {
            live_trace_<<(events.size()-1)<<",decrypt,0,0,0,0,"<<input<<",0\n";live_trace_.flush();
        }
        auto p = backend_.decrypt(c); decrypted_values.push_back(p.values); return p;
    }
    std::vector<std::complex<double>> decryptComplex(const CipherVector&) const override {
        throw std::runtime_error("Unexpected complex decryption in the earthwork slope path");
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
    std::size_t cipherId(const CipherVector& c) const { return id(c); }
    mutable std::map<std::size_t, std::vector<double>> encrypted_inputs;
    mutable std::vector<std::size_t> decrypted_ids;
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
        // Historical edges already store numeric IDs. An allocation may reuse
        // an expired handle's address, in which case this map entry receives a
        // fresh ID. Live candidate/input handles cannot reuse that address.
        // Do not retain hundreds of large intermediates solely for tracing.
        lifetimes_.push_back(c.native);
        nodes_[c.native.get()] = {output,operation};
        events.push_back({operation,output,a,b});
        if(live_trace_.is_open()) {
            live_trace_<<(events.size()-1)<<','<<operation<<",0,0,0,"<<output<<','<<a<<','<<b<<'\n';live_trace_.flush();
        }
        return c;
    }
    const IHeBackend& backend_;
    mutable bool decrypt_started_ = false;
    mutable std::size_t next_id_ = 1;
    mutable std::map<const void*,std::pair<std::size_t,std::string>> nodes_;
    mutable std::vector<std::weak_ptr<void>> lifetimes_;
    mutable std::ofstream live_trace_;
    mutable fs::path failure_path_;
    mutable std::string case_name_;
};

// The oracle does not call production geometry helpers, encrypted interpolation,
// or finalization. Use extended precision coordinates and original DEM lookups.
struct XY { long double x = 0, y = 0; };
XY local(const GeoPoint& origin, const GeoPoint& p) {
    const long double rad = std::acos(-1.0L) / 180;
    return {(static_cast<long double>(p.lon)-origin.lon)*rad*6371008.8L*std::cos(origin.lat*rad),
            (static_cast<long double>(p.lat)-origin.lat)*rad*6371008.8L};
}
GeoPoint geographic(const GeoPoint& origin, XY p) {
    const long double rad = std::acos(-1.0L) / 180;
    return {static_cast<double>(origin.lon+p.x/(6371008.8L*std::cos(origin.lat*rad)*rad)),
            static_cast<double>(origin.lat+p.y/(6371008.8L*rad))};
}
long double separation(XY a, XY b) { return std::hypot(a.x-b.x,a.y-b.y); }
bool inside(const std::vector<XY>& ring, XY p) {
    // Winding-number test, independent of the production ray-crossing helper.
    int winding = 0;
    for (std::size_t i=0;i<ring.size();++i) {
        const XY a=ring[i],b=ring[(i+1)%ring.size()];
        const long double cross=(b.x-a.x)*(p.y-a.y)-(b.y-a.y)*(p.x-a.x);
        if (a.y<=p.y && b.y>p.y && cross>0) ++winding;
        if (a.y>p.y && b.y<=p.y && cross<0) --winding;
    }
    return winding != 0;
}
struct ReferenceGeometry {
    GeoPoint origin;
    std::size_t columns=0,rows=0;
    long double min_x=0,min_y=0,dx=0,dy=0,area=0;
    std::vector<XY> nodes,boundary;
    std::vector<long double> weights;
    std::vector<std::array<std::size_t,4>> cells;
};
ReferenceGeometry referenceGeometry(const EarthworkSlopeQuery& q) {
    ReferenceGeometry g; g.origin=q.polygon.front();
    Route points=q.polygon;
    if (points.size()>1 && points.front().lon==points.back().lon && points.front().lat==points.back().lat)
        points.pop_back();
    std::vector<XY> ring;
    for (const auto& p:points) ring.push_back(local(g.origin,p));
    long double max_x=ring.front().x,max_y=ring.front().y;
    g.min_x=max_x;g.min_y=max_y;
    for (XY p:ring) {
        g.min_x=std::min(g.min_x,p.x);g.min_y=std::min(g.min_y,p.y);
        max_x=std::max(max_x,p.x);max_y=std::max(max_y,p.y);
    }
    g.columns=static_cast<std::size_t>(std::ceil((max_x-g.min_x)/q.sample_interval_m));
    g.rows=static_cast<std::size_t>(std::ceil((max_y-g.min_y)/q.sample_interval_m));
    g.dx=(max_x-g.min_x)/g.columns;g.dy=(max_y-g.min_y)/g.rows;
    using Index=std::pair<std::size_t,std::size_t>;
    std::map<Index,long double> weights;
    std::vector<std::array<Index,4>> cells;
    for (std::size_t y=0;y<g.rows;++y) for (std::size_t x=0;x<g.columns;++x) {
        if (!inside(ring,{g.min_x+(x+0.5L)*g.dx,g.min_y+(y+0.5L)*g.dy})) continue;
        std::array<Index,4> vertices{{{x,y},{x+1,y},{x,y+1},{x+1,y+1}}};
        for (const auto& vertex:vertices) weights[vertex]+=g.dx*g.dy/4;
        cells.push_back(vertices);g.area+=g.dx*g.dy;
    }
    std::map<Index,std::size_t> indices;
    for (const auto& item:weights) {
        indices[item.first]=g.nodes.size();
        g.nodes.push_back({g.min_x+item.first.first*g.dx,g.min_y+item.first.second*g.dy});
        g.weights.push_back(item.second);
    }
    for (const auto& cell:cells)
        g.cells.push_back({indices.at(cell[0]),indices.at(cell[1]),indices.at(cell[2]),indices.at(cell[3])});
    const double interval=q.boundary_interval_m==0?q.sample_interval_m:q.boundary_interval_m;
    for (std::size_t i=0;i<ring.size();++i) {
        const XY a=ring[i],b=ring[(i+1)%ring.size()];
        const auto count=static_cast<std::size_t>(std::ceil(separation(a,b)/interval));
        for (std::size_t j=0;j<count;++j) {
            const long double t=static_cast<long double>(j)/count;
            g.boundary.push_back({a.x+(b.x-a.x)*t,a.y+(b.y-a.y)*t});
        }
    }
    return g;
}

std::vector<std::size_t> compareGeometry(const EarthworkSlopeGeometry& actual,
                                        const ReferenceGeometry& expected) {
    constexpr long double position_tolerance_m=1e-6L;
    require(actual.grid_columns==expected.columns && actual.grid_rows==expected.rows,
            "Grid shape differs from independent ceil(bbox/interval) reference");
    require(actual.points.size()==expected.nodes.size() &&
            actual.point_weights_m2.size()==expected.weights.size(),"FourCorners node count mismatch");
    require(actual.boundary_points.size()==expected.boundary.size(),"Boundary sampling count mismatch");
    require(std::abs(actual.grid_dx_m-expected.dx)<position_tolerance_m &&
            std::abs(actual.grid_dy_m-expected.dy)<position_tolerance_m,"Grid spacing mismatch");
    require(std::abs(actual.integration_area_m2-expected.area)<1e-6L,"Quadrature area mismatch");
    std::vector<std::size_t> mapping;std::vector<bool> used(expected.nodes.size(),false);
    for (std::size_t i=0;i<actual.points.size();++i) {
        const XY p=local(expected.origin,actual.points[i]);
        std::size_t best=0;long double nearest=std::numeric_limits<long double>::infinity();
        for (std::size_t j=0;j<expected.nodes.size();++j) {
            const auto distance=separation(p,expected.nodes[j]);
            if(distance<nearest) {nearest=distance;best=j;}
        }
        require(nearest<position_tolerance_m && !used[best],"Grid node missing, duplicated or displaced");
        used[best]=true;mapping.push_back(best);
        require(std::abs(actual.point_weights_m2[i]-expected.weights[best])<1e-6L,
                "Shared-node area/4 accumulation mismatch");
    }
    using Cell=std::array<std::size_t,4>;
    std::vector<Cell> actual_cells,expected_cells=expected.cells;
    for (const auto& cell:actual.cells) {
        Cell mapped;
        for(std::size_t j=0;j<4;++j) {require(cell[j]<mapping.size(),"Invalid cell node index");mapped[j]=mapping[cell[j]];}
        std::sort(mapped.begin(),mapped.end());actual_cells.push_back(mapped);
    }
    for(auto& cell:expected_cells)std::sort(cell.begin(),cell.end());
    std::sort(actual_cells.begin(),actual_cells.end());std::sort(expected_cells.begin(),expected_cells.end());
    require(actual_cells==expected_cells,"Center-in-polygon selected cells differ from independent winding test");
    for(std::size_t i=0;i<actual.boundary_points.size();++i)
        require(separation(local(expected.origin,actual.boundary_points[i]),expected.boundary[i])<position_tolerance_m,
                "Boundary point order or edge subdivision mismatch");
    return mapping;
}

struct Oracle {
    std::vector<double> ground,boundary_ground,upper,lower,delta;
    std::vector<bool> infeasible;
    std::size_t invalid_points=0,invalid_boundary=0,invalid_cells=0;
    double max_gap=0,cut=0,fill=0,net=0,moved=0,cut_area=0,fill_area=0,unchanged_area=0;
    bool valid=true;
};
Oracle referenceModel(const IDtmProvider& provider,const EarthworkSlopeQuery& q,const ReferenceGeometry& g) {
    Oracle r;
    const long double slope=std::tan(static_cast<long double>(q.slope_angle_degrees)*std::acos(-1.0L)/180);
    for (XY b:g.boundary) r.boundary_ground.push_back(provider.elevationAt(geographic(g.origin,b)));
    std::vector<XY> targets=g.nodes;targets.insert(targets.end(),g.boundary.begin(),g.boundary.end());
    for (std::size_t i=0;i<targets.size();++i) {
        const XY p=targets[i];
        const double z=provider.elevationAt(geographic(g.origin,p));
        long double upper=std::numeric_limits<long double>::infinity(),lower=-upper;
        for (std::size_t b=0;b<g.boundary.size();++b) {
            const long double rise=slope*separation(p,g.boundary[b]);
            upper=std::min(upper,r.boundary_ground[b]+rise);
            lower=std::max(lower,r.boundary_ground[b]-rise);
        }
        const double delta=static_cast<double>(std::max(lower,std::min(static_cast<long double>(q.design_height_m),upper))-z);
        const bool infeasible=lower-upper>q.envelope_tolerance_m;
        r.ground.push_back(z);r.upper.push_back(static_cast<double>(upper));r.lower.push_back(static_cast<double>(lower));
        r.delta.push_back(delta);r.infeasible.push_back(infeasible);
        r.max_gap=std::max(r.max_gap,static_cast<double>(lower-upper));
        if (infeasible) {r.valid=false;if(i<g.nodes.size())++r.invalid_points;else ++r.invalid_boundary;}
        if(i<g.nodes.size()) {
            const double weight=static_cast<double>(g.weights[i]);
            r.cut+=std::max(0.0,-delta)*weight;r.fill+=std::max(0.0,delta)*weight;
            if(delta>q.area_threshold_m)r.fill_area+=weight;
            else if(delta<-q.area_threshold_m)r.cut_area+=weight;
            else r.unchanged_area+=weight;
        }
    }
    for(const auto& cell:g.cells)
        if(r.infeasible[cell[0]]||r.infeasible[cell[1]]||r.infeasible[cell[2]]||r.infeasible[cell[3]])++r.invalid_cells;
    r.net=r.fill-r.cut;r.moved=r.fill+r.cut;return r;
}

const Event& outputEvent(const AuditedBackend& he,std::size_t id) {
    for(const auto& event:he.events)if(event.output==id)return event;
    throw std::runtime_error("Ciphertext event not found");
}
const std::vector<double>& decoded(const AuditedBackend& he,const CipherVector& c) {
    const auto id=he.cipherId(c);
    auto found=std::find(he.decrypted_ids.begin(),he.decrypted_ids.end(),id);
    require(found!=he.decrypted_ids.end(),"A final candidate was not decrypted");
    require(std::count(he.decrypted_ids.begin(),he.decrypted_ids.end(),id)==1,"A final candidate was decrypted twice");
    return he.decrypted_values.at(static_cast<std::size_t>(found-he.decrypted_ids.begin()));
}

struct CaseResult {
    std::string name;std::size_t nodes=0,boundary=0,cells=0,batches=0,encryptions=0,decryptions=0;
    std::size_t slots=0,columns=0,rows=0;
    std::size_t invalid_points=0,invalid_boundary=0,invalid_cells=0;
    bool envelope_valid=false;
    double angle=0,height=0,requested_interval=0,dx=0,dy=0,area=0;
    double max_height_error=0,max_candidate_error=0,max_gap_error=0;
    double plain_cut=0,plain_fill=0,plain_net=0,he_cut=0,he_fill=0,he_net=0,max_volume_error=0,seconds=0;
};

CaseResult checkCase(const std::string& name,const IDtmProvider& provider,const EarthworkSlopeQuery& q,
                     AuditedBackend& he,const fs::path& output,int expected_kind=0) {
    std::cerr<<"BEGIN_CASE "<<name<<" slots="<<he.slotCount()<<'\n';
    const ReferenceGeometry geometry=referenceGeometry(q);
    const Oracle oracle=referenceModel(provider,q,geometry);
    if(expected_kind==1) require(oracle.valid && oracle.cut>1 && oracle.fill<1e-7,"Cut-only fixture is invalid");
    if(expected_kind==2) require(oracle.valid && oracle.fill>1 && oracle.cut<1e-7,"Fill-only fixture is invalid");
    if(expected_kind==3) require(oracle.valid && oracle.cut>1 && oracle.fill>1,"Mixed cut/fill fixture is invalid");
    if(expected_kind==4) require(!oracle.valid && oracle.invalid_boundary>0,"Infeasible-boundary fixture is valid");
    he.reset();he.startCase(name,output);GuardedProvider guarded(provider);TerrainAnalysis analysis(guarded,he);
    const auto start=std::chrono::steady_clock::now();
    const auto encrypted=analysis.earthworkSlopeEncrypted(q);
    require(he.decrypt_calls==0,"Encrypted earthwork API performed an intermediate decrypt");
    const auto mapping=compareGeometry(encrypted.geometry,geometry);
    const auto n=geometry.nodes.size(),boundary=geometry.boundary.size(),targets=n+boundary,slots=he.slotCount();
    const auto batches=(targets+slots-1)/slots;
    require(encrypted.batches.size()==batches,"Incorrect candidate batching");
    require(he.encrypt_calls==batches*(5+4*boundary),"Unexpected encryption count (possible re-encryption or unencrypted H)");
    std::size_t next=0;
    for(const auto& batch:encrypted.batches) {
        const auto count=std::min(slots,targets-next);
        require(batch.begin==next && batch.sample_count==count,"Batch target indices differ or omit targets");
        require(batch.upper_candidates.size()==boundary && batch.lower_candidates.size()==boundary,
                "Boundary candidate array omitted a boundary point");
        const Event& area_product=outputEvent(he,he.cipherId(batch.target_minus_ground));
        require(area_product.operation=="mulPlain","U is not area-weighted before final decryption");
        const Event& difference=outputEvent(he,area_product.input_a);
        require(difference.operation=="sub","U does not contain encrypted H-Z subtraction");
        const Event& target=outputEvent(he,difference.input_a);
        require(target.operation=="encrypt","Design height H was not encrypted at the input boundary");
        const auto& design=he.encrypted_inputs.at(target.output);
        require(design.size()>=count,"Encrypted design-height batch too short");
        for(std::size_t i=0;i<count;++i)require(design[i]==q.design_height_m,"Ciphertext design height differs from caller H");
        next+=count;
    }
    const auto result=analysis.finalizeEarthworkSlope(encrypted);
    require(he.decrypt_calls==batches*(2*boundary+1),"Final candidate decrypt count differs from (2B+1)*batches");
    require(!he.events.empty() && he.events.back().operation=="decrypt","Final backend call is not a terminal decrypt");
    compareGeometry(result.geometry,geometry);
    require(result.envelope_valid==oracle.valid,"Envelope validity differs from independent original-DEM model");
    require(result.infeasible_point_count==oracle.invalid_points &&
            result.infeasible_boundary_point_count==oracle.invalid_boundary &&
            result.infeasible_cell_count==oracle.invalid_cells,"Infeasible node/boundary/cell diagnostics differ");
    require(!result.status.empty(),"Missing envelope status");
    CaseResult summary;summary.name=name;summary.nodes=n;summary.boundary=boundary;summary.cells=geometry.cells.size();
    summary.slots=slots;summary.columns=geometry.columns;summary.rows=geometry.rows;
    summary.height=q.design_height_m;summary.requested_interval=q.sample_interval_m;
    summary.dx=static_cast<double>(geometry.dx);summary.dy=static_cast<double>(geometry.dy);
    summary.batches=batches;summary.encryptions=he.encrypt_calls;summary.decryptions=he.decrypt_calls;
    summary.envelope_valid=oracle.valid;summary.angle=q.slope_angle_degrees;summary.area=static_cast<double>(geometry.area);
    summary.invalid_points=oracle.invalid_points;summary.invalid_boundary=oracle.invalid_boundary;summary.invalid_cells=oracle.invalid_cells;
    summary.max_gap_error=std::abs(result.max_envelope_gap_m-oracle.max_gap);
    require(summary.max_gap_error<1e-6,"Envelope gap diagnostic exceeds 1e-6 m tolerance");
    summary.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream samples(output/(name+"_samples.csv"));require(static_cast<bool>(samples),"Cannot create sample CSV");
    samples<<std::setprecision(std::numeric_limits<double>::max_digits10)
           <<"target,type,lon,lat,weight_m2,plain_ground_m,plain_upper_m,plain_lower_m,plain_delta_m,decoded_delta_m,abs_error_m,infeasible\n";
    const long double tangent=std::tan(static_cast<long double>(q.slope_angle_degrees)*std::acos(-1.0L)/180);
    for(const auto& batch:encrypted.batches) {
        const auto& u=decoded(he,batch.target_minus_ground);
        require(u.size()>=batch.sample_count,"Decoded U is shorter than active batch");
        for(std::size_t i=0;i<batch.sample_count;++i) {
            const auto target=batch.begin+i;
            const auto ref=target<n?mapping[target]:target;
            const XY position=target<n?geometry.nodes[ref]:geometry.boundary[target-n];
            const double weight=target<n?static_cast<double>(geometry.weights[ref])
                                       :static_cast<double>(geometry.dx*geometry.dy/4);
            require(std::isfinite(u[i]),"Non-finite decoded design candidate");
            double upper=std::numeric_limits<double>::infinity(),lower=-upper;
            summary.max_candidate_error=std::max(summary.max_candidate_error,std::abs(u[i]/weight-(q.design_height_m-oracle.ground[ref])));
            for(std::size_t b=0;b<boundary;++b) {
                const auto& f=decoded(he,batch.upper_candidates[b]);const auto& g=decoded(he,batch.lower_candidates[b]);
                require(f.size()>i && g.size()>i,"Decoded boundary candidate is shorter than active batch");
                require(std::isfinite(f[i]) && std::isfinite(g[i]),"Non-finite decoded boundary candidate");
                const double rise=static_cast<double>(tangent*separation(position,geometry.boundary[b]));
                summary.max_candidate_error=std::max(summary.max_candidate_error,
                    std::abs(f[i]/weight-(oracle.boundary_ground[b]+rise-oracle.ground[ref])));
                summary.max_candidate_error=std::max(summary.max_candidate_error,
                    std::abs(g[i]/weight-(oracle.boundary_ground[b]-rise-oracle.ground[ref])));
                upper=std::min(upper,f[i]);lower=std::max(lower,g[i]);
            }
            const double decoded_delta=std::max(lower,std::min(u[i],upper))/weight;
            summary.max_height_error=std::max(summary.max_height_error,std::abs(decoded_delta-oracle.delta[ref]));
            if(oracle.valid && target<n) {
                require(result.point_height_changes_m.size()==n,"Final result omits per-node diagnostic values");
                require(std::abs(result.point_height_changes_m[target]-decoded_delta)<1e-10,
                        "Final point delta differs from final-only decoded candidate min/max");
            }
            const GeoPoint p=geographic(geometry.origin,position);
            samples<<target<<','<<(target<n?"quadrature":"boundary_diagnostic")<<','<<p.lon<<','<<p.lat<<','<<weight<<','
                   <<oracle.ground[ref]<<','<<oracle.upper[ref]<<','<<oracle.lower[ref]<<','<<oracle.delta[ref]<<','
                   <<decoded_delta<<','<<std::abs(decoded_delta-oracle.delta[ref])<<','<<oracle.infeasible[ref]<<'\n';
        }
    }
    require(summary.max_candidate_error<1e-6 && summary.max_height_error<1e-6,
            "Candidate/height error exceeds predeclared 1e-6 m tolerance");
    if(oracle.valid) {
        const double volume_tolerance=std::max(1.0,summary.area)*1e-6;
        const auto volume_check=[&](double actual,double expected) {
            require(std::isfinite(actual),"Valid envelope returned non-finite volume");
            summary.max_volume_error=std::max(summary.max_volume_error,std::abs(actual-expected));
            require(std::abs(actual-expected)<volume_tolerance,"Volume error exceeds area * 1e-6 m tolerance");
        };
        volume_check(result.cut_volume_m3,oracle.cut);volume_check(result.fill_volume_m3,oracle.fill);
        volume_check(result.net_volume_m3,oracle.net);volume_check(result.absolute_moved_volume_m3,oracle.moved);
        // The implementation sums in long double then rounds each returned
        // aggregate separately. Re-subtracting rounded double fields can differ
        // by one ulp at billion-cubic-metre magnitudes; bound only that rounding.
        const double identity_roundoff = 2 * std::numeric_limits<double>::epsilon() *
            (std::abs(result.fill_volume_m3) + std::abs(result.cut_volume_m3) + 1.0);
        require(std::abs(result.net_volume_m3-(result.fill_volume_m3-result.cut_volume_m3))<=identity_roundoff,
                "Net sign differs from fill minus cut");
        require(std::abs(result.cut_area_m2-oracle.cut_area)<1e-6 &&
                std::abs(result.fill_area_m2-oracle.fill_area)<1e-6 &&
                std::abs(result.unchanged_area_m2-oracle.unchanged_area)<1e-6,"Area classification differs from independent threshold rule");
        summary.plain_cut=oracle.cut;summary.plain_fill=oracle.fill;summary.plain_net=oracle.net;
        summary.he_cut=result.cut_volume_m3;summary.he_fill=result.fill_volume_m3;summary.he_net=result.net_volume_m3;
    } else {
        require(std::isnan(result.cut_volume_m3) && std::isnan(result.fill_volume_m3) &&
                std::isnan(result.net_volume_m3) && std::isnan(result.absolute_moved_volume_m3) &&
                std::isnan(result.cut_area_m2) && std::isnan(result.fill_area_m2) && std::isnan(result.unchanged_area_m2),
                "Infeasible envelope was reported as a valid earthwork volume/area");
    }
    std::ofstream trace(output/(name+"_trace.csv"));require(static_cast<bool>(trace),"Cannot create trace CSV");
    trace<<"event,operation,output_cipher,input_a,input_b\n";
    for(std::size_t i=0;i<he.events.size();++i) {
        const auto& event=he.events[i];trace<<i<<','<<event.operation<<','<<event.output<<','<<event.input_a<<','<<event.input_b<<'\n';
    }
    return summary;
}

RasterDtmProvider planeProvider(double east_slope=0,double north_slope=0) {
    RasterGrid grid;grid.origin={0,0};grid.width=96;grid.height=96;grid.cell_size_m=5;
    for(std::size_t y=0;y<grid.height;++y)for(std::size_t x=0;x<grid.width;++x)
        grid.elevations_m.push_back(100+east_slope*5*x+north_slope*5*y);
    DtmMetadata metadata;metadata.name="Independent analytic plane fixture";metadata.source="validation fixture";
    metadata.crs="local equirectangular over geographic coordinates";metadata.resolution_m=5;
    return RasterDtmProvider(metadata,std::move(grid));
}
EarthworkSlopeQuery rectangle(const RasterDtmProvider& provider,double x,double y,double width,double height,
                              double H,double angle,double step=11,double boundary_step=100) {
    EarthworkSlopeQuery q;
    q.polygon={provider.pointFromMeters(x,y),provider.pointFromMeters(x+width,y),
               provider.pointFromMeters(x+width,y+height),provider.pointFromMeters(x,y+height)};
    q.design_height_m=H;q.slope_angle_degrees=angle;q.sample_interval_m=step;q.boundary_interval_m=boundary_step;
    return q;
}
void writeCase(std::ostream& out,const CaseResult& r) {
    out<<"{\"name\":"<<jsonString(r.name)<<",\"passed\":true,\"envelope_valid\":"<<(r.envelope_valid?"true":"false")
       <<",\"nodes\":"<<r.nodes<<",\"boundary_points\":"<<r.boundary<<",\"cells\":"<<r.cells<<",\"batches\":"<<r.batches
       <<",\"slots\":"<<r.slots<<",\"grid_columns\":"<<r.columns<<",\"grid_rows\":"<<r.rows
       <<",\"design_height_m\":"<<r.height<<",\"requested_grid_interval_m\":"<<r.requested_interval
       <<",\"grid_dx_m\":"<<r.dx<<",\"grid_dy_m\":"<<r.dy
       <<",\"encrypt_calls\":"<<r.encryptions<<",\"decrypt_calls\":"<<r.decryptions<<",\"angle_degrees\":"<<r.angle
       <<",\"integration_area_m2\":"<<r.area<<",\"max_height_error_m\":"<<r.max_height_error
       <<",\"max_candidate_error_m\":"<<r.max_candidate_error<<",\"max_envelope_gap_error_m\":"<<r.max_gap_error
       <<",\"max_volume_error_m3\":"<<r.max_volume_error<<",\"infeasible_nodes\":"<<r.invalid_points
       <<",\"infeasible_boundary_points\":"<<r.invalid_boundary<<",\"infeasible_cells\":"<<r.invalid_cells
       <<",\"elapsed_seconds\":"<<r.seconds;
    if(r.envelope_valid) {
        out<<",\"plain_cut_m3\":"<<r.plain_cut<<",\"plain_fill_m3\":"<<r.plain_fill<<",\"plain_net_m3\":"<<r.plain_net
           <<",\"he_cut_m3\":"<<r.he_cut<<",\"he_fill_m3\":"<<r.he_fill<<",\"he_net_m3\":"<<r.he_net
           <<",\"cut_abs_error_m3\":"<<std::abs(r.he_cut-r.plain_cut)
           <<",\"fill_abs_error_m3\":"<<std::abs(r.he_fill-r.plain_fill)
           <<",\"net_abs_error_m3\":"<<std::abs(r.he_net-r.plain_net);
        const auto relative=[&](const char* field,double actual,double expected) {
            out<<",\""<<field<<"\":";
            if(expected==0) out<<"null";
            else out<<100*std::abs(actual-expected)/std::abs(expected);
        };
        relative("cut_relative_error_percent",r.he_cut,r.plain_cut);
        relative("fill_relative_error_percent",r.he_fill,r.plain_fill);
        relative("net_relative_error_percent",r.he_net,r.plain_net);
    } else {
        out<<",\"plain_cut_m3\":null,\"plain_fill_m3\":null,\"plain_net_m3\":null,\"he_cut_m3\":null,\"he_fill_m3\":null,\"he_net_m3\":null"
           <<",\"cut_abs_error_m3\":null,\"fill_abs_error_m3\":null,\"net_abs_error_m3\":null"
           <<",\"cut_relative_error_percent\":null,\"fill_relative_error_percent\":null,\"net_relative_error_percent\":null";
    }
    out<<'}';
}

std::size_t invalidInputs(const RasterDtmProvider& provider,AuditedBackend& he) {
    const auto base=rectangle(provider,80,80,20,20,105,45);
    std::size_t passed=0;
    const auto reject=[&](const EarthworkSlopeQuery& query) {
        he.reset();GuardedProvider guarded(provider);TerrainAnalysis analysis(guarded,he);bool rejected=false;
        try{(void)analysis.earthworkSlopeEncrypted(query);}catch(const std::exception&){rejected=true;}
        require(rejected,"Invalid query was not rejected");
        require(he.encrypt_calls==0 && he.decrypt_calls==0,"Invalid query reached encryption or decryption");++passed;
    };
    auto q=base;q.design_height_m=std::numeric_limits<double>::quiet_NaN();reject(q);
    q=base;q.slope_angle_degrees=std::numeric_limits<double>::infinity();reject(q);
    q=base;q.slope_angle_degrees=90;reject(q);
    q=base;q.slope_angle_degrees=-1;reject(q);
    q=base;q.sample_interval_m=0;reject(q);
    q=base;q.boundary_interval_m=-1;reject(q);
    q=base;q.polygon[0].lon=std::numeric_limits<double>::quiet_NaN();reject(q);
    q=base;q.polygon.resize(2);reject(q);
    q=base;q.max_grid_cells=1;reject(q);
    q=base;q.max_boundary_points=3;reject(q);
    q=base;q.max_candidate_ciphertexts=1;reject(q);
    return passed;
}
} // namespace

int main(int argc,char** argv) {
    if(argc<3 || argc>4) {std::cerr<<"usage: earthwork_slope_validation <taiwan.tdb> <output_dir> [runs=1]\n";return 2;}
    try {
        const auto taiwan=loadValidationTdbDataset(argv[1]);
        const fs::path output=argv[2];fs::create_directories(output);
        const int runs=argc==4?std::stoi(argv[3]):1;require(runs>0 && runs<=100,"runs must be in [1,100]");
        auto flat=planeProvider(),mixed=planeProvider(0.12,0.04),steep=planeProvider(2,0);
        std::vector<std::vector<CaseResult>> all;std::size_t invalid_count=0;
        std::cout<<std::setprecision(std::numeric_limits<double>::max_digits10);
        for(int run=1;run<=runs;++run) {
            OpenFheBackend backend(16);AuditedBackend he(backend);
            const fs::path run_dir=output/("run_"+std::to_string(run));fs::create_directories(run_dir);
            std::vector<CaseResult> results;
            for(const double angle:{30.0,45.0,60.0}) {
                const std::string suffix=std::to_string(static_cast<int>(angle));
                // A 12 m target offset activates the cone bound at 30/45 deg,
                // while 60 deg reaches H at interior nodes: angle must matter.
                results.push_back(checkCase("flat_cut_"+suffix,flat,rectangle(flat,80,80,20,20,88,angle),he,run_dir,1));
                results.push_back(checkCase("flat_fill_"+suffix,flat,rectangle(flat,80,80,20,20,112,angle),he,run_dir,2));
                results.push_back(checkCase("plane_both_"+suffix,mixed,rectangle(mixed,80,80,20,20,114.4,angle),he,run_dir,3));
                const double width=std::min(30.0,taiwan.widthMeters()*0.15),height=std::min(24.0,taiwan.heightMeters()*0.15);
                const double x=taiwan.widthMeters()*0.35,y=taiwan.heightMeters()*0.35;
                const double target=taiwan.elevationAt(taiwan.pointFromMeters(x+width/2,y+height/2));
                results.push_back(checkCase("taiwan_"+suffix,taiwan,rectangle(taiwan,x,y,width,height,target,angle,
                                                                          std::max(width,height)/2.2,std::max(width,height)*2),he,run_dir));
            }
            results.push_back(checkCase("multiple_batches_shared_nodes",mixed,rectangle(mixed,80,80,40,40,116,45,11,12),he,run_dir,3));
            auto concave=rectangle(flat,80,80,30,30,104,45,8,100);
            concave.polygon={flat.pointFromMeters(80,80),flat.pointFromMeters(110,80),flat.pointFromMeters(110,90),
                             flat.pointFromMeters(90,90),flat.pointFromMeters(90,110),flat.pointFromMeters(80,110)};
            results.push_back(checkCase("concave_selected_cells",flat,concave,he,run_dir,2));
            results.push_back(checkCase("infeasible_boundary",steep,rectangle(steep,80,80,20,20,280,30),he,run_dir,4));
            {
                // Match final_analysis_demo.cpp's paper-region construction
                // exactly, using a separate 2048-slot context to fit one batch.
                const double span=std::max(2.0,std::min(taiwan.widthMeters(),taiwan.heightMeters())*0.55);
                const double x0=span*0.12,y0=span*0.12;
                const auto point=[&](double x,double y){return taiwan.pointFromMeters(x,y);};
                EarthworkSlopeQuery paper;
                paper.polygon={point(x0+span*0.11,y0+span*0.11),point(x0+span*0.78,y0+span*0.11),
                               point(x0+span*0.78,y0+span*0.78),point(x0+span*0.11,y0+span*0.78)};
                paper.design_height_m=14.76;paper.slope_angle_degrees=45;
                paper.sample_interval_m=std::max(2.0,span/18.0);
                OpenFheBackend paper_backend(2048);AuditedBackend paper_he(paper_backend);
                auto checked=checkCase("paper_taiwan_45",taiwan,paper,paper_he,run_dir);
                require(checked.batches==1,"Paper-region validation does not fit its required one batch");
                results.push_back(std::move(checked));
            }
            invalid_count=invalidInputs(flat,he);
            // Exercise the convenience wrapper independently of the split API.
            he.reset();GuardedProvider guard(flat);TerrainAnalysis wrapped(guard,he);
            const auto q=rectangle(flat,80,80,20,20,104,45);
            const auto smoke=wrapped.earthworkSlope(q);
            const auto ref=referenceGeometry(q);const auto oracle=referenceModel(flat,q,ref);
            require(smoke.envelope_valid && std::abs(smoke.fill_volume_m3-oracle.fill)<ref.area*1e-6L,
                    "Convenience wrapper differs from independent finalization baseline");
            const auto nb=(ref.nodes.size()+ref.boundary.size()+he.slotCount()-1)/he.slotCount();
            require(he.decrypt_calls==nb*(2*ref.boundary.size()+1),"Convenience wrapper performed extra decryption");
            std::ofstream json(run_dir/"results.json");require(static_cast<bool>(json),"Cannot create run JSON");
            json<<std::setprecision(std::numeric_limits<double>::max_digits10)<<"{\"run\":"<<run<<",\"fresh_keys\":true,\"slots_by_case\":true,\"invalid_inputs_passed\":"<<invalid_count<<",\"cases\":[";
            for(std::size_t i=0;i<results.size();++i) {
                if(i) json<<',';
                writeCase(json,results[i]);
            }
            json<<"]}\n";
            std::cout<<"{\"run\":"<<run<<",\"cases_passed\":"<<results.size()<<",\"invalid_inputs_passed\":"<<invalid_count<<"}\n";
            all.push_back(std::move(results));
        }
        std::ofstream summary(output/"summary.json");require(static_cast<bool>(summary),"Cannot create summary JSON");
        summary<<std::setprecision(std::numeric_limits<double>::max_digits10)
               <<"{\"passed\":true,\"runs\":"<<runs<<",\"numeric_tolerance_m\":1e-6,\"volume_tolerance\":\"integration_area_m2 * 1e-6 m\",\"intermediate_decrypts\":0,\"he_operations_after_first_decrypt\":0,\"direct_plain_elevation_in_production\":false,\"design_height_encrypted\":true,\"net_convention\":\"fill - cut\",\"baseline\":\"Independent winding-number cell selection, shared FourCorners weights, independent equal boundary subdivisions and original DEM elevationAt, followed by max(g,min(H,f))-Z.\",\"scope\":\"Same discrete boundary-envelope model only; no claim of equality to PPT GIS buffer/difference/TIN implementation or field survey truth.\",\"cases\":[";
        bool first=true;
        for(std::size_t r=0;r<all.size();++r)for(const auto& c:all[r]) {
            if(!first) summary<<',';
            first=false;
            summary<<"{\"run\":"<<(r+1)<<",\"result\":";
            writeCase(summary,c);
            summary<<'}';
        }
        summary<<"]}\n";return 0;
    } catch(const std::exception& e) {std::cerr<<"earthwork slope validation failed: "<<e.what()<<'\n';return 1;}
}


#include "openfhe_dtm/EarthworkTheta.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/TdbDataset.hpp"

#include <algorithm>
#include <cmath>
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

namespace {
constexpr long double radius = 6371008.8L;
const long double pi = std::acos(-1.0L);
using XY = std::pair<long double, long double>;

void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}
void near(long double actual, long double expected, long double tolerance,
          const std::string& message) {
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, message);
}
GeoPoint globalPoint(long double east, long double north) {
    return {static_cast<double>(east / radius * 180 / pi),
            static_cast<double>(north / radius * 180 / pi)};
}
XY offsets(const GeoPoint& origin, const GeoPoint& point) {
    return {(static_cast<long double>(point.lon) - origin.lon) * pi / 180 * radius *
                std::cos(static_cast<long double>(origin.lat) * pi / 180),
            (static_cast<long double>(point.lat) - origin.lat) * pi / 180 * radius};
}
GeoPoint displaced(const GeoPoint& origin, long double east, long double north) {
    return {static_cast<double>(origin.lon + east / radius /
                 std::cos(static_cast<long double>(origin.lat) * pi / 180) * 180 / pi),
            static_cast<double>(origin.lat + north / radius * 180 / pi)};
}
long double groundAt(long double east, long double north) {
    // Exactly bilinear, so an independent analytic value is available without
    // calling the production interpolation or its plaintext elevation API.
    return 100 + 0.07L * east - 0.03L * north + 0.001L * east * north;
}

class AnalyticDtm final : public IDtmProvider {
public:
    explicit AnalyticDtm(const IDtmProvider* source = nullptr) : source_(source) {}
    const DtmMetadata& metadata() const override { return metadata_; }
    InterpolationStencil interpolationStencil(const GeoPoint& p) const override {
        points.push_back(p);
        if (source_) return source_->interpolationStencil(p);
        const auto [x, y] = offsets({0, 0}, p);
        const auto x0 = std::floor(x), y0 = std::floor(y);
        return {static_cast<double>(groundAt(x0, y0)),
                static_cast<double>(groundAt(x0 + 1, y0)),
                static_cast<double>(groundAt(x0, y0 + 1)),
                static_cast<double>(groundAt(x0 + 1, y0 + 1)),
                static_cast<double>(x - x0), static_cast<double>(y - y0)};
    }
    double elevationAt(const GeoPoint&) const override {
        throw std::runtime_error("Production path called plaintext elevationAt");
    }
    DtmTile loadTile(int, int) const override {
        throw std::runtime_error("Unexpected tile lookup in the plane test");
    }
    mutable Route points;
private:
    DtmMetadata metadata_;
    const IDtmProvider* source_;
};

struct Node {
    std::string operation;
    const void* input = nullptr;
    std::vector<double> plain;
};

// Real CKKS operations, with a graph and a strict terminal-decryption barrier.
class AuditedBackend final : public IHeBackend {
public:
    explicit AuditedBackend(const IHeBackend& backend) : backend_(backend) {}
    std::size_t slotCount() const override { return backend_.slotCount(); }
    PlainVector encode(const std::vector<double>& v, int l = 0) const override {
        before();
        return backend_.encode(v, l);
    }
    CipherVector encrypt(const PlainVector& p) const override {
        before(); ++encryptions;
        return record(backend_.encrypt(p), "encrypt", nullptr, p.values);
    }
    PlainVector decrypt(const CipherVector& c) const override {
        require(nodes.at(c.native.get()).operation == "sumSlots",
                "Decryption must receive the final encrypted slot sum");
        require(decryptions == 0, "More than one production decryption");
        ++decryptions;
        return backend_.decrypt(c);
    }
    CipherVector add(const CipherVector& a, const CipherVector& b) const override {
        before();
        return record(backend_.add(a, b), "add", a.native.get());
    }
    CipherVector mulPlain(const CipherVector& a, const PlainVector& p) const override {
        before();
        return record(backend_.mulPlain(a, p), "mulPlain", a.native.get(), p.values);
    }
    CipherVector subPlain(const CipherVector& a, const PlainVector& p) const override {
        before();
        return record(backend_.subPlain(a, p), "subPlain", a.native.get(), p.values);
    }
    CipherVector plainSub(const PlainVector& p, const CipherVector& a) const override {
        before();
        return record(backend_.plainSub(p, a), "plainSub", a.native.get(), p.values);
    }
    CipherVector sumSlots(const CipherVector& a) const override {
        before(); ++sums;
        return record(backend_.sumSlots(a), "sumSlots", a.native.get());
    }
    bool supportsBootstrap() const override { return backend_.supportsBootstrap(); }
    PlainVector encodeComplex(const std::vector<std::complex<double>>&, int = 0) const override {
        throw std::runtime_error("Unexpected complex encoding");
    }
    std::vector<std::complex<double>> decryptComplex(const CipherVector&) const override {
        throw std::runtime_error("Unexpected complex decryption");
    }
    CipherVector sub(const CipherVector&, const CipherVector&) const override {
        throw std::runtime_error("Unexpected ciphertext subtraction; expected plainSub");
    }
    CipherVector mul(const CipherVector&, const CipherVector&) const override {
        throw std::runtime_error("Unexpected ciphertext multiplication");
    }
    CipherVector mulScalar(const CipherVector&, std::complex<double>) const override {
        throw std::runtime_error("Unexpected scalar multiplication");
    }
    CipherVector conjugate(const CipherVector&) const override {
        throw std::runtime_error("Unexpected conjugation");
    }
    CipherVector rotate(const CipherVector&, int) const override {
        throw std::runtime_error("Unexpected explicit rotation; expected sumSlots");
    }
    CipherVector bootstrap(const CipherVector&) const override {
        throw std::runtime_error("Unexpected bootstrap");
    }
    CipherVector evaluateFunction(const CipherVector&, double, double, std::size_t,
                                 double (*)(double)) const override {
        throw std::runtime_error("Unexpected function evaluation");
    }
    CipherVector evaluateChebyshevSeries(const CipherVector&, const std::vector<double>&,
                                         double, double) const override {
        throw std::runtime_error("Unexpected polynomial evaluation");
    }
    void saveCipher(const CipherVector&, const std::filesystem::path&) const override {
        throw std::runtime_error("Unexpected ciphertext persistence");
    }
    CipherVector loadCipher(const std::filesystem::path&) const override {
        throw std::runtime_error("Unexpected ciphertext load");
    }
    mutable std::size_t encryptions = 0, decryptions = 0, sums = 0;
    mutable std::map<const void*, Node> nodes;
private:
    void before() const {
        require(decryptions == 0, "HE or encoding operation after final decryption");
    }
    CipherVector record(CipherVector c, const char* operation, const void* input,
                        std::vector<double> plain = {}) const {
        require(c.native != nullptr, "Null ciphertext");
        // Preserve lifetimes, so allocator pointer reuse cannot corrupt the graph.
        handles_.push_back(c.native);
        nodes[c.native.get()] = {operation, input, std::move(plain)};
        return c;
    }
    const IHeBackend& backend_;
    mutable std::vector<std::shared_ptr<void>> handles_;
};

bool windingInside(const std::vector<XY>& polygon, long double x, long double y) {
    int winding = 0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto [ax, ay] = polygon[i];
        const auto [bx, by] = polygon[(i + 1) % polygon.size()];
        const long double cross = (bx - ax) * (y - ay) - (x - ax) * (by - ay);
        if (ay <= y && by > y && cross > 0) ++winding;
        if (ay > y && by <= y && cross < 0) --winding;
    }
    return winding != 0;
}
struct Reference {
    Route points;
    std::vector<long double> design;
    long double net = 0, east = 0, north = 0;
    bool negative_projection = false, positive_projection = false;
};
Reference reference(const EarthworkThetaQuery& q, const IDtmProvider* source = nullptr) {
    Reference result;
    std::vector<XY> polygon;
    for (const auto& point : q.polygon) polygon.push_back(offsets(q.polygon.front(), point));
    const auto edge = polygon.at(1);
    const long double length = std::hypot(edge.first, edge.second);
    result.east = edge.first / length; result.north = edge.second / length;
    long double xmin = polygon[0].first, xmax = xmin;
    long double ymin = polygon[0].second, ymax = ymin;
    for (const auto& p : polygon) {
        xmin = std::min(xmin, p.first); xmax = std::max(xmax, p.first);
        ymin = std::min(ymin, p.second); ymax = std::max(ymax, p.second);
    }
    const long double step = q.sample_interval_m;
    const long double tangent = std::tan(static_cast<long double>(q.slope_angle_degrees) * pi / 180);
    for (long double y = ymin + step / 2; y <= ymax; y += step) {
        for (long double x = xmin + step / 2; x <= xmax; x += step) {
            if (!windingInside(polygon, x, y)) continue;
            const long double d = x * result.east + y * result.north;
            result.negative_projection |= d < 0; result.positive_projection |= d > 0;
            const long double design = q.design_height_m + d * tangent;
            const auto point = displaced(q.polygon.front(), x, y);
            const auto global = offsets({0, 0}, point);
            result.points.push_back(point); result.design.push_back(design);
            long double ground = groundAt(global.first, global.second);
            if (source) {
                // Validation only: independent long-double arithmetic directly
                // on original DEM corners, never decrypted production heights.
                const auto s = source->interpolationStencil(point);
                const long double dx = s.dx, dy = s.dy;
                ground = s.z00_m * (1-dx) * (1-dy) + s.z10_m * dx * (1-dy) +
                         s.z01_m * (1-dx) * dy + s.z11_m * dx * dy;
            }
            result.net += (design - ground) * step * step;
        }
    }
    return result;
}

EarthworkThetaQuery makeQuery(std::initializer_list<XY> points, double theta = 0,
                             double height = 103, double interval = 5) {
    EarthworkThetaQuery q;
    for (const auto& p : points) q.polygon.push_back(globalPoint(p.first, p.second));
    q.design_height_m = height; q.slope_angle_degrees = theta; q.sample_interval_m = interval;
    return q;
}

void check(const char* name, const EarthworkThetaQuery& q, const IHeBackend& real,
           std::size_t expected_count, bool signed_projection = false, bool one_shot = false,
           const IDtmProvider* source = nullptr) {
    AnalyticDtm dtm(source);
    AuditedBackend he(real);
    const auto expected = reference(q, source);
    require(expected.points.size() == expected_count, std::string(name) + ": oracle point count");
    if (signed_projection)
        require(expected.negative_projection && expected.positive_projection,
                "Signed-projection test must exercise both sides of H");
    EarthworkThetaResult actual;
    if (one_shot) {
        actual = earthworkTheta(q, dtm, he);
    } else {
        const auto encrypted = earthworkThetaEncrypted(q, dtm, he);
        require(he.decryptions == 0, "Intermediate decryption");
        require(std::isnan(encrypted.geometry.net_volume_m3), "Encrypted result exposed a scalar volume");
        const auto& sum = he.nodes.at(encrypted.net_volume.native.get());
        require(sum.operation == "sumSlots", "Missing ciphertext slot sum");
        const auto& area = he.nodes.at(sum.input);
        require(area.operation == "mulPlain", "Missing area multiplication before sum");
        const auto& difference = he.nodes.at(area.input);
        require(difference.operation == "plainSub", "Expected public-design-minus-encrypted-ground subtraction");
        require(area.plain.size() == real.slotCount(), "Area mask length differs from slots");
        require(difference.plain.size() == real.slotCount(), "Encoded design length differs from slots");
        const double cell_area = q.sample_interval_m * q.sample_interval_m;
        for (std::size_t i = 0; i < real.slotCount(); ++i) {
            const bool valid = i < expected.points.size();
            require(area.plain[i] == (valid ? cell_area : 0.0),
                    "Area sign or padding changed; expected valid +Ac, padding zero");
            near(difference.plain[i], valid ? expected.design[i] : 0, 1e-9L,
                 "Design must be H + signed projection * tan(theta), with zero padding");
        }
        actual = finalizeEarthworkTheta(encrypted, he);
    }
    require(he.encryptions == 4 && he.sums == 1 && he.decryptions == 1,
            "Expected four corner encryptions, one EvalSum and one terminal decryption");
    require(actual.sample_count == expected_count && dtm.points.size() == expected_count,
            "Sampling or stencil count changed");
    near(actual.sample_area_m2, q.sample_interval_m * q.sample_interval_m, 1e-12L, "Cell area changed");
    near(actual.total_area_m2, expected_count * actual.sample_area_m2, 1e-9L, "Total area changed");
    near(actual.direction_east, expected.east, 1e-12L, "East direction is wrong");
    near(actual.direction_north, expected.north, 1e-12L, "North direction is wrong");
    near(actual.origin.lon, q.polygon.front().lon, 1e-15L, "H origin longitude changed");
    near(actual.origin.lat, q.polygon.front().lat, 1e-15L, "H origin latitude changed");
    for (std::size_t i = 0; i < expected_count; ++i) {
        const auto difference = offsets(expected.points[i], dtm.points[i]);
        near(difference.first, 0, 1e-6L, "Sample east coordinate mismatch");
        near(difference.second, 0, 1e-6L, "Sample north coordinate mismatch");
    }
    const long double tolerance = std::max(1e-6L, actual.total_area_m2 * 1e-8L);
    near(actual.net_volume_m3, expected.net, tolerance, std::string(name) + ": encrypted fill-minus-cut mismatch");
    std::cout << "PASS " << name << " n=" << actual.sample_count
              << " net_fill_minus_cut_m3=" << actual.net_volume_m3
              << " reference_m3=" << static_cast<double>(expected.net)
              << " error_m3=" << static_cast<double>(std::abs(actual.net_volume_m3 - expected.net)) << '\n';
}

void reject(const char* name, const EarthworkThetaQuery& q, const IHeBackend& real) {
    AnalyticDtm dtm;
    AuditedBackend he(real);
    bool rejected = false;
    try { (void)earthworkThetaEncrypted(q, dtm, he); }
    catch (const std::invalid_argument&) { rejected = true; }
    catch (const std::length_error&) { rejected = true; }
    catch (const std::runtime_error&) { rejected = true; }
    require(rejected, std::string("Invalid input accepted: ") + name);
    require(he.decryptions == 0 && he.encryptions == 0, "Invalid input reached cryptographic computation");
    require(dtm.points.empty(), "Invalid input reached terrain interpolation");
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc <= 2, "Usage: theta_plane_validation [taiwan.tdb]");
        std::cout << std::setprecision(17);
        OpenFheBackend small(16, {}, false);
        const auto rectangle = makeQuery({{0,0},{20,0},{20,20},{0,20}});
        check("theta_zero_full_slots", rectangle, small, 16);
        auto q = rectangle; q.slope_angle_degrees = 15;
        check("east_positive", q, small, 16);
        q.slope_angle_degrees = -15; check("east_negative", q, small, 16);
        check("north", makeQuery({{0,0},{0,20},{20,20},{20,0}}, 20), small, 16);
        check("reversed_edge", makeQuery({{20,0},{0,0},{0,20},{20,20}}, 20), small, 16);
        check("diagonal", makeQuery({{0,0},{12,9},{3,21},{-9,12}}, 20), small, 10);
        check("signed_projection", makeQuery({{0,0},{10,0},{10,20},{-20,20},{-20,-10},{0,-10}},
                                               12, 103, 10), small, 8, true);
        check("one_point", makeQuery({{0,0},{5,0},{5,5},{0,5}}, 0, 100), small, 1);
        check("partial_padding", makeQuery({{0,0},{15,0},{15,10},{0,10}}, -10), small, 6);
        check("one_padding_slot", makeQuery({{0,0},{15,0},{15,25},{0,25}}, 10), small, 15);
        check("one_shot_wrapper", makeQuery({{0,0},{5,0},{5,5},{0,5}}, 10), small, 1, false, true);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        for (double height : {nan, inf, -inf}) { q = rectangle; q.design_height_m = height; reject("height", q, small); }
        for (double angle : {nan, inf, -inf, -90.0, 90.0, 91.0}) {
            q = rectangle; q.slope_angle_degrees = angle; reject("theta", q, small);
        }
        for (double interval : {nan, inf, 0.0, -1.0}) {
            q = rectangle; q.sample_interval_m = interval; reject("interval", q, small);
        }
        q = rectangle; q.polygon.resize(2); reject("too_few_vertices", q, small);
        q = rectangle; q.polygon[1] = q.polygon[0]; reject("coincident_first_edge", q, small);
        reject("slot_overflow", makeQuery({{0,0},{25,0},{25,20},{0,20}}), small);
        q = rectangle; q.sample_interval_m = 100; reject("no_samples", q, small);

        OpenFheBackend large(2048, {}, false);
        check("144_points_2048_slots", makeQuery({{0,0},{60,0},{60,60},{0,60}}, 0, 14.76), large, 144);
        check("144_points_theta", makeQuery({{0,0},{60,0},{60,60},{0,60}}, 25, 102), large, 144);
        if (argc == 2) {
            const auto terrain = loadValidationTdbDataset(argv[1]);
            const double span = std::max(2.0, std::min(terrain.widthMeters(), terrain.heightMeters()) * 0.55);
            const double start = span * 0.12;
            auto point = [&](double x, double y) { return terrain.pointFromMeters(x, y); };
            EarthworkThetaQuery taiwan;
            taiwan.polygon = {point(start+span*0.11, start+span*0.11),
                              point(start+span*0.78, start+span*0.11),
                              point(start+span*0.78, start+span*0.78),
                              point(start+span*0.11, start+span*0.78)};
            taiwan.design_height_m = 14.76;
            taiwan.sample_interval_m = std::max(2.0, span/18.0);
            check("taiwan_original_grid_theta_zero", taiwan, large, 144, false, false, &terrain);
            taiwan.slope_angle_degrees = 10;
            check("taiwan_original_grid_theta_ten", taiwan, large, 144, false, false, &terrain);
        }
        std::cout << "PASS all " << (argc == 2 ? 15 : 13)
                  << " numerical/trace cases and 17 invalid inputs; net = fill - cut\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "theta plane validation failed: " << e.what() << '\n';
        return 1;
    }
}

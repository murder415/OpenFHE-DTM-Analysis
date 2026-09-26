#pragma once
// Rotation-based bilinear interpolation over the encrypted packed DEM.
//
// The DEM provider encrypts whole 128x128 stitches (2x2 tiles of 64, row-major,
// one DEM cell per slot, z / kDemScale) two per 32,768-slot ciphertext, the same
// layout as EncryptedDemStore. The analysis side never sees corner heights: for a
// query point in cell (x0, y0) it only knows the public slot of that cell and the
// fractions (dx, dy). Rotating the stitch ciphertext by 1, W and W+1 (W = 128)
// brings the right, upper and diagonal neighbours onto that slot; four plaintext
// weight masks then give the bilinear height in the point's slot.
//
// Stitches overlap by one tile (stitch (tx, ty) covers tiles tx..tx+1, ty..ty+1)
// and a point uses the stitch of its own tile, so the local cell is at most 63
// and all four corners always lie in that one stitch: no cross-stitch case.
#include "openfhe_dtm/EncryptedDem.hpp"

#include <map>
#include <utility>
#include <vector>

namespace openfhe_dtm {

struct SightSurfaceQuery {
    GeoPoint base;
    GeoPoint direction;
    double angle_degrees = 15.0;   // regulation-plane angle theta
    double max_distance_m = 90.0;
    double sample_interval_m = 10.0;
};

struct IndexedBatch {
    CipherVector ciphertext;          // value of point points[k] sits in slot slots[k]
    std::vector<std::size_t> slots;
    std::vector<std::size_t> points;  // indices into the query route
};

struct IndexedCipher {
    std::size_t point_count = 0;
    std::vector<IndexedBatch> batches;
};

struct PackedElevationStats {
    std::size_t packs = 0;            // encrypted packed ciphertexts (provider side)
    double provider_encrypt_s = 0.0;  // time spent building + encrypting packs
    std::size_t rotations = 0;        // server-side rotations
    std::size_t batches = 0;          // interpolation batches before merging
};

class PackedElevation {
public:
    PackedElevation(const RasterDtmProvider& dtm, const IHeBackend& he);

    // Provider side: encrypt every packed stitch the route needs (cached).
    void encryptFor(const Route& route);
    // Server side: interpolation using only rotations and plaintext weights.
    IndexedCipher interpolate(const Route& route) const;
    // Adds batches whose slot sets are disjoint so fewer ciphertexts are decrypted.
    IndexedCipher mergeDisjoint(const IndexedCipher& in) const;
    // Owner side: decrypt and return values in route order.
    std::vector<double> decrypt(const IndexedCipher& c, std::size_t* decryptions = nullptr) const;

    const PackedElevationStats& stats() const { return stats_; }
    const IHeBackend& backend() const { return he_; }

    struct Cell {
        int tx = 0, ty = 0;
        std::size_t slot = 0;
        double dx = 0.0, dy = 0.0;
    };
    Cell locate(const GeoPoint& p) const;

private:
    using PackKey = std::pair<int, int>;  // (packLeftX0(tx), ty)
    static PackKey keyOf(const Cell& c) { return {packLeftX0(c.tx), c.ty}; }
    std::vector<double> stitchValues(int tx, int ty) const;

    const RasterDtmProvider& dtm_;
    const IHeBackend& he_;
    std::map<PackKey, CipherVector> packs_;
    mutable PackedElevationStats stats_;
};

struct PackedSightResult {
    std::vector<double> distances_m;
    std::vector<double> buildable_height_m;
    std::size_t decryptions = 0;
};

// Sight surface on the packed path: G(d) by rotation interpolation, z_B isolated
// by a mask and broadcast by a slot sum, r = z_B + d tan(theta), h = r - G, and
// only h is decrypted.
PackedSightResult sightSurfacePacked(const SightSurfaceQuery& query, PackedElevation& elevation);

} // namespace openfhe_dtm

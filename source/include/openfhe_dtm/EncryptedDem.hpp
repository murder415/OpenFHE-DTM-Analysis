#pragma once
#include "openfhe_dtm/DtmProvider.hpp"
#include "openfhe_dtm/HeBackend.hpp"
#include <filesystem>
#include <map>
#include <optional>
#include <tuple>
#include <vector>
namespace openfhe_dtm {
inline constexpr std::size_t kDemTileWidth=64, kDemStitchWidth=128;
inline constexpr std::size_t kDemStitchSlots=16384, kDemPackSlots=32768;
inline constexpr double kDemScale=1024.0;
inline int packSlotOffset(int x){return (x&2)?16384:0;} inline int packLeftX0(int x){return x&~2;}
inline int packStitchBit(int x){return (x&2)?1:0;} inline int packImagChannel(int y){return (y&2)?1:0;}
inline int pack16BaseX0(int x){return x&~2;} inline int pack16BaseY0(int y){return y&~2;}
inline int pack16StitchBit(int x,int y){return ((x&2)?1:0)|((y&2)?2:0);}
struct StitchAddress{int level=0,x0=0,y0=0;};
struct EncryptedDemStitch{int level=0,x0=0,y0=0;std::size_t slot_offset=0,valid_width=128,valid_height=128;CipherVector ciphertext;};
struct LevelEncryptStats{std::size_t total=0,encrypted=0,skipped=0;};
class EncryptedDemStore{
public:
 EncryptedDemStore(const IDtmProvider&,const IHeBackend&,std::filesystem::path,bool compress=false);
 EncryptedDemStitch loadOrEncrypt(int,int,int); EncryptedDemStitch load(int,int,int)const;
 std::optional<EncryptedDemStitch> tryLoad(int,int,int)const; bool has(int,int,int)const;
 std::vector<StitchAddress> list(int level=-1)const; LevelEncryptStats encryptLevel(int);
 void clearCache()const; bool compress()const{return compress_;}
private:
 using SplitKey=std::tuple<int,int,int>; bool stitchHasData(int,int)const;
 std::pair<std::size_t,std::size_t> validShape(int,int)const; std::vector<double> stitchValues(int,int)const;
 std::filesystem::path findPack(int,int,int,int*)const,packPath(int,int,int,int)const;
 CipherVector createPack(int,int,int,int*)const; std::pair<CipherVector,CipherVector> splitPack16(const CipherVector&)const;
 const IDtmProvider& dtm_;const IHeBackend& he_;std::filesystem::path directory_;bool compress_;
 mutable std::map<SplitKey,std::pair<CipherVector,CipherVector>> split_cache_;
};
}

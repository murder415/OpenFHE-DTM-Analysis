#include "openfhe_dtm/EncryptedDem.hpp"
#include "openfhe_dtm/OpenFheBackend.hpp"
#include "openfhe_dtm/SlopeGrid.hpp"
#include "openfhe_dtm/TdbDataset.hpp"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
double expected(const openfhe_dtm::IDtmProvider&d,int x0,int y0){auto t=d.loadTile(x0,y0);return t.elevations_m[0]/openfhe_dtm::kDemScale;}
void check(const openfhe_dtm::IHeBackend&he,const openfhe_dtm::IDtmProvider&dtm,
           const openfhe_dtm::EncryptedDemStitch&s){auto v=he.decrypt(s.ciphertext).values;double want=expected(dtm,s.x0,s.y0),got=v[s.slot_offset];if(std::abs(got-want)>1e-6)throw std::runtime_error("packed stitch value mismatch at "+std::to_string(s.x0)+","+std::to_string(s.y0)+" got="+std::to_string(got)+" expected="+std::to_string(want));}
}
int main(int argc,char**argv){using namespace openfhe_dtm;if(argc!=3){std::cerr<<"usage: packed_store_validation <dem.tdb> <output-dir>\n";return 2;}try{
 std::filesystem::path root(argv[2]);OpenFheBackend he(kDemPackSlots,root/"keys");auto uk=loadValidationTdbDataset(argv[1]);
 std::vector<std::complex<double>> probe(kDemPackSlots);probe[0]={1,2};auto probe_back=he.decryptComplex(he.encrypt(he.encodeComplex(probe)));std::cout<<"complex probe="<<probe_back[0]<<"\n";
 EncryptedDemStore normal(uk,he,root/"pack",false);auto ns=normal.encryptLevel(0);check(he,uk,normal.load(0,0,0));
 const std::size_t stitches_x=(uk.widthCells()+127)/128,stitches_y=(uk.heightCells()+127)/128;
 const int last_x=static_cast<int>(2*(stitches_x-1)),last_y=static_cast<int>(2*(stitches_y-1));
 if(!normal.has(0,last_x,last_y)||normal.list(0).size()!=stitches_x*stitches_y)
  throw std::runtime_error("2-stitch population mask failed");
 RasterGrid grid;grid.origin={-1.47,53.38};grid.width=256;grid.height=256;grid.cell_size_m=1;grid.nodata=-9999;
 grid.elevations_m.resize(256*256);for(std::size_t y=0;y<256;++y)for(std::size_t x=0;x<256;++x)grid.elevations_m[y*256+x]=100+0.1*x+0.2*y;
 RasterDtmProvider synthetic(DtmMetadata{"pack16 fixture","analytic","local","m",1,-9999},std::move(grid));
 EncryptedDemStore compressed(synthetic,he,root/"pack16",true);auto cs=compressed.encryptLevel(0);
 auto raw=he.decryptComplex(he.loadCipher(root/"pack16"/"pack16_0_0_0_p15.ct"));
 std::cout<<"pack16 raw slot0="<<raw[0].real()<<"+"<<raw[0].imag()<<"i\n";
 for(int y:{0,2})for(int x:{0,2})check(he,synthetic,compressed.load(0,x,y));
 if(compressed.list(0).size()!=4)throw std::runtime_error("pack16 population/list failed");
 compressed.clearCache();
 check(he,synthetic,compressed.load(0,2,2));
 SlopeGridAnalysis slope(he);auto right=compressed.load(0,2,0);auto encrypted_slope=slope.computeEncrypted(right,1,1);
 auto right_actual=slope.decrypt(encrypted_slope),right_expected=SlopeGridAnalysis::originalApproximationPlain(synthetic,2,0,1,1);
 auto parity=SlopeGridAnalysis::accuracy(right_actual,right_expected);
 std::cout<<"right-half slope first actual/expected="<<right_actual[0]<<"/"<<right_expected[0]<<"\n";
 std::cout<<"right-half slope parity mae="<<parity.mae_degrees<<" max="<<parity.max_abs_error_degrees<<"\n";
 if(parity.mae_degrees>1e-3)throw std::runtime_error("right-half slope offset/mask parity failed");
 std::cout<<"pack: total="<<ns.total<<" encrypted="<<ns.encrypted<<" skipped="<<ns.skipped<<" mask/list="<<normal.list(0).size()<<"\n";
 std::cout<<"pack16: total="<<cs.total<<" encrypted="<<cs.encrypted<<" skipped="<<cs.skipped<<" split-stitches="<<compressed.list(0).size()<<"\n";
 std::cout<<"pack offsets: left="<<packSlotOffset(0)<<" right="<<packSlotOffset(2)<<"; real/imag split verified\n";return 0;
 }catch(const std::exception&e){std::cerr<<"packed store validation failed: "<<e.what()<<"\n";return 1;}}

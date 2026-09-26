# OpenFHE DTM Analysis

OpenFHE CKKS(32,768 슬롯)를 사용하는 **C++17 암호문 지형 분석 코드**다. Python은 필요하지 않다.

- 암호화한 **패킹 DEM**을 회전해 고도를 보간
- 시곡면 가용 높이
- 계획고 H와 굴착각 φ를 입력받는 **경사면 반영 순 토공량**

원본 DEM과 비밀키는 발주처가 보유한다. 발주처가 DEM을 한 번 패킹해 암호화하면, 서버는 암호문을 회전하고 평문 가중치를 곱해 계산한다. 마지막 결과만 복호한다.

## 고도 보간 (`PackedElevation`)

64×64 타일 네 장을 128×128 묶음(stitch)으로 합쳐, 두 묶음을 한 암호문(32,768 슬롯)에 담는다. 한 슬롯에 DEM 셀 하나가 들어간다. 묶음은 타일 하나씩 겹치므로 어떤 점의 네 모서리도 항상 같은 암호문 안에 있다.

점 p가 셀 (x0, y0)에 있을 때 서버는 그 셀의 슬롯 t와 셀 안 상대 위치 (dx, dy)만 평문으로 안다.

```text
c10 = Rot(C, 1), c01 = Rot(C, 128), c11 = Rot(C, 129)   # 오른쪽·위·대각 이웃을 슬롯 t로
c_z = C⊙w00 + c10⊙w10 + c01⊙w01 + c11⊙w11              # w는 슬롯 t에만 쌍선형 가중치
```

회전은 묶음 암호문마다 3번이며 그 묶음의 모든 점이 함께 쓴다. 표본 i의 값은 슬롯 t_i에 있다. 쓰는 슬롯이 다른 묶음 결과끼리 더해 복호 횟수를 줄인다.

## 시곡면 (`sightSurfacePacked`)

기준점 b에서 방향 q로 거리 L을 간격 Δ로 등분한 표본의 땅 높이 c_g를 보간한다. 기준점 슬롯만 1인 마스크를 곱하고 모든 슬롯을 더해(EvalSum) z_b를 모든 슬롯에 복사한다. 여기에 d_i·tanθ를 더한 기준면에서 c_g를 빼고, 마지막에 한 번 복호해 음수를 0으로 바꾼다. 지형의 가림은 판정하지 않는다.

## 경사면 반영 순 토공량 (`cutSlopeEarthworkPacked`)

영역 안은 계획고 H로 평평하게 만들고, 영역 경계선에서 바깥쪽으로 거리 W 안의 외곽 영역에는 굴착각 φ의 경사면을 둔다. 꼭짓점 p_i에서 경계선까지의 최단 거리를 ρ_i라 하면 경사면 높이는 s_i = H + ρ_i·tanφ다.

두 영역을 한 변이 Δ_v 이하인 셀로 나누고, 중심이 두 영역 안에 있는 셀의 꼭짓점을 표본으로 쓴다. 셀 하나의 넓이 A_c는 네 꼭짓점이 A_c / 4씩 나눠 맡는다(점고법). 영역 안 가중치는 A, 외곽 가중치는 A′다.

```text
c_z ← InterpEnc(P)
s   ← H + ρ·tanφ
c_δ ← Enc(Encode0(H,…,H)) − c_z
c_u ← EvalSum(c_δ ⊙ Encode0(A))            # 영역 안 성토−절토 합계
c_e ← (c_z − Encode0(s)) ⊙ Encode0(A′)     # 외곽 꼭짓점별 (땅−경사면)×면적
u ← Dec(c_u)[1], e ← Dec(c_e)
V ← u − Σ max(0, e_i)                       # 양수(경사면 위로 솟은 땅)만 굴착량
```

순 토공량은 **성토 − 절토**다. H·φ·W·좌표·가중치는 공개 평문이고 DEM 고도가 암호화 대상이다. 외곽 영역은 깎는 경사면만 계산하며 성토 비탈은 다루지 않는다.

## 구성

| 경로 | 용도 |
| --- | --- |
| `source/include/openfhe_dtm/PackedElevation.hpp`, `source/src/elevation/PackedElevation.cpp` | 패킹 DEM 암호화, 회전 보간, 시곡면 |
| `source/include/openfhe_dtm/CutSlope.hpp`, `source/src/earthwork/CutSlope.cpp` | 경사면 반영 순 토공량 |
| `source/include/openfhe_dtm/EncryptedDem.hpp`, `source/src/dem/EncryptedDem.cpp` | 패킹 암호문 저장·복원 |
| `source/include/openfhe_dtm/OpenFheBackend.hpp`, `source/src/OpenFheBackend.cpp` | OpenFHE CKKS 백엔드(회전 키 ±1…±128, 129) |
| `source/include/openfhe_dtm/TdbDataset.hpp`, `DtmProvider.hpp` | DEM 입력 |
| `analysis_demo.cpp` | TDB 입력으로 시곡면·순 토공량 실행 |
| `tests/cut_slope_validation.cpp` | 합성 경사 평면에서 보간·시곡면·토공량을 독립 평문과 비교, 잘못된 입력 거부 |
| `source/example/make_validation_tdb.cpp`, `convert_neh_to_validation_tdb.cpp` | 합성 DEM 생성, TWD97 N/E/H CSV 변환 |

## 빌드

검증 환경은 Windows 10, GCC 16.1.0 MinGW, C++17, OpenFHE 1.5.1이다. CMake 3.20 이상과 컴파일러 ABI가 호환되는 OpenFHE 개발 설치가 필요하다.

```sh
cmake -S . -B build -DOpenFHE_DIR="<OpenFHEConfig.cmake 디렉터리>" -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

## 실행

```sh
build/make_validation_tdb.exe sample.tdb
build/analysis_demo.exe sample.tdb            # H=14.76 m, 굴착각 45°
build/analysis_demo.exe "<dem.tdb>" 20 30     # H=20 m, 굴착각 30°
```

DEM 안에 사각형 영역과 시선을 만들어 시곡면과 순 토공량을 계산하고, 사용한 패킹 암호문 수·회전 횟수·복호 횟수를 출력한다. 대만 DEM·키·실행 파일은 저장소에 포함하지 않는다. 현재 단일 프로세스 구현이며 서버 간 비밀키·입력 분리는 별도 구현 대상이다.

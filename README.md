# OpenFHE DTM Analysis

OpenFHE CKKS를 사용하는 **C++17 지형 분석 코드**다. Python은 필요하지 않다.

- 고도 보간과 시곡면 가용 높이 계산
- 기준 높이 H와 평면 기울기 θ를 입력받는 **평면 순 토공량** 계산
- 암호문 DEM 타일 패킹·저장·복원 및 격자 경사도 라이브러리

현재 논문의 토공량에 대응하는 API는 `earthworkTheta()`이고 실행 예제는 **`theta_plane_demo`**다. 영역의 첫 점에서 둘째 점으로 향하는 방향으로 평면을 기울이며, 순 토공량은 **성토량−절토량**이다. 이전의 영역 안쪽 경계 사면 모델 `earthworkSlope()`는 별도 API로 보존되어 있다.

## 평면 토공량의 계산

첫 점을 원점으로 하는 동·북 좌표에서 첫 변의 단위 방향을 `(u_x,u_y)`라 한다.

```text
d_i = x_i*u_x + y_i*u_y
s_i = H + d_i*tan(theta)
V   = sum((s_i - z_i)*A_c)
```

`d_i`는 그 방향으로 잰 부호 있는 투영거리다. H는 첫 점에서의 평면 높이이며, θ는 수평면에 대한 기울기다. API·명령행에서는 각도를 **도** 단위로 입력하고 내부에서 라디안으로 바꾼다. `-90<θ<90`이며 θ=0이면 수평면이다. 양수 θ에서는 첫 점→둘째 점 방향으로 높아지고 음수에서는 낮아진다. 첫 두 점의 순서를 바꾸면 방향과 H의 기준점도 바뀐다.

영역의 외접 사각형에 일정 간격의 정사각형 격자를 만들고, 중심이 영역 안에 있는 셀을 선택한다. 셀 중심에서 고도를 보간하고, 각 높이 차에 `A_c=표본 간격²`를 곱한다. 경계 셀도 선택되면 전체 면적을 사용한다. 모든 표본을 한 묶음으로 처리하므로 `1<=n<=S`여야 한다. 표본이 슬롯보다 많으면 오류를 반환한다.

```text
[순 토공량] 보간 암호문 c_z 사용
6: s[1,…,n] ← H+d*tan(theta)
7: A_i ← A_c (i≤n), A_i ← 0 (n<i≤S)
8: c_delta ← Encode_0(s)−c_z
9: c_v ← c_delta⊙Encode_0(A)
10: c_V ← EvalSum(c_v,S)
11: V ← Dec_sk(c_V)[1]
```

공개 설계고에서 암호화된 지반고를 빼며 결과는 암호문이다. `plainSub()`는 OpenFHE의 `EvalSub(평문, 암호문)`을 호출한다. 양의 면적을 곱하고 암호문 안의 값을 합산한 뒤 **최종 합계만 한 번 복호**한다. 빈 슬롯은 면적 가중치 0을 곱해 합산에서 제외한다. 성토는 양수, 절토는 음수이므로 최종 값은 성토−절토다.

좌표·θ·H·설계고·면적은 공개 평문이며 DEM 고도가 암호화 대상이다. 이 API는 순량만 계산하고 개별 절토량·성토량이나 지반고를 출력하지 않는다. 같은 DEM·표본·면적으로 구한 평문 결과와의 수치 차이가 실제 지형에 대한 측량 정확도를 뜻하지는 않는다.

## 굴착 경사면 토공량 (CutSlope)

`openfhe_dtm/CutSlope.hpp`의 `cutslope::cutSlopeEarthwork`는 계획고 H로 영역을 평탄화하고, 영역 밖 공개 폭 W의 띠에 굴착각 θ의 경사면을 둔다.

- 영역과 띠 위에 간격 이하의 격자를 만들고, 셀 면적 A_c/4를 공유 꼭짓점에 더한 가중치(영역 안 A, 띠 A′)를 쓴다(점고법).
- 영역 밖 꼭짓점은 경계선까지의 최단거리 d 하나로 경사면 높이 s = H + d·tanθ를 정한다.
- 암호문으로 Σ(H−z)A를 슬롯 합산하고, (z−s)A′ 벡터를 만든다. 묶음마다 두 암호문만 복호한다.
- 순량 = 영역 안(성토−절토) − 경사면 굴착(양수만 합산).

H·θ·W·좌표·가중치는 공개 평문이고, 지반고는 보간 암호문이다. `tests/cut_slope_validation.cpp`는 합성 경사 평면에서 독립 평문 합계와 비교하고 잘못된 입력 거부를 검사한다.

## 구성과 API

| 경로 또는 API | 용도 |
| --- | --- |
| `source/include/openfhe_dtm/EarthworkTheta.hpp` | 평면 θ 입력과 암호문·최종 결과 자료형 |
| `source/src/earthwork/ThetaPlane.cpp` | 현재 논문의 θ 평면 토공량 구현 |
| `theta_plane_demo.cpp` | θ 평면 토공량 실행 예제 |
| `tests/theta_plane_validation.cpp` | 실제 OpenFHE와 독립 평문 기준, 연산·복호 시점 검사 |
| `source/example/make_validation_tdb.cpp` | 합성 DEM 입력 생성 |
| `source/example/convert_neh_to_validation_tdb.cpp` | TWD97 N/E/H CSV 입력 변환 |

`earthworkThetaEncrypted()`는 최종 합계 암호문을 반환하고, `finalizeEarthworkTheta()`가 이를 한 번 복호한다. `earthworkTheta()`는 두 단계를 연속 실행한다. 고도 보간의 `elevationEncrypted()`는 암호문을 반환한다. 시곡면의 `sightSurfaceEncrypted()`는 높이 차까지 암호문으로 계산하고, `sightSurface()`는 마지막 높이 차를 복호하여 0과 비교한다. 시곡면은 지형의 가림을 판정하지 않는다.

기존 모델은 다음과 같이 구분한다. **평면 논문의 재현에는 첫 행을 사용한다.**

| API / 예제 | 모델과 복호 시점 | 부호 |
| --- | --- | --- |
| `earthworkTheta()` / `theta_plane_demo` | H·θ 평면, 암호문 합산 후 최종 합계 복호 | 성토−절토 |
| `earthworkSlope()` / `final_analysis_demo` | 영역 안쪽 경계 사면, 최종 후보 복호 후 min/max·집계 | 성토−절토 |
| 기존 `TerrainAnalysis::earthwork()` | `a_x,a_y` 평면을 사용하는 이전 API | 절토−성토 |

## 빌드

검증 환경은 Windows 10, GCC 16.1.0 MinGW, C++17, OpenFHE 1.5.1이다. CMake 3.20 이상과 컴파일러 ABI가 호환되는 OpenFHE 개발 설치가 필요하다. 저장소 최상위에서 실행한다.

```sh
cmake -S . -B build -DOpenFHE_DIR="<OpenFHEConfig.cmake 디렉터리>" -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

Windows에서 MinGW를 명시하는 예는 다음과 같다. 생성기를 바꾸면 새 빌드 폴더를 사용한다.

```sh
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="<MinGW bin>/g++.exe" -DCMAKE_MAKE_PROGRAM="<MinGW bin>/mingw32-make.exe" -DOpenFHE_DIR="<OpenFHE CMake 디렉터리>" -DCMAKE_BUILD_TYPE=Release
```

## 실행과 검증

```sh
build/make_validation_tdb.exe sample.tdb
build/theta_plane_demo.exe sample.tdb 100 15 10
build/theta_plane_demo.exe "<dem.tdb>" 14.76 0 915.75
```

인자는 **텍스트 TDB 경로, 기준 높이 H(m), 평면 기울기 θ(도), 선택적인 표본 간격(m)** 순서다. 실행 예제는 2,048슬롯과 설정 깊이 3을 사용하며, DEM 안의 사각형 구역을 구성한다. 임의 다각형은 `EarthworkThetaQuery::polygon`으로 전달한다. GeoTIFF를 직접 읽거나 별도 타일 패킹 시험을 실행하지 않는다.

```sh
build/convert_neh_to_validation_tdb.exe input.csv output.tdb
build/theta_plane_validation.exe "<대만 DEM.tdb>"
```

기본 CTest는 외부 DEM 없이 13개 수치·연산 추적 사례와 17개 잘못된 입력을 검사한다. 선택적인 대만 TDB 인자를 주면 기존 144점 배치의 θ=0°, 10° 비교를 추가해 15개 수치·추적 사례가 된다. 테스트는 독립 `long double` 평문 기준과 비교하며, 네 꼭짓점 벡터의 초기 암호화 4회, 중간 복호 0회, 암호문 합산 1회, 마지막 복호 1회를 확인한다. 정상·역방향·음의 θ·부호 있는 거리·빈 슬롯·슬롯 초과를 포함한다.

이 테스트는 현재 `Encode_0(s)−c_z`와 양의 면적을 사용하는 구현을 검증한다. 논문의 기존 4회 측정은 수학적으로 같은 `(z−s)×(−A)`를 사용한 앞선 구현 기록이다. 그 기록을 현재 소스의 4회 재측정 결과로 표시하지 않으며, 기본 검증 실행은 논문의 반복 성능 측정을 대체하지 않는다.

대만 DEM·키·실행 파일·Python 도구는 저장소에 포함하지 않는다. Linux 계열에서는 실행 파일의 `.exe`를 생략한다. 다른 도구 모음에서의 전체 빌드는 확인하지 않았다. 현재 단일 프로세스 구현이며 서버 간 비밀키·입력 분리는 별도 구현 대상이다.

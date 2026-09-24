# OpenFHE DTM Analysis

CKKS를 이용한 고도 보간, 시곡면 가용 높이, 계획고 **H**와 굴착 사면각 **θ**를 입력받는 토공량 분석 연구 코드다. OpenFHE 1.5.1에서 검증했다. 시곡면과 토공량은 중간 복호 없이 후보를 계산하고, 마지막 비교 단계에서 복호한다.

OpenFHE CKKS terrain-analysis prototype with encrypted interpolation and sightline differences, plus boundary-constrained earthwork controlled by target height and slope angle. Decryption occurs at the final comparison stage; min/max and earthwork aggregation are plaintext owner-side operations.

## 계산 범위

| 기능 | 암호문 계산 | 최종 평문 처리 |
| --- | --- | --- |
| 고도 조회 | 네 모서리 고도의 쌍선형 보간 | 조회 고도 복호 |
| 시곡면 | 보간 → 기준면 → 기준면과 지반고의 차 | 최종 차분 복호 후 max(0, 차이) |
| H·θ 토공량 | 계획고·경계 사면과 지반고의 차, 면적 가중 후보 | 모든 후보 생성 후 복호 → min/max → 절토·성토 분리·집계 |

토공량의 경계 사면은 선택 영역 **안쪽**에 남긴다. 경계점 b의 고도 z_b와 수평거리 d_ib에서 상한은 `min_b(z_b+d_ib tanθ)`, 하한은 `max_b(z_b−d_ib tanθ)`이고 완성면은 `max(하한,min(H,상한))`이다. 새 순 토공량은 **성토−절토**다. 원래 `earthwork()`의 방향별 설계 평면과 절토−성토 정의는 비교용으로 보존했다.

암호문 보간 API는 `elevationEncrypted()`, 최종 차분 시곡면은 `sightSurfaceEncrypted()`/`sightSurface()`, 새 토공량은 `earthworkSlopeEncrypted()`/`finalizeEarthworkSlope()`/`earthworkSlope()`다. 타입과 입력 설명은 [Analysis.hpp](source/include/openfhe_dtm/Analysis.hpp), 사면 구현은 [Excavation.cpp](source/src/earthwork/Excavation.cpp)에 있다.

## 빌드

검증 환경은 Windows 10, GCC 16.1.0 MinGW, C++17 Release, OpenFHE 1.5.1이다. CMake 3.20 이상과 ABI가 호환되는 OpenFHE 개발 설치가 필요하다. 아래 경로를 자신의 설치에 맞게 지정한다. `source/`가 아닌 저장소 최상위에서 구성해야 새 검증 프로그램과 데모가 포함된다.

```sh
cmake -S . -B build -DOpenFHE_DIR="<OpenFHEConfig.cmake 디렉터리>" -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

Windows에서 MinGW 도구를 명시하려면 다음 형태를 사용할 수 있다.

```sh
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="<MinGW bin>/g++.exe" -DCMAKE_MAKE_PROGRAM="<MinGW bin>/mingw32-make.exe" -DOpenFHE_DIR="<OpenFHE CMake 디렉터리>" -DCMAKE_BUILD_TYPE=Release
```

컴파일러나 생성기를 바꿀 때는 새 빌드 폴더를 사용한다. `quadmath`를 사용하는 진단 프로그램이 있어 검증된 GCC 환경을 권장한다. 다른 도구 모음에서의 전체 빌드는 확인하지 않았다.

## 실행

```sh
build/final_analysis_demo.exe "<dem.tdb>" 14.76 45 915.75
```

인자는 **텍스트 TDB 경로, 절대 계획고 H(m), 사면각 θ(도), 선택적인 요청 격자 간격(m)** 순서다. 각도 범위는 `0<θ<90`이다. GeoTIFF를 직접 읽지 않는다. CLI는 DEM 안의 재현 가능한 사각형 구역을 사용하며, 임의 다각형은 `EarthworkSlopeQuery::polygon`으로 전달한다. 새 토공량 결과의 `envelope_valid`가 false이면 유효한 체적으로 사용하지 않는다.

데이터 없이 CTest의 합성 지형 검증을 실행할 수 있다. 합성 텍스트 TDB 생성기는 `build/upstream/make_validation_tdb.exe <output.tdb>`이다. 이 합성 데이터의 결과를 대만 실험 결과로 해석하지 않는다.

## 검증과 기록

새 시곡면·토공량 검증은 별도 CLI다. **CTest만 실행하면 아래 두 새 검증까지 실행되는 것은 아니다.**

```sh
build/sightline_validation.exe "<taiwan.tdb>" "runs/sightline" 4
build/earthwork_slope_validation.exe "<taiwan.tdb>" "runs/earthwork" 4
```

Linux 계열에서는 `.exe`를 생략한다. 재실행 기록은 `runs/`에 저장하여 제공된 `logs/`를 덮어쓰지 않는 예시다. Windows 프로세스 시간·메모리 측정 도구는 `benchmark_process.py`이며 인자는 `--help`로 확인한다.

| 검증 | 기록 |
| --- | --- |
| 시곡면 12조건×4회, 48건 통과 | [시곡면 요약](logs/sightline/summary.json) |
| 토공량 16조건×4회, 64건 통과 | [토공량 요약](logs/earthwork/summary.json) |
| 잘못된 입력 11종×4회 거부 | `logs/earthwork/run_*/results.json` |
| 기존 회귀 테스트 6/6 | [CTest 기록](logs/ctest/1_ctest_final.txt) |
| 분석·패킹의 독립 프로세스 각 4회 | `logs/benchmarks/` |

실행 순서 추적으로 모든 후보 생성 전 중간 복호가 없고 첫 복호 이후 암호문 연산이 없는지 확인했다. 원본 DEM과 별도 기하 계산의 독립 평문 기준을 비교했다. 토공량은 배치마다 `2×경계점 수+1`개의 암호문을 마지막 단계에서 복호하며, 복호 호출이 한 번이라는 뜻은 아니다.

대만 H=14.76m·θ=45° 시험은 169셀·196적분 정점·52경계 검사점, 2,048슬롯·설정 깊이3이다. 4회 평균 순량(성토−절토)은 약 −7.419041×10⁸m³이고, 절토·성토·순량의 독립 평문 대비 최대 절대차는 약 5.44×10⁻⁵m³다. 같은 이산 모델 간 계산 차이이며 실측 정확도는 아니다.

## 대만 DEM 입력

대만 시험에 사용한 원본 fixture는 이 공개 저장소에 포함하지 않는다. 정확히 같은 기록을 재현하려면 [dataset_manifest.json](dataset_manifest.json)의 SHA-256과 일치하는 텍스트 TDB를 별도로 준비한다. 이미 제공된 로컬 배포 패키지에는 `data/taiwan.tdb`가 있다.

출처는 Copernicus DEM GLO-30 Public 2021, Taipei/New Taipei의 1000×1000 추출 격자다. [제품 DOI](https://doi.org/10.5270/ESA-c5d3d65), [공개 데이터 안내](https://registry.opendata.aws/copernicus-dem/)에서 원 제품을 확인할 수 있다. 다른 데이터 파일로 실행한 결과는 제공된 대만 수치와 일치하지 않을 수 있다.

## 적용 범위와 공개 파일

현재 단일 프로세스 prototype이며 backend가 비밀키를 함께 보유한다. 실제 서버 간 입력·키 분리, 접근 통제는 구현하지 않았다. 복호 후보에서 고도를 역산할 수 있으므로 최종 비교는 원본을 볼 권한이 있는 소유자 내부에서 수행한다. θ·좌표·거리·격자·가중치는 공개다. 시곡면은 지형 가림을 판정하지 않는다.

사면은 유한 경계 표본의 포락면이고 적분은 셀 정점 가중 근사다. 경계 부분 셀·혼합 절성토 셀의 정확한 분할, DEM 측량오차, 해상도에 따른 수렴, 실제 시공량은 검증하지 않았다.

검증된 C++·Python·CMake 파일은 로컬 배포본과 바이트 단위로 같다. 공개 로그의 개인 PC 절대경로만 `<LOCAL_...>` 표기로 치환했으며 측정값과 소스·실행파일·데이터 해시는 보존했다. 현재 공개 파일의 해시는 `public_manifest.json`에 기록한다. 문서별 AI 검토 메모, 논문 파일, 실행파일, 키·캐시·DEM 원본은 업로드하지 않는다.

별도의 저장소 전체 재사용 라이선스는 지정하지 않았다. 제3자 코드에 포함된 저작권·라이선스 고지는 그대로 유지했다. OpenFHE 예제의 BSD 2-Clause 고지는 [해당 소스](source/example/openfhe_official_ckks.cpp)에 있다.

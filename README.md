# OpenFHE DTM Analysis

OpenFHE CKKS를 사용하는 **C++17 지형 분석 코드**다. Python은 필요하지 않다.

- 고도 보간과 시곡면 가용 높이 계산
- 계획고 H와 굴착 사면각 θ를 입력받는 토공량 계산
- 암호문 DEM 타일 패킹·저장·복원 및 격자 경사도 라이브러리

시곡면은 최종 높이 차이를 복호한 뒤 0과 비교한다. 토공량은 모든 사면 후보를 암호문으로 계산한 뒤 복호하여 min/max 비교와 집계를 수행한다. 새 토공량은 선택 영역 안쪽에 사면을 남기며 **순 토공량=성토량−절토량**이다.

## 구성

| 경로 | 용도 |
| --- | --- |
| `source/include/openfhe_dtm/` | 분석 API와 자료형 |
| `source/src/` | 고도·시곡면·토공량·패킹·경사도 구현 |
| `final_analysis_demo.cpp` | 시곡면과 H·θ 토공량 실행 예제 |
| `source/example/make_validation_tdb.cpp` | 합성 DEM 입력 생성 |
| `source/example/convert_neh_to_validation_tdb.cpp` | TWD97 N/E/H CSV 입력 변환 |

## 빌드

검증 환경: Windows 10, GCC 16.1.0 MinGW, C++17, OpenFHE 1.5.1. CMake 3.20 이상과 컴파일러 ABI가 호환되는 OpenFHE 개발 설치가 필요하다. 저장소 최상위에서 실행한다.

```sh
cmake -S . -B build -DOpenFHE_DIR="<OpenFHEConfig.cmake 디렉터리>" -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
```

Windows에서 MinGW를 명시하는 구성 명령은 다음과 같다. 생성기를 바꾸면 새 빌드 폴더를 사용한다.

```sh
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="<MinGW bin>/g++.exe" -DCMAKE_MAKE_PROGRAM="<MinGW bin>/mingw32-make.exe" -DOpenFHE_DIR="<OpenFHE CMake 디렉터리>" -DCMAKE_BUILD_TYPE=Release
```

## 실행

```sh
build/final_analysis_demo.exe "<dem.tdb>" 14.76 45 915.75
```

인자는 **텍스트 TDB 경로, 절대 계획고 H(m), 사면각 θ(도), 선택적인 요청 격자 간격(m)** 순서다. `0<θ<90`이며 GeoTIFF를 직접 읽지는 않는다. 예제는 DEM 안의 사각형 구역을 사용한다. 임의 다각형은 [EarthworkSlopeQuery](source/include/openfhe_dtm/Analysis.hpp)의 `polygon`으로 지정한다. 예제 프로그램은 별도 타일 패킹 시험을 실행하지 않는다.

합성 입력으로 실행하려면 다음과 같이 사용할 수 있다. 대만 실험 데이터는 포함하지 않는다.

```sh
build/make_validation_tdb.exe sample.tdb
build/final_analysis_demo.exe sample.tdb 100 45 10
```

TWD97 N/E/H CSV 입력이 있으면 다음 변환기를 사용할 수 있다.

```sh
build/convert_neh_to_validation_tdb.exe input.csv output.tdb
```

Linux 계열에서는 실행 파일의 `.exe`를 생략한다. 다른 도구 모음에서의 전체 빌드는 확인하지 않았다.

## API와 적용 범위

`elevationEncrypted()`는 암호문 보간 결과를 반환한다. `sightSurface()`와 `earthworkSlope()`는 최종 복호·비교까지 실행한다. 암호문 단계만 필요하면 `sightSurfaceEncrypted()` 또는 `earthworkSlopeEncrypted()`를 사용한다. 후자의 결과는 `finalizeEarthworkSlope()`에 전달한다. `envelope_valid=false`이면 정상 토공량으로 사용하지 않는다. 기존 평면 토공량 `earthwork()`는 별도 API이며 부호 정의가 절토−성토다.

현재 단일 프로세스 구현이며 실제 서버 간 비밀키·입력 분리는 구현하지 않았다. 최종 후보는 원본을 볼 권한이 있는 소유자 내부에서 복호한다. 좌표·거리·격자·각도·면적 가중치는 공개다. 시곡면은 지형의 가림을 판정하지 않으며 토공량은 격자와 경계 표본을 이용한 근사다.

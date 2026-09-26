cmake_minimum_required(VERSION 3.20)
project(openfhe_dtm_core LANGUAGES CXX)

find_package(OpenFHE CONFIG REQUIRED)
add_library(openfhe_dtm_core STATIC
    src/DtmProvider.cpp
    src/OpenFheBackend.cpp
    src/TdbDataset.cpp
    src/elevation/PackedElevation.cpp
    src/earthwork/CutSlope.cpp
    src/dem/EncryptedDem.cpp
)
target_compile_features(openfhe_dtm_core PUBLIC cxx_std_17)
target_compile_definitions(openfhe_dtm_core PUBLIC OPENFHE_DTM_ENABLE_OPENFHE=1)
target_include_directories(openfhe_dtm_core PUBLIC
    include
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
)
separate_arguments(OpenFHE_CXX_FLAGS_LIST NATIVE_COMMAND "${OpenFHE_CXX_FLAGS}")
target_compile_options(openfhe_dtm_core PUBLIC ${OpenFHE_CXX_FLAGS_LIST})
target_link_libraries(openfhe_dtm_core PUBLIC ${OpenFHE_LIBRARIES})

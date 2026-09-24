cmake_minimum_required(VERSION 3.20)

project(openfhe_dtm_port LANGUAGES CXX)

option(OPENFHE_DTM_ENABLE_OPENFHE "Build with OpenFHE" ON)

include(CTest)

add_library(openfhe_dtm_core STATIC
    src/Analysis.cpp
    src/DtmProvider.cpp
    src/OpenFheBackend.cpp
    src/TdbDataset.cpp
    src/elevation/Elevation.cpp
    src/earthwork/Excavation.cpp
    src/dem/EncryptedDem.cpp
    src/slope/SlopeGrid.cpp
)

target_compile_features(openfhe_dtm_core PUBLIC cxx_std_17)
target_include_directories(openfhe_dtm_core PUBLIC include)

if(NOT OPENFHE_DTM_ENABLE_OPENFHE)
    message(FATAL_ERROR "This workspace is OpenFHE-only. Configure with OPENFHE_DTM_ENABLE_OPENFHE=ON.")
endif()

find_package(OpenFHE CONFIG REQUIRED)
separate_arguments(OpenFHE_CXX_FLAGS_LIST NATIVE_COMMAND "${OpenFHE_CXX_FLAGS}")
target_compile_definitions(openfhe_dtm_core PUBLIC OPENFHE_DTM_ENABLE_OPENFHE=1)
target_include_directories(openfhe_dtm_core PUBLIC
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
)
target_compile_options(openfhe_dtm_core PUBLIC ${OpenFHE_CXX_FLAGS_LIST})
target_link_libraries(openfhe_dtm_core PUBLIC ${OpenFHE_LIBRARIES})

add_executable(openfhe_dtm_demo example/demo.cpp)
target_link_libraries(openfhe_dtm_demo PRIVATE openfhe_dtm_core)

add_executable(synthetic_validation example/synthetic_validation.cpp)
target_link_libraries(synthetic_validation PRIVATE openfhe_dtm_core)

add_executable(make_validation_tdb example/make_validation_tdb.cpp)
target_link_libraries(make_validation_tdb PRIVATE openfhe_dtm_core)

add_executable(convert_neh_to_validation_tdb example/convert_neh_to_validation_tdb.cpp)
target_link_libraries(convert_neh_to_validation_tdb PRIVATE openfhe_dtm_core)

add_executable(slope_grid_validation example/slope_grid_validation.cpp)
target_link_libraries(slope_grid_validation PRIVATE openfhe_dtm_core)
add_executable(packed_store_validation example/packed_store_validation.cpp)
target_link_libraries(packed_store_validation PRIVATE openfhe_dtm_core)
add_executable(bootstrap_validation example/bootstrap_validation.cpp)
target_link_libraries(bootstrap_validation PRIVATE openfhe_dtm_core)
add_executable(manifest_validation example/manifest_validation.cpp)
target_link_libraries(manifest_validation PRIVATE openfhe_dtm_core)
add_executable(experiments_validation example/experiments_validation.cpp)
target_link_libraries(experiments_validation PRIVATE openfhe_dtm_core)
add_executable(plain_baseline example/plain_baseline.cpp)
target_link_libraries(plain_baseline PRIVATE openfhe_dtm_core)
add_executable(ckks_precision_probe example/ckks_precision_probe.cpp)
target_link_libraries(ckks_precision_probe PRIVATE openfhe_dtm_core)
add_executable(ckks_honest_probe example/ckks_honest_probe.cpp)
target_link_libraries(ckks_honest_probe PRIVATE openfhe_dtm_core quadmath)
add_executable(ckks_stage_probe example/ckks_stage_probe.cpp)
target_link_libraries(ckks_stage_probe PRIVATE openfhe_dtm_core quadmath)
add_executable(openfhe_official_ckks example/openfhe_official_ckks.cpp)
target_link_libraries(openfhe_official_ckks PRIVATE openfhe_dtm_core)
add_executable(ckks_encode_probe example/ckks_encode_probe.cpp)
target_link_libraries(ckks_encode_probe PRIVATE openfhe_dtm_core quadmath)

# Windows does not search the OpenFHE or MinGW installation directories when an
# executable is launched directly from Explorer.  Keep every required runtime
# DLL beside each generated executable so build-tree programs work without a
# manually modified PATH.
if(WIN32)
    get_filename_component(OPENFHE_DTM_COMPILER_BIN "${CMAKE_CXX_COMPILER}" DIRECTORY)
    set(OPENFHE_DTM_RUNTIME_DLLS
        "${OpenFHE_LIBDIR}/libOPENFHEcore.dll"
        "${OpenFHE_LIBDIR}/libOPENFHEpke.dll"
        "${OpenFHE_LIBDIR}/libOPENFHEbinfhe.dll"
        "${OPENFHE_DTM_COMPILER_BIN}/libgcc_s_seh-1.dll"
        "${OPENFHE_DTM_COMPILER_BIN}/libstdc++-6.dll"
        "${OPENFHE_DTM_COMPILER_BIN}/libwinpthread-1.dll"
    )
    foreach(runtime_dll IN LISTS OPENFHE_DTM_RUNTIME_DLLS)
        if(NOT EXISTS "${runtime_dll}")
            message(FATAL_ERROR "Required runtime DLL was not found: ${runtime_dll}")
        endif()
    endforeach()

    function(openfhe_dtm_copy_runtime_dlls target_name)
        foreach(runtime_dll IN LISTS OPENFHE_DTM_RUNTIME_DLLS)
            add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${runtime_dll}" "$<TARGET_FILE_DIR:${target_name}>"
                VERBATIM)
        endforeach()
    endfunction()

    foreach(executable_target
            openfhe_dtm_demo
            synthetic_validation
            make_validation_tdb
            convert_neh_to_validation_tdb
            slope_grid_validation
            packed_store_validation
            bootstrap_validation
            manifest_validation
            experiments_validation
            plain_baseline
            ckks_precision_probe
            ckks_honest_probe
            ckks_stage_probe
            openfhe_official_ckks
            ckks_encode_probe)
        openfhe_dtm_copy_runtime_dlls(${executable_target})
    endforeach()
endif()

if(BUILD_TESTING)
    add_test(NAME synthetic_validation COMMAND synthetic_validation)
    add_test(NAME demo_synthetic COMMAND openfhe_dtm_demo --synthetic)
    add_test(NAME ckks_bootstrap COMMAND bootstrap_validation)
    add_test(NAME key_manifest COMMAND manifest_validation
        "${CMAKE_CURRENT_BINARY_DIR}/manifest-test-keys")
    add_test(NAME converter_rejects_unknown_option
        COMMAND convert_neh_to_validation_tdb input.csv output.tdb --unknown)
    set_tests_properties(converter_rejects_unknown_option PROPERTIES WILL_FAIL TRUE)
    add_test(NAME converter_rejects_bad_bounds
        COMMAND convert_neh_to_validation_tdb input.csv output.tdb --bounds bad 1 2 3)
    set_tests_properties(converter_rejects_bad_bounds PROPERTIES WILL_FAIL TRUE)
endif()

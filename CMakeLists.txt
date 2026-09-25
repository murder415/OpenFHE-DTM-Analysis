cmake_minimum_required(VERSION 3.20)
project(OpenFHE_DTM_Analysis LANGUAGES CXX)

find_package(OpenFHE CONFIG REQUIRED)
add_subdirectory(source)

function(add_dtm_executable target source_file)
    add_executable(${target} ${source_file})
    target_link_libraries(${target} PRIVATE openfhe_dtm_core)
    target_compile_features(${target} PRIVATE cxx_std_17)
    if(WIN32)
        get_filename_component(compiler_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
        foreach(dll
            "${OpenFHE_LIBDIR}/libOPENFHEcore.dll"
            "${OpenFHE_LIBDIR}/libOPENFHEpke.dll"
            "${OpenFHE_LIBDIR}/libOPENFHEbinfhe.dll"
            "${compiler_bin}/libgcc_s_seh-1.dll"
            "${compiler_bin}/libstdc++-6.dll"
            "${compiler_bin}/libwinpthread-1.dll"
            "${compiler_bin}/libgomp-1.dll"
            "${compiler_bin}/libquadmath-0.dll")
            if(EXISTS "${dll}")
                add_custom_command(TARGET ${target} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${dll}" "$<TARGET_FILE_DIR:${target}>" VERBATIM)
            endif()
        endforeach()
    endif()
endfunction()

add_dtm_executable(final_analysis_demo final_analysis_demo.cpp)
add_dtm_executable(make_validation_tdb source/example/make_validation_tdb.cpp)
add_dtm_executable(convert_neh_to_validation_tdb source/example/convert_neh_to_validation_tdb.cpp)

add_dtm_executable(theta_plane_demo theta_plane_demo.cpp)
include(CTest)
if(BUILD_TESTING)
    add_dtm_executable(theta_plane_validation tests/theta_plane_validation.cpp)
    add_test(NAME theta_plane_validation COMMAND theta_plane_validation)
endif()
if(BUILD_TESTING)
    add_dtm_executable(cut_slope_validation tests/cut_slope_validation.cpp)
    add_test(NAME cut_slope_validation COMMAND cut_slope_validation)
endif()

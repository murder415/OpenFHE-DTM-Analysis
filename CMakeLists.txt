cmake_minimum_required(VERSION 3.20)
project(final_decrypt_validation LANGUAGES CXX)
include(CTest)
find_package(OpenFHE CONFIG REQUIRED)
add_subdirectory(source upstream)

function(add_review_validation target source_file)
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

add_review_validation(sightline_validation sightline_validation.cpp)
add_review_validation(final_analysis_demo final_analysis_demo.cpp)
add_review_validation(earthwork_slope_validation earthwork_slope_validation.cpp)

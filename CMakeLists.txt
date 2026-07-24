cmake_minimum_required(VERSION 3.20)
project(vcftools_ng VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(HTSLIB_ROOT "" CACHE PATH "HTSlib installation or build directory")

find_path(
  HTSLIB_INCLUDE_DIR
  NAMES htslib/vcf.h
  HINTS "${HTSLIB_ROOT}" ENV HTSLIB_ROOT
  PATH_SUFFIXES include
  REQUIRED
)

find_library(
  HTSLIB_LIBRARY
  NAMES hts libhts
  HINTS "${HTSLIB_ROOT}" ENV HTSLIB_ROOT
  PATH_SUFFIXES lib
  REQUIRED
)

find_package(Threads REQUIRED)

add_executable(
  vcftools-ng
  src/main.cpp
)

target_include_directories(vcftools-ng PRIVATE "${HTSLIB_INCLUDE_DIR}")
target_link_libraries(vcftools-ng PRIVATE "${HTSLIB_LIBRARY}" Threads::Threads)
target_compile_options(vcftools-ng PRIVATE -Wall -Wextra -Wpedantic)

get_filename_component(HTSLIB_LIBRARY_DIR "${HTSLIB_LIBRARY}" DIRECTORY)
set_target_properties(
  vcftools-ng
  PROPERTIES
    BUILD_RPATH "${HTSLIB_LIBRARY_DIR}"
    INSTALL_RPATH "${HTSLIB_LIBRARY_DIR}"
)

include(CTest)
if(BUILD_TESTING)
  add_test(
    NAME differential_compatibility
    COMMAND
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential.sh"
      "$<TARGET_FILE:vcftools-ng>"
      "${CMAKE_CURRENT_SOURCE_DIR}"
  )
endif()

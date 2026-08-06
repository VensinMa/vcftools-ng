cmake_minimum_required(VERSION 3.20)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" VCFTOOLS_NG_VERSION
     LIMIT_COUNT 1)
project(vcftools_ng VERSION "${VCFTOOLS_NG_VERSION}" LANGUAGES CXX)

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
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
include(GNUInstallDirs)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/version.h"
  @ONLY
)

add_executable(
  vcftools-ng
  src/fast_counts.cpp
  src/input_source.cpp
  src/main.cpp
  src/output_transaction.cpp
  src/run_logger.cpp
)

target_include_directories(
  vcftools-ng
  PRIVATE
    "${HTSLIB_INCLUDE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/generated"
)
target_link_libraries(
  vcftools-ng
  PRIVATE
    "${HTSLIB_LIBRARY}"
    Threads::Threads
    ${LAPACK_LIBRARIES}
    ZLIB::ZLIB
)
target_compile_options(vcftools-ng PRIVATE -Wall -Wextra -Wpedantic)

get_filename_component(HTSLIB_LIBRARY_DIR "${HTSLIB_LIBRARY}" DIRECTORY)
set_target_properties(
  vcftools-ng
  PROPERTIES
    BUILD_RPATH "${HTSLIB_LIBRARY_DIR}"
)

install(
  TARGETS vcftools-ng
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
install(
  FILES LICENSE NOTICE README.md README.zh-CN.md
  DESTINATION "${CMAKE_INSTALL_DOCDIR}"
)

include(CTest)
if(BUILD_TESTING)
  find_program(BCFTOOLS_EXECUTABLE NAMES bcftools)
  add_test(
    NAME help_output
    COMMAND
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/help_output.sh"
      "$<TARGET_FILE:vcftools-ng>"
  )
  if(BCFTOOLS_EXECUTABLE)
    add_test(
      NAME run_logging
      COMMAND
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_logging.sh"
        "$<TARGET_FILE:vcftools-ng>"
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    set_tests_properties(
      run_logging PROPERTIES
      ENVIRONMENT "BCFTOOLS=${BCFTOOLS_EXECUTABLE}"
    )
  else()
    message(STATUS "bcftools not found: run_logging test is disabled")
  endif()

  add_test(
    NAME reliability_regression
    COMMAND
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/reliability_regression.sh"
      "$<TARGET_FILE:vcftools-ng>"
      "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/osmanthus412.23chr_100k.bcf")
    add_test(
      NAME input_backend_differential
      COMMAND
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_backends.sh"
        "$<TARGET_FILE:vcftools-ng>"
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(
      NAME differential_compatibility
      COMMAND
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential.sh"
        "$<TARGET_FILE:vcftools-ng>"
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(
      NAME adaptive_counts
      COMMAND
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/adaptive_counts.sh"
        "$<TARGET_FILE:vcftools-ng>"
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(
      NAME adaptive_index_policy
      COMMAND
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/adaptive_index_policy.sh"
        "$<TARGET_FILE:vcftools-ng>"
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
  else()
    message(STATUS "large local fixtures absent: real-data tests are disabled")
  endif()
endif()

cmake_minimum_required(VERSION 3.20)
set_property(
  DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/VERSION"
)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" VCFTOOLS_NG_VERSION
     LIMIT_COUNT 1)
project(vcftools_ng VERSION "${VCFTOOLS_NG_VERSION}" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(HTSLIB_ROOT "" CACHE PATH "HTSlib installation or build directory")
option(VCFTOOLS_NG_ENABLE_IPO "Enable supported release IPO/LTO" ON)
set(VCFTOOLS_NG_PGO "OFF" CACHE STRING "PGO mode: OFF, GENERATE, or USE")
set_property(CACHE VCFTOOLS_NG_PGO PROPERTY STRINGS OFF GENERATE USE)
set(VCFTOOLS_NG_PGO_DIR "${CMAKE_BINARY_DIR}/pgo" CACHE PATH
    "Directory containing compiler profile data")

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
  src/exact_hwe.cpp
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

if(VCFTOOLS_NG_ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT VCFTOOLS_NG_IPO_SUPPORTED
                      OUTPUT VCFTOOLS_NG_IPO_ERROR)
  if(VCFTOOLS_NG_IPO_SUPPORTED)
    set_property(TARGET vcftools-ng PROPERTY
                 INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set_property(TARGET vcftools-ng PROPERTY
                 INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
    set_property(TARGET vcftools-ng PROPERTY
                 INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL TRUE)
  else()
    message(STATUS "IPO/LTO disabled: ${VCFTOOLS_NG_IPO_ERROR}")
  endif()
endif()

string(TOUPPER "${VCFTOOLS_NG_PGO}" VCFTOOLS_NG_PGO_MODE)
if(NOT VCFTOOLS_NG_PGO_MODE MATCHES "^(OFF|GENERATE|USE)$")
  message(FATAL_ERROR "VCFTOOLS_NG_PGO must be OFF, GENERATE, or USE")
endif()
if(NOT VCFTOOLS_NG_PGO_MODE STREQUAL "OFF")
  file(MAKE_DIRECTORY "${VCFTOOLS_NG_PGO_DIR}")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(VCFTOOLS_NG_PGO_MODE STREQUAL "GENERATE")
      set(VCFTOOLS_NG_PGO_FLAGS "-fprofile-generate=${VCFTOOLS_NG_PGO_DIR}")
    else()
      set(VCFTOOLS_NG_PGO_FLAGS
          "-fprofile-use=${VCFTOOLS_NG_PGO_DIR}" -fprofile-correction)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(VCFTOOLS_NG_PGO_MODE STREQUAL "GENERATE")
      set(VCFTOOLS_NG_PGO_FLAGS
          "-fprofile-instr-generate=${VCFTOOLS_NG_PGO_DIR}/default-%p.profraw")
    else()
      set(VCFTOOLS_NG_PGO_FLAGS
          "-fprofile-instr-use=${VCFTOOLS_NG_PGO_DIR}/merged.profdata")
    endif()
  else()
    message(FATAL_ERROR
            "PGO is currently supported only with GCC and Clang")
  endif()
  target_compile_options(vcftools-ng PRIVATE ${VCFTOOLS_NG_PGO_FLAGS})
  target_link_options(vcftools-ng PRIVATE ${VCFTOOLS_NG_PGO_FLAGS})
endif()

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
  add_executable(
    resource_budget_test
    tests/resource_budget.cpp
    src/input_source.cpp
  )
  target_include_directories(
    resource_budget_test
    PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}/src"
      "${HTSLIB_INCLUDE_DIR}"
  )
  target_link_libraries(
    resource_budget_test PRIVATE "${HTSLIB_LIBRARY}" Threads::Threads
  )
  target_compile_options(
    resource_budget_test PRIVATE -Wall -Wextra -Wpedantic
  )
  add_test(NAME resource_budget COMMAND resource_budget_test)

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
    # ASan/LSan runtimes create their own helper thread, so /proc task counts
    # no longer measure the application's strict budget. Keep this gate on
    # normal Linux builds; sanitizer jobs still exercise all planning,
    # compatibility, and failure paths.
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
       NOT CMAKE_CXX_FLAGS MATCHES "fsanitize")
      add_test(
        NAME thread_budget
        COMMAND
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/thread_budget.sh"
          "$<TARGET_FILE:vcftools-ng>"
          "${BCFTOOLS_EXECUTABLE}"
      )
    endif()
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

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/osmanthus205.gatk.23chr_1k.vcf.gz"
     AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/osmanthus205.gatk.23chr_1k.bcf")
    add_test(
      NAME gatk_fixture_compatibility
      COMMAND
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/gatk_fixture_compatibility.sh"
        "$<TARGET_FILE:vcftools-ng>"
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
  else()
    message(STATUS "GATK 23k fixtures absent: GATK compatibility test is disabled")
  endif()

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

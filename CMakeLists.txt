cmake_minimum_required(VERSION 3.20)
project(Quaternion VERSION 1.0.0 DESCRIPTION "Header-only C++23 quaternions modeled on std::complex" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_library(quaternion INTERFACE)
add_library(Quaternion::quaternion ALIAS quaternion)
target_compile_features(quaternion INTERFACE cxx_std_23)
target_include_directories(quaternion INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

option(QUATERNION_BUILD_EXAMPLES "Build examples" ${PROJECT_IS_TOP_LEVEL})
option(QUATERNION_BUILD_BENCHMARKS "Build the benchmark suite (bench/)" OFF)
include(CTest)

if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
    foreach(t core math io special special2 annexg)
        add_executable(test_${t} tests/test_${t}.cpp)
        target_link_libraries(test_${t} PRIVATE quaternion)
        if(MSVC)
            target_compile_options(test_${t} PRIVATE /W4 /permissive- /utf-8)
        else()
            target_compile_options(test_${t} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Werror)
        endif()
        add_test(NAME ${t} COMMAND test_${t})
    endforeach()
    # Evaluator for the nightly mpmath fuzz job (fuzz/fuzz.py, .github/workflows/fuzz.yml);
    # built with the tests so that it keeps compiling on every toolchain.
    add_executable(quaternion_fuzz_eval fuzz/fuzz_eval.cpp)
    target_link_libraries(quaternion_fuzz_eval PRIVATE quaternion)
    if(NOT MSVC)
        target_compile_options(quaternion_fuzz_eval PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Werror)
    endif()
endif()

if(QUATERNION_BUILD_EXAMPLES)
    add_executable(example_basic examples/basic.cpp)
    target_link_libraries(example_basic PRIVATE quaternion)
    if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
        add_test(NAME example_basic COMMAND example_basic)
    endif()
endif()

if(QUATERNION_BUILD_BENCHMARKS)
    add_executable(quaternion_bench bench/bench.cpp)
    target_link_libraries(quaternion_bench PRIVATE quaternion)
endif()

install(TARGETS quaternion EXPORT QuaternionTargets)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT QuaternionTargets NAMESPACE Quaternion::
        DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Quaternion)
configure_package_config_file(cmake/QuaternionConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/QuaternionConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Quaternion)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/QuaternionConfigVersion.cmake
    COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/QuaternionConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/QuaternionConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Quaternion)

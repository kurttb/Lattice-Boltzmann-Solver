cmake_minimum_required(VERSION 3.12)

project(LBM VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(hello)

# add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../kokkos" "${CMAKE_CURRENT_BINARY_DIR}/kokkos")
# # find_package(Kokkos REQUIRED)

set(Kokkos_ROOT "/home/kurttb/ners570_project/kokkos/build" CACHE PATH "Root directory of Kokkos installation")
set(Kokkos_COMPILE_LAUNCHER "${Kokkos_ROOT}/../bin/kokkos_launch_compiler" CACHE FILEPATH "Path to kokkos_launch_compiler script")
set(Kokkos_NVCC_WRAPPER "${Kokkos_ROOT}/../bin/nvcc_wrapper" CACHE FILEPATH "Path to nvcc_wrapper script")

# Switch to find_package to use the pre-built library
# add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../kokkos" "${CMAKE_CURRENT_BINARY_DIR}/kokkos")
find_package(Kokkos REQUIRED)


add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)

hello(${CMAKE_CURRENT_SOURCE_DIR})

# Define the exectuable
add_executable(lbm_example
    LBM_example.cpp
)



MESSAGE(STATUS "LBM_INCLUDE_DIR=${LBM_INCLUDE_DIR}")

# Declare what directories to include
target_include_directories(lbm_example PUBLIC ${LBM_INCLUDE_DIR})

find_package(Kokkos REQUIRED)

add_executable(kokkos_implementation 
    kokkos_implementation.cpp
    ../src/vtk_writer.cpp
)
target_link_libraries(kokkos_implementation PRIVATE Kokkos::kokkos)
target_include_directories(kokkos_implementation PUBLIC ${LBM_INCLUDE_DIR})

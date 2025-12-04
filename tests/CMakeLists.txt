include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# D2Q9Problem Test
add_executable(D2Q9Problem_test D2Q9Problem_test.cpp)
target_link_libraries(D2Q9Problem_test PRIVATE GTest::gtest_main LBM Kokkos::kokkos)
# We need to include the root include directory so tests can find headers
target_include_directories(D2Q9Problem_test PRIVATE ${CMAKE_SOURCE_DIR}/include)

include(GoogleTest)
gtest_discover_tests(D2Q9Problem_test)

# VTK Writer Test
add_executable(vtk_test vtk_test.cpp)
target_link_libraries(vtk_test PRIVATE GTest::gtest_main LBM)
target_include_directories(vtk_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
gtest_discover_tests(vtk_test)

# Kernels Test
add_executable(Kernels_test Kernels_test.cpp)
target_link_libraries(Kernels_test PRIVATE GTest::gtest_main Kokkos::kokkos)
target_include_directories(Kernels_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
gtest_discover_tests(Kernels_test)
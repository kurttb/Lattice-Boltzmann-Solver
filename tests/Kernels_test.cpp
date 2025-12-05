#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include "EquilibriumKokkos.hpp"
#include "CollisionKokkos.hpp"
#include "StreamingKokkos.hpp"
#include "BCsKokkos.hpp"
#include "GridTypes.hpp"
#include "ComputeStateKokkos.hpp"

// Define vec2_const if not defined (hack for testing if missing in header)
// But better to fix the header if it's missing.
// Let's assume it might be missing and see if it compiles.

/*
 * KernelTest
 * 
 * Motivation:
 * The LBM solver relies on several core kernels (Equilibrium, Collision, Streaming, BCs) 
 * that are executed in parallel on the GPU/CPU via Kokkos. Testing these in isolation 
 * ensures the physics at the microscopic level is correct before assembling the full solver.
 * 
 * Mechanism:
 * - EquilibriumRestTest: Verifies that a fluid at rest (u=0) produces the correct Maxwellian 
 *   weights. This is the ground state of the simulation.
 * - EquilibriumMovingTest: Verifies that a moving fluid produces an anisotropic distribution 
 *   where particles moving in the flow direction have higher population than those opposing it.
 * - CollisionConservationTest: The collision operator must conserve mass (and momentum). 
 *   This test sums the distribution function before and after collision to ensure mass is constant.
 * - StreamingTest: Verifies the advection step. A particle placed at one node with a specific 
 *   velocity direction must appear at the correct neighbor node in the next step.
 * - BounceBackTest: Verifies the no-slip boundary condition. Particles hitting a wall must 
 *   be reflected back in the opposite direction.
 */

TEST(KernelTest, EquilibriumRestTest) {
    int N = 1;
    Kokkos::View<float*> rho("rho", N);
    Kokkos::View<float*> ux("ux", N);
    Kokkos::View<float*> uy("uy", N);
    Kokkos::View<float**> f("f", N, 9);

    auto rho_h = Kokkos::create_mirror_view(rho);
    auto ux_h = Kokkos::create_mirror_view(ux);
    auto uy_h = Kokkos::create_mirror_view(uy);

    rho_h(0) = 1.0f;
    ux_h(0) = 0.0f;
    uy_h(0) = 0.0f;

    Kokkos::deep_copy(rho, rho_h);
    Kokkos::deep_copy(ux, ux_h);
    Kokkos::deep_copy(uy, uy_h);

    Kokkos::parallel_for("EqTest", N, CalcEq(rho, ux, uy, f));
    Kokkos::fence();

    auto f_h = Kokkos::create_mirror_view(f);
    Kokkos::deep_copy(f_h, f);

    // Check weights for rest
    float w[9] = {4.0f/9.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f};
    
    for(int k=0; k<9; ++k) {
        EXPECT_NEAR(f_h(0, k), w[k], 1e-5);
    }
}

TEST(KernelTest, EquilibriumMovingTest) {
    int N = 1;
    Kokkos::View<float*> rho("rho", N);
    Kokkos::View<float*> ux("ux", N);
    Kokkos::View<float*> uy("uy", N);
    Kokkos::View<float**> f("f", N, 9);

    auto rho_h = Kokkos::create_mirror_view(rho);
    auto ux_h = Kokkos::create_mirror_view(ux);
    auto uy_h = Kokkos::create_mirror_view(uy);

    rho_h(0) = 1.0f;
    ux_h(0) = 0.1f; // Moving East
    uy_h(0) = 0.0f;

    Kokkos::deep_copy(rho, rho_h);
    Kokkos::deep_copy(ux, ux_h);
    Kokkos::deep_copy(uy, uy_h);

    Kokkos::parallel_for("EqMovingTest", N, CalcEq(rho, ux, uy, f));
    Kokkos::fence();

    auto f_h = Kokkos::create_mirror_view(f);
    Kokkos::deep_copy(f_h, f);

    // Expect East (1) > West (5)
    EXPECT_GT(f_h(0, 1), f_h(0, 5));
    // Expect NE (2) > SW (6)
    EXPECT_GT(f_h(0, 2), f_h(0, 6));
    // Expect SE (8) > NW (4)
    EXPECT_GT(f_h(0, 8), f_h(0, 4));
    // Expect North (3) == South (7) (Symmetry)
    EXPECT_NEAR(f_h(0, 3), f_h(0, 7), 1e-5);
}

TEST(KernelTest, CollisionConservationTest) {
    int N = 1;
    Kokkos::View<float*> rho("rho", N);
    Kokkos::View<float*> ux("ux", N);
    Kokkos::View<float*> uy("uy", N);
    Kokkos::View<float**> f("f", N, 9);

    auto rho_h = Kokkos::create_mirror_view(rho);
    auto ux_h = Kokkos::create_mirror_view(ux);
    auto uy_h = Kokkos::create_mirror_view(uy);
    auto f_h = Kokkos::create_mirror_view(f);

    rho_h(0) = 1.0f;
    ux_h(0) = 0.1f;
    uy_h(0) = 0.0f;
    
    // Initialize f with some values (e.g. equilibrium)
    float w[9] = {4.0f/9.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f};
    for(int k=0; k<9; ++k) f_h(0, k) = w[k]; // Just approximation

    Kokkos::deep_copy(rho, rho_h);
    Kokkos::deep_copy(ux, ux_h);
    Kokkos::deep_copy(uy, uy_h);
    Kokkos::deep_copy(f, f_h);

    float omega = 1.0f; // Relaxation parameter

    Kokkos::parallel_for("CollTest", N, ComputeCollision(rho, ux, uy, f, omega, 0.0f));
    Kokkos::fence();

    Kokkos::deep_copy(f_h, f);

    // Check mass conservation: sum(f) should be rho
    float sum_f = 0.0f;
    for(int k=0; k<9; ++k) {
        sum_f += f_h(0, k);
    }
    EXPECT_NEAR(sum_f, rho_h(0), 1e-5);
}

TEST(KernelTest, StreamingTest) {
    int Nx = 3;
    int Ny = 3;
    int N = Nx * Ny;
    Kokkos::View<float**> f("f", N, 9);
    Kokkos::View<float**> fstream("fstream", N, 9);

    auto f_h = Kokkos::create_mirror_view(f);
    
    // Set a particle moving East (k=1) at center (1,1) -> index 4
    // Center is i=1, j=1. n = 1 + 1*3 = 4.
    f_h(4, 1) = 1.0f; 
    
    Kokkos::deep_copy(f, f_h);
    Kokkos::deep_copy(fstream, 0.0f);

    Kokkos::parallel_for("StreamTest", N, ComputeStreaming(f, fstream, Nx, Ny));
    Kokkos::fence();

    auto fstream_h = Kokkos::create_mirror_view(fstream);
    Kokkos::deep_copy(fstream_h, fstream);

    // Expected destination: (1+1, 1) = (2, 1) -> index 2 + 1*3 = 5.
    EXPECT_EQ(fstream_h(5, 1), 1.0f);
    
    // Check source is empty (in fstream)
    EXPECT_EQ(fstream_h(4, 1), 0.0f);
}

TEST(KernelTest, BounceBackTest) {
    int Nx = 1;
    int Ny = 1;
    int N = Nx * Ny;
    Kokkos::View<float**> f("f", N, 9);
    
    auto f_h = Kokkos::create_mirror_view(f);
    
    // Simulate Bottom Wall
    // Incident: 8 (SE), 7 (S), 6 (SW)
    // Reflected: 4 (NW), 3 (N), 2 (NE)
    
    int f_inc[3] = {8, 7, 6};
    int f_ref[3] = {4, 3, 2};
    
    // Set incident values
    f_h(0, 8) = 0.8f;
    f_h(0, 7) = 0.7f;
    f_h(0, 6) = 0.6f;
    
    Kokkos::deep_copy(f, f_h);
    
    // Run BounceBack on the single node (0,0)
    Kokkos::parallel_for("BounceBackTest", 
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {1, 1}),
        BounceBack(f, f_inc, f_ref, Nx)
    );
    Kokkos::fence();
    
    Kokkos::deep_copy(f_h, f);
    
    // Check reflections
    EXPECT_EQ(f_h(0, 4), 0.8f); // 8 -> 4
    EXPECT_EQ(f_h(0, 3), 0.7f); // 7 -> 3
    EXPECT_EQ(f_h(0, 2), 0.6f); // 6 -> 2
}

TEST(KernelTest, GridIndexingTest) {
    LBM::CartesianGrid2D grid;
    grid.Nx = 10;
    grid.Ny = 5;
    
    // Test (0,0)
    EXPECT_EQ(grid.idx(0, 0), 0);
    
    // Test (1,0) -> 1
    EXPECT_EQ(grid.idx(1, 0), 1);
    
    // Test (0,1) -> Nx = 10
    EXPECT_EQ(grid.idx(0, 1), 10);
    
    // Test (Nx-1, Ny-1) -> 9 + 4*10 = 49
    EXPECT_EQ(grid.idx(9, 4), 49);
}

TEST(KernelTest, ForceApplicationTest) {
    int N = 1;
    Kokkos::View<float*> rho("rho", N);
    Kokkos::View<float*> ux("ux", N);
    Kokkos::View<float*> uy("uy", N);
    Kokkos::View<float**> f("f", N, 9);

    auto f_h = Kokkos::create_mirror_view(f);
    
    // Initialize f to rest (rho=1, u=0)
    float w[9] = {4.0f/9.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f};
    for(int k=0; k<9; ++k) f_h(0, k) = w[k];
    
    Kokkos::deep_copy(f, f_h);
    
    float Fx = 0.1f;
    float Fy = 0.0f;
    float tau = 1.0f;
    
    // We need to initialize rho, ux, uy because ComputeState writes to them, but it also reads f.
    // Actually ComputeState calculates rho, ux, uy FROM f, and then adds force.
    // So we don't need to init rho, ux, uy.
    
    Kokkos::parallel_for("ForceTest", N, ComputeState(rho, ux, uy, f, Fx, Fy, tau));
    Kokkos::fence();
    
    auto ux_h = Kokkos::create_mirror_view(ux);
    auto rho_h = Kokkos::create_mirror_view(rho);
    Kokkos::deep_copy(ux_h, ux);
    Kokkos::deep_copy(rho_h, rho);
    
    // Expected ux = (0 / 1) + (Fx * tau / 1) = 0.1
    EXPECT_NEAR(ux_h(0), 0.1f, 1e-5);
    EXPECT_NEAR(rho_h(0), 1.0f, 1e-5);
}

TEST(KernelTest, TangentVelocityBCTest) {
    int Nx = 1;
    int Ny = 1;
    int N = Nx * Ny;
    Kokkos::View<float**> f("f", N, 9);
    
    auto f_h = Kokkos::create_mirror_view(f);
    
    // Simulate Top Wall (Moving East)
    // Incident: 2 (NE), 3 (N), 4 (NW) -> These are going INTO the wall (North-ward)
    // Reflected: 6 (SW), 7 (S), 8 (SE) -> These are coming FROM the wall (South-ward)
    
    int f_inc[3] = {2, 3, 4};
    int f_ref[3] = {6, 7, 8};
    float U_wall = 0.1f;
    float cs2 = 1.0f/3.0f;
    
    // Set incident values (going into wall)
    f_h(0, 2) = 0.2f;
    f_h(0, 3) = 0.3f;
    f_h(0, 4) = 0.4f;
    
    // Set other values to 0 for density calc
    f_h(0, 0) = 0.0f;
    f_h(0, 1) = 0.0f;
    f_h(0, 5) = 0.0f;
    f_h(0, 6) = 0.0f;
    f_h(0, 7) = 0.0f;
    f_h(0, 8) = 0.0f;
    
    Kokkos::deep_copy(f, f_h);
    
    Kokkos::parallel_for("TanVelTest", 
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {1, 1}),
        TangentVelocity(f, f_inc, f_ref, U_wall, cs2, Nx)
    );
    Kokkos::fence();
    
    Kokkos::deep_copy(f_h, f);
    
    // Calculate expected rho (sum of all f's we set)
    float rho_expected = 0.2f + 0.3f + 0.4f; 
    
    // Check 6 (SW) from 2 (NE)
    // f_ref = f_inc - 2*w*rho*U/cs2
    // w[6] = 1/36
    float w6 = 1.0f/36.0f;
    float expected_6 = 0.2f - 2.0f * w6 * rho_expected * U_wall / cs2;
    
    EXPECT_NEAR(f_h(0, 6), expected_6, 1e-5);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    Kokkos::initialize(argc, argv);
    int result = RUN_ALL_TESTS();
    Kokkos::finalize();
    return result;
}

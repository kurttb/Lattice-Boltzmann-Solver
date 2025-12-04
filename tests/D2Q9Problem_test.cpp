#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include "D2Q9Problem.hpp"

/*
 * D2Q9ProblemTest
 * 
 * Motivation:
 * These tests verify the high-level functionality of the D2Q9Problem class, which acts as the 
 * main driver for the simulation. It ensures that the class can be instantiated, configured, 
 * and executed without runtime errors.
 * 
 * Mechanism:
 * - InitializationTest: Checks if the constructor runs successfully.
 * - SetParametersTest: Verifies that setting physics parameters (viscosity, ICs, forces) 
 *   does not cause crashes or invalid states.
 * - RunSimulationTest: Performs a minimal integration test by running a single time step 
 *   of the simulation with boundary conditions applied. This ensures the full pipeline 
 *   (Initialize -> Collision -> Streaming -> BCs) is connected correctly.
 */

TEST(D2Q9ProblemTest, InitializationTest) {
    int Nx = 10;
    int Ny = 10;
    LBM::D2Q9Problem prob(Nx, Ny);
    
    // Just verify we can instantiate it
    SUCCEED();
}

TEST(D2Q9ProblemTest, SetParametersTest) {
    int Nx = 10;
    int Ny = 10;
    LBM::D2Q9Problem prob(Nx, Ny);

    prob.setViscosity(0.1f);
    prob.setIC(1.0f, 0.0f, 0.0f);
    prob.setNumTimeSteps(10);
    prob.setForces(0.0f, 0.0f);

    // Verify no crash
    SUCCEED();
}

TEST(D2Q9ProblemTest, RunSimulationTest) {
    int Nx = 10;
    int Ny = 10;
    LBM::D2Q9Problem prob(Nx, Ny);

    prob.setViscosity(0.1f);
    prob.setIC(1.0f, 0.0f, 0.0f);
    prob.setNumTimeSteps(1); // Run 1 step
    
    // Set BCs (required for runSimulation usually)
    prob.setBC("Top", "Periodic");
    prob.setBC("Bottom", "Periodic");
    prob.setBC("Left", "Periodic");
    prob.setBC("Right", "Periodic");

    prob.runSimulation();

    SUCCEED();
}

TEST(D2Q9ProblemTest, ZeroTimeStepsTest) {
    int Nx = 10;
    int Ny = 10;
    LBM::D2Q9Problem prob(Nx, Ny);

    float rho0 = 2.5f;
    float ux0 = 0.1f;
    float uy0 = -0.1f;

    prob.setIC(rho0, ux0, uy0);
    prob.setNumTimeSteps(0);
    
    // Even with 0 steps, we need BCs set to avoid crashes if they are checked during init (though they aren't currently)
    prob.setBC("Top", "Periodic");
    prob.setBC("Bottom", "Periodic");
    prob.setBC("Left", "Periodic");
    prob.setBC("Right", "Periodic");

    prob.runSimulation();

    // Verify that the state matches the initial condition
    std::vector<float> rho = prob.getRho();
    std::vector<float> ux = prob.getUx();
    std::vector<float> uy = prob.getUy();

    ASSERT_EQ(rho.size(), Nx * Ny);
    
    for(size_t i=0; i<rho.size(); ++i) {
        EXPECT_FLOAT_EQ(rho[i], rho0);
        EXPECT_FLOAT_EQ(ux[i], ux0);
        EXPECT_FLOAT_EQ(uy[i], uy0);
    }
}

TEST(D2Q9ProblemTest, InvalidConstructorTest) {
    EXPECT_THROW(LBM::D2Q9Problem(0, 10), std::invalid_argument);
    EXPECT_THROW(LBM::D2Q9Problem(10, 0), std::invalid_argument);
    EXPECT_THROW(LBM::D2Q9Problem(0, 0), std::invalid_argument);
}

TEST(D2Q9ProblemTest, InvalidViscosityTest) {
    LBM::D2Q9Problem prob(10, 10);
    EXPECT_THROW(prob.setViscosity(-0.1f), std::invalid_argument);
    EXPECT_NO_THROW(prob.setViscosity(0.0f)); // 0 should be allowed (inviscid limit)
    EXPECT_NO_THROW(prob.setViscosity(0.1f));
}

TEST(D2Q9ProblemTest, InvalidDensityTest) {
    LBM::D2Q9Problem prob(10, 10);
    EXPECT_THROW(prob.setIC(-1.0f, 0.0f, 0.0f), std::invalid_argument);
    EXPECT_THROW(prob.setIC(0.0f, 0.0f, 0.0f), std::invalid_argument);
    EXPECT_NO_THROW(prob.setIC(1.0f, 0.0f, 0.0f));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    Kokkos::initialize(argc, argv);
    int result = RUN_ALL_TESTS();
    Kokkos::finalize();
    return result;
}

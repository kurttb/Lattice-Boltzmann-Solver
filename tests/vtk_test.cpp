#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <fstream>
#include "VtkWriter.hpp" // Corrected filename case

/*
 * VtkWriterTest
 * 
 * Motivation:
 * The simulation results are useless if they cannot be visualized. This test ensures that 
 * the VTK writer component correctly generates an output file given valid data.
 * 
 * Mechanism:
 * - WriteOutputTest: Creates dummy density and velocity fields, calls the WriteVtk function, 
 *   and then checks the filesystem to confirm that the file was created and is readable. 
 *   It cleans up the file afterwards to keep the test environment clean.
 */

TEST(VtkWriterTest, WriteOutputTest) {
    int Nx = 10;
    int Ny = 5;
    int N = Nx * Ny;
    std::vector<float> rho(N, 1.0f);
    std::vector<float> ux(N, 5.0f);
    std::vector<float> uy(N, 0.0f);

    std::string file_path = "test_output.vtk";
    std::string pv_title = "LBM Field";

    // Call the function
    LBM::WriteVtk(rho, ux, uy, Nx, Ny, file_path, pv_title);

    // Check if file exists
    std::ifstream f(file_path.c_str());
    EXPECT_TRUE(f.good());
    
    // Clean up
    if (f.good()) {
        f.close();
        std::remove(file_path.c_str());
    }
}

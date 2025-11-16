#include <vector> 
#include <fstream>
#include "vtk_writer.hpp"


int main() {

	int Nx = 10;
	int Ny = 5;
	int N = Nx*Ny;
	std::vector<double> rho(N, 1.0); 
	std::vector<double> ux(N, 5.0); 
	std::vector<double> uy(N, 0.0); 

	std::string file_path = "paraview.vtk";
	std::string pv_title = "LBM Field";

	write_vtk(rho, ux, uy, Nx, Ny, file_path, pv_title);

	return 0;
}

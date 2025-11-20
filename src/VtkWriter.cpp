#include <vector>
#include <string>
#include <fstream>
#include "VtkWriter.hpp"

namespace LBM {

	void WriteVtk(const std::vector<double>& rho, 
				   const std::vector<double>& ux, 
				   const std::vector<double>& uy, 
				   const size_t Nx, const size_t Ny, 
				   const std::string& file_path, 
				   const std::string& pv_title) {


		// Create file 
		std::ofstream file (file_path);

		// Write file type, format and grid type
		file << "# vtk DataFile Version 3.0" << std::endl;
		file << pv_title << std::endl;
		file << "ASCII" << std::endl;
		file << "DATASET STRUCTURED_POINTS" << std::endl;
		file << "\n";

		// Write dimensions, origin, and grid spacing
		file << "DIMENSIONS" << " " << Nx << " " << Ny << " " << 1 << std::endl;
		file << "ORIGIN" << " " << 0 << " " << 0 << " " << 0 << std::endl;

		// Write grid spacing
		file << "SPACING" << " " << 1 << " " << 1 << " " << 1 << std::endl;
		file << "\n";

		// Determine the number of points and write 
		size_t N = Nx*Ny;
		file << "POINT_DATA" << " " << N << std::endl;
		file << "\n";

		// Header for density section
		file << "SCALARS" << " " << "rho" << " " << "double" << " " << 1 << std::endl;
		file << "LOOKUP_TABLE" << " " << "default" << std::endl;

		// Write the density
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				file << rho[i] << " ";
			}
			file << "\n";
		}
		file << "\n";


		// Write Velocity
		file << "VECTORS" << " " << "velocity" << " " << "double" << std::endl;
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				file << ux[i + Nx*j] << " " << uy[i + Nx*j] << " " << 0.0 << std::endl;
			}
			file << "\n";
		}


		// Close file
		file.close();
	}
}

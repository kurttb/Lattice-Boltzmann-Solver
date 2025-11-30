#ifndef VTK_WRITER
#define VTK_WRITER

#include <vector>
#include <string>

namespace LBM {
	void WriteVtk(const std::vector<float>& rho, 
				   const std::vector<float>& ux, 
				   const std::vector<float>& uy, 
				   const size_t Nx, const size_t Ny, 
				   const std::string& file_path, 
				   const std::string& pv_title);

}

#endif

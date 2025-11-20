cmake_minimum_required(VERSION 3.23)

project(LBM VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Add Library
add_library(LBM SHARED)


# Add source files
target_sources(LBM
	PRIVATE
		src/BoundaryConditions.cpp
		src/Collisions.cpp
		src/ComputeState.cpp
		src/D2Q9Problem.cpp
		src/Streaming.cpp
		src/VtkWriter.cpp

	PUBLIC
		FILE_SET HEADERS
		BASE_DIRS
			include
		FILES
			include/BoundaryConditions.hpp
			include/Collisions.hpp
			include/ComputeState.hpp
			include/D2Q9Problem.hpp
			include/GridTypes.hpp
			include/Streaming.hpp
			include/VtkWriter.hpp
)

# Add compiler flag macro
macro(add_compiler_flag FLAG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
endmacro()

# Add compiler flags
add_compiler_flag("-O3")
add_compiler_flag("-ffast-math")
add_compiler_flag("-march=native")
add_compiler_flag("-fopenmp")


# Write Directory
set_target_properties(LBM 
	PROPERTIES
		LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
)

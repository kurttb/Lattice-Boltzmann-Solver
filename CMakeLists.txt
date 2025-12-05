cmake_minimum_required(VERSION 3.15)

project(LBM VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Default build type to release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
#message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
#message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")


# Determine whether build will be static or shared
option(BUILD_SHARED "Build shared library" OFF) # Static by default
if(BUILD_SHARED)
	set(LBM_LIB_TYPE SHARED)
else()
	set(LBM_LIB_TYPE STATIC)
endif()
message(STATUS "BUILD TYPE: ${LBM_LIB_TYPE}") # Build type status



# Define Variables for Kokkos Backend
set(Kokkos_ROOT "/home/kurttb/kokkos-cuda-install" CACHE PATH "Root directory of Kokkos installation")

# Find Kokkos
find_package(Kokkos REQUIRED)


# Set Location for Libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)


# Add Library and its source files
add_library(LBM ${LBM_LIB_TYPE}
		src/D2Q9Problem.cpp
		src/VtkWriter.cpp
)

# Link LBM library to kokkos
target_link_libraries(LBM
	PRIVATE
		Kokkos::kokkos
)

# Add include directory
target_include_directories(LBM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Set location for exectuables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)


# Add exectuables
file(GLOB files "examples/*_case.cpp")
foreach(file ${files})
	
	get_filename_component(name ${file} NAME_WE) # Get name without extension
	add_executable(${name}) # Make an executable

	# Add source file for target
	target_sources(${name}
		PRIVATE
			${file}
	)

	# Link against the LBM library
	target_link_libraries(${name}
		PRIVATE
			LBM
	)

	# -I flag for the specific target
	target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

endforeach()

# Add an install target
include(GNUInstallDirs)

install(TARGETS LBM
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Add an uninstall function




# Add compiler flag macro
macro(add_compiler_flag FLAG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
endmacro()


# Add compiler flags
add_compiler_flag("-Wall")
add_compiler_flag("-Wextra")
add_compiler_flag("-Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops")

enable_testing()
add_subdirectory(tests)


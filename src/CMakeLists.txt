#------------------------------------------------------------------#
# CMakeLists for PELLET code build.
# Shouldn't need to make any changes here!
#------------------------------------------------------------------#

# Set minimum required cmake version
cmake_minimum_required (VERSION 3.5)

# Project name
project (pellet Fortran)

# Disable in source builds 
set(CMAKE_DISABLE_IN_SOURCE_BUILD TRUE)

# Check if trying to build in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are disabled. Please create a build directory and use cmake .. OR run cmake . -B <build-dir> from the source directory.")
endif()

# Create a separate subdirectory for the module files
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# CONFIG_DIR is for the location of user.cmake 
# (or default.cmake) and saves having to specify it explicitly
# when loading.
set(CONFIG_DIR "${CMAKE_SOURCE_DIR}/config")

# Check whether a user.cmake file exists. If yes, use
# that to build. If no, use defaults.cmake.
# The user(defaults).cmake path needs to be in reference to wherever
# the CMakeLists.txt file is. 
if(EXISTS "${CONFIG_DIR}/user.cmake")
    message(STATUS "Loading user-defined settings from user.cmake")
    include("${CONFIG_DIR}/user.cmake")
endif()

# Print build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Print NetCDF header and library paths 
message(STATUS "NetCDF include path is ${NETCDF_INC_PATH}")
message(STATUS "NetCDF library path is ${NETCDF_LIB_PATH}")

# Make sure to link the library and set the include path
link_directories(${NETCDF_LIB_PATH})
include_directories(${NETCDF_INC_PATH} ${CMAKE_Fortran_MODULE_DIRECTORY})

# Searching for the required libraries
find_library(NETCDF_F_LIB NAMES netcdff ${NETCDF_F_LIB_NAME} PATHS ${NETCDF_LIB_PATH})
find_library(NETCDF_C_LIB NAMES netcdf ${NETCDF_C_LIB_NAME} PATHS ${NETCDF_LIB_PATH})

# Setting the library path for both the c++ and Fortran versions 
# of NetCDF 
set(NETCDF_LIB ${NETCDF_F_LIB} ${NETCDF_C_LIB})

# Print whether the NetCDF libraries were found.
if(NETCDF_F_LIB AND NETCDF_C_LIB)
    message(STATUS "NetCDF Fortran library found: ${NETCDF_F_LIB}")
    message(STATUS "NetCDF C library found: ${NETCDF_C_LIB}")
else()
    message(FATAL_ERROR "Failed to find required NetCDF libraries.")
endif()

# Print for successful library linking
message(STATUS "Linking against NetCDF libraries: ${NETCDF_LIB}")

# Set module filenames
set (pellet_sources
    pellet_mod.f90
    ajax_mod.f90
    fluxav_mod.f90
    fluxav_vmec.f
    linear1_mod.f90
    mirth_mod.f90
    pellet_dr.f90
    prl_mod_intel.f90
    setup_ajax.f90
    spec_kind_mod.f90
    spline1_mod.f90
    track_mod.f90
    write_mod.f90
    x_mod.f90
    drifts_mod.f90
)

# Extra files not used immediately
#    fluxav_vmec.f
#    prldata_mod.f90

# Create the executable -- now you should be good to go!
add_executable (xpellet ${pellet_sources})

# add_custom_target(clean-all
#     COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_SOURCE_DIR}/build
#     COMMENT "Removing build directory"
# )

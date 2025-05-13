cmake_minimum_required (VERSION 2.8.12)

project (pellet Fortran)

set(CONFIG_DIR "${CMAKE_SOURCE_DIR}/config")

# The user.cmake path needs to be in reference to wherever
# the CMakeLists.txt file is. 
if(EXISTS "${CONFIG_DIR}/user.cmake")
    message(STATUS "Loading user-defined settings from user.cmake")
    include("${CONFIG_DIR}/user.cmake")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

message(STATUS "NetCDF include path is ${NETCDF_INC_PATH}")
message(STATUS "NetCDF library path is ${NETCDF_LIB_PATH}")

link_directories(${NETCDF_LIB_PATH})
include_directories(${NETCDF_INC_PATH})

find_library(NETCDF_F_LIB NAMES netcdff ${NETCDF_F_LIB_NAME} PATHS ${NETCDF_LIB_PATH})
find_library(NETCDF_C_LIB NAMES netcdf ${NETCDF_C_LIB_NAME} PATHS ${NETCDF_LIB_PATH})

set(NETCDF_LIB ${NETCDF_F_LIB} ${NETCDF_C_LIB})

if(NETCDF_F_LIB AND NETCDF_C_LIB)
    message(STATUS "NetCDF Fortran library found: ${NETCDF_F_LIB}")
    message(STATUS "NetCDF C library found: ${NETCDF_C_LIB}")
else()
    message(FATAL_ERROR "Failed to find required NetCDF libraries.")
endif()

message(STATUS "Linking against NetCDF libraries: ${NETCDF_LIB}")

set (pellet_sources
    pellet_mod.f90
    ajax_mod.f90
    fluxav_mod.f90
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
)

# Extra files not used immediately
#    fluxav_vmec.f
#    prldata_mod.f90

add_executable (xpellet ${pellet_sources})

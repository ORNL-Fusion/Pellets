cmake_minimum_required (VERSION 2.8.12)

project (pellet Fortran)

find_library (NETCDF_F_LIB NAMES netcdff ${NETCDF_F_LIB_NAME})
find_library (NETCDF_C_LIB NAMES netcdf ${NETCDF_C_LIB_NAME})

set (NETCDF_LIB ${NETCDF_F_LIB} ${NETCDF_C_LIB})

include_directories(/usr/local/include)
#include_directories(${NETCDF_INC_PATH})

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

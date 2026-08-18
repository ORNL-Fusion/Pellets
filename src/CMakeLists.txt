set(pellet_sources
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
    drifts_mod.f90
    eqdsk_mod.f90
)

# Extra files not currently used
#    fluxav_vmec.f
#    prldata_mod.f90

add_executable(xpellet)

target_sources(xpellet PRIVATE
    ${pellet_sources}
)

set_target_properties(xpellet PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
        "${PROJECT_BINARY_DIR}"

    Fortran_MODULE_DIRECTORY
        "${PROJECT_BINARY_DIR}/modules"
)

target_link_libraries(xpellet PRIVATE netCDF::netcdff)
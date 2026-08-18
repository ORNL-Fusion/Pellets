#------------------------------------------------------------------#
# CMakeLists for PELLET code build.
#------------------------------------------------------------------#

# Set minimum required cmake version
cmake_minimum_required (VERSION 3.19)

# Project name
project (pellet Fortran)

# Disable in source builds 
set(CMAKE_DISABLE_IN_SOURCE_BUILD TRUE)

#-------------------------------------------------------------------------------
#  Setup build types
#-------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
else()
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()
set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
              Debug
              Release
              MinSizeRel
              RelWithDebInfo
)

find_package(netCDF CONFIG REQUIRED)
find_package(netCDF-Fortran CONFIG REQUIRED)

add_subdirectory(src)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(xpellet PRIVATE
        $<$<CONFIG:Debug>:-Og>
        $<$<CONFIG:Debug>:-Wall>
        $<$<CONFIG:Debug>:-Wextra>
        $<$<CONFIG:Debug>:-fcheck=all>
        $<$<CONFIG:Debug>:-fbacktrace>
        $<$<CONFIG:Debug>:-ffpe-trap=invalid,zero,overflow>

        $<$<CONFIG:RelWithDebInfo>:-O2>
        $<$<CONFIG:RelWithDebInfo>:-Wall>
        $<$<CONFIG:RelWithDebInfo>:-Wextra>

        $<$<CONFIG:Release>:-O2>
    )
endif()
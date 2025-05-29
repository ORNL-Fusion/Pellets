#------------------------------------------------------------------#
# CMake defaults to build the PELLET code.
#
# Pretty minimal given that the code only uses NetCDF :).
#
# To add your personal touch, copy this file to user.cmake and  
# make your desired modifications (e.g., NetCDF file path...). 
#------------------------------------------------------------------#

# Set build type
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type: Debug, Release, etc...")

# Check build type and set flags accordingly 
if(CMAKE_BUILD_TYPE MATCHES "[Dd]ebug")
  set(FORTRAN_FLAGS "-Og -Wall -Wextra -fcheck=all" CACHE STRING "Fortran compiler flags")
else()
  set(FORTRAN_FLAGS "-O2" CACHE STRING "Fortran compiler flags")
endif()

# Set path for the NETCDF headers and libs.
set(NETCDF_INC_PATH "/usr/local/include" CACHE PATH "Path to NetCDF include directory")
set(NETCDF_LIB_PATH "/usr/local/lib" CACHE PATH "Path to NetCDF library directory")

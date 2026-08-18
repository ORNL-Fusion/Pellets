# PELLET
## ORNL implementation of the PELLET model

PELLET uses 0D ablation models to calculate the pellet ablation rate, and can track the pellet path up to 3D. 

## STRUCTURE

```
Pellets/
|--- CMakeLists.txt
|--- CMakePresets.json
|--- CMakeUserPresets.json.example
|--- documentation/
|--- src/
|    |--- CMakeLists.txt
|    |___ *.f90
```

## DEPENDENCIES

PELLET requires gfortran and NetCDF. 

### macOS

Using a package manager like [Homebrew](https://brew.sh) is usually the easiest way. Once installed, run

```
brew install gfortran netcdf
```

For Macs with an Intel chip, brew's preferred install location is `/usr/local`. For Apple Silicon macs, it's `/opt/homebrew`.

### Debian based *nix flavours

```
sudo apt install gfortran netcdf
```

This method will likely require root access. If you don't have it, you might have to download and build from source.

## BUILD

To build using CMake, copy `CMakeUserPresets.json.example` to `CMakeUserPresets.json`, update the `CMAKE_PREFIX_PATH` and run from the top directory,

```
cmake --preset build
cmake --build build
```

which puts executable `xpellet` in the build directory.

> [!NOTE]
> The presets file is just a simple example...
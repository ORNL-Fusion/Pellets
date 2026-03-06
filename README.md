# PELLET
## ORNL implementation of the PELLET model

PELLET uses 0D ablation models to calculate the pellet ablation rate, and can track the pellet path up to 3D. 


> [!NOTE]
> The `rlb-dev` branch will no longer be updated. Active development of the PELLET code to include drift terms (the purpose of the original branch) has moved to the branch `feature/include-gradb-drifts`

## STRUCTURE

```
Pellets/
|--- src/
|--- Python/
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

To build using CMake, copy the `defaults.cmake` to a new `user.cmake` file and make any changes you need to e.g., the library paths. Once done, run `cmake . -B build` in the `src/` directory.

> [!NOTE]
> In source builds are not supported, so you must create a "build" directory.

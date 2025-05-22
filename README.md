# PELLET
## ORNL implementation of the PELLET model

PELLET uses 0D ablation models to calculate the pellet ablation rate, and can track the pellet path up to 3D. 

## Structure

```
Pellets/
|--- src/
|	|---
```

## BUILD

PELLET requires NetCDF. 

To build using CMake, copy the `defaults.cmake` to a new `user.cmake` file and make any changes you need to e.g., the library paths. Once done, run `cmake . ` in the `src/` directory.
# Overkit

Overset meshes are one of several approaches to representing complex geometry in CFD and other
  types of physical simulations.
An overset mesh consists of a set of overlapping component meshes (often structured grids) that
  each represent a portion of the simulation domain.
The component meshes communicate with other parts of the domain by interpolating from neighboring
  grids at their outer edges.
This is often simpler than trying to create a single unstructured mesh to represent the
  entire domain, and when used with structured grids also enables the use of efficient and
  accurate numerical methods.

Overkit automates tasks commonly required for setting up overset meshes, such as:

1. Hole cutting of boundary geometry and overlapped coarse grids
2. Determination of donors and receivers for inter-grid communication
3. Generation of interpolation weights

Disclaimer: Overkit is still in experimental stages of development, so it's likely that some things 
  won't work as well as they could.
It's also likely that the API will change from version to version.
You have been warned. :)

# Installation

## Requirements

Overkit requires working Fortran and C compilers (tested on GNU >= 4.4, Intel >= 12) and
  CMake.
It currently has no other third-party library/tool dependencies.

## Building

Overkit uses CMake, so the build process is similar to that of other CMake-based build systems.
Here's an example build procedure:
```
  cd ~/overkit-src
  mkdir build && cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/overkit
  make
  make install
```

See CMake's documentation for more details.

Use the flag `-DEXAMPLES=On` when configuring to build some example programs that use Overkit to
  generate overset grids.

# Documentation

Documentation will remain somewhat sparse for the time being as I'm planning to do some
  restructuring of the API.
In the meantime, the best way to learn how it works is to look at the code in the `examples`
  directory and browse the source code corresponding to the routines that are used.

# Issues

There are currently several known issues:

* The disjoint fringes option isn't fully implemented for cubic interpolation. Use fringe padding
  to compensate -- a value of 2 or 3 should be enough.
* Boundary hole cutting functionality is rather fragile at the moment. If it fails, try turning it
  off and roughly pre-blanking out the points that overlap with the boundary region (blanking a
  square region enclosing the boundary is often enough).
* Donor selection needs to be improved; it occasionally gives unexpected results near junctions of
  three or more grids.

# Acknowledgements

Overkit is one of several software projects produced by the University of Illinois Center for
  Exascale Simulation of Plasma-Coupled Combustion (XPACC): http://xpacc.illinois.edu

This material is based in part upon work supported by the Department of Energy, National Nuclear
  Security Administration, under Award Number DE-NA0002374.

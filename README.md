_(Looking for the newer parallel C++ version? It's [here](https://github.com/majosm/overkit).)_

# Overkit

Overset meshes are a method for representing complex geometry in computational fluid dynamics
and other types of PDE-based simulations. An overset mesh is a composite mesh made up of a set
of overlapping component grids (structured or unstructured), in which the individual grids
exchange data with each other via interpolation.

![An overset mesh](misc/readme/blobs.png)

Overset meshes have several nice properties, including:

1. **Simplicity:** Component grids are (often) structured, so operating on them is similar to
   operating on Cartesian grids
2. **Flexibility:** Overset meshes' modular nature make them convenient for moving geometry
   problems and adaptive mesh refinement, as well as simulations in which different physical
   phenomena or numerical methods are relevant in different regions
3. **Efficiency:** Structured grids enable the use of efficient and accurate numerical methods.

However, some processing is required to convert a set of overlapping grids into a working overset
mesh. Overkit is designed to automate the tasks involved in this processing, such as:

1. Hole cutting of boundary geometry and overlapped coarse grids
2. Determination of donors and receivers for inter-grid communication
3. Generation of interpolation weights.

Overkit is also capable of generating interpolation data for mesh remapping (more on this below).

_Note:_ Overkit is currently designed for structured grids with node-centered data. Support for
cell-centered data and unstructured grids may be added in future versions.

# Installation

## Requirements

Overkit requires Fortran/C compilers and CMake. It currently has no other third-party library/tool
dependencies.

## Building

Overkit uses CMake, so the build process is similar to that of other CMake-based build systems.
An example build procedure is:

```bash
  cd ~/overkit-src
  mkdir build && cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/overkit
  make
  make install
```

See CMake's documentation for more details.

### Examples

Use the flag **`-DEXAMPLES=ON`** when configuring to build some example programs that use
Overkit. Source code for the examples can be found in `examples/`; the built programs are placed
in `<cmake-build-dir>/examples/`. Run `<example-program-name> --help` for information on how to
use each one.

### OpenMP

Overkit has optional OpenMP support that can be enabled by configuring with **`-DOPENMP=ON`**.

# Usage

What follows is a brief(ish) Overkit user's guide.

_Note:_ It may be helpful to look at some of the examples' sources while reading this.

## Overset basics and terminology

In an overset mesh, each component grid has **_fringes_** along its non-boundary edges. A
fringe is made up of one or more layers of **_receiver_** points, which are points that receive
data from other grids. Each receiver has a corresponding **_donor_** on the sending grid which
consists of the set of points that makes up the interpolation stencil and each point's
associated interpolation weight. Overkit calls the combination of a donor and a receiver a
**_connection_**. Most of the time, a pair of grids in an overset mesh will communicate in both
directions, so an interface between them will have two adjacent fringes (one on each grid).

In some cases when assembling an overset mesh, suitable donors may not be found for certain
receiver points (e.g. if there is insufficient overlap between grids); these receiver points
are referred to as **_orphan_** points.

Overset grids often have regions in which solutions are not computed called **_holes_**. These
result from several processes during assembly (details in the assembly section below).

In a typical explicit solver the solution procedure for an overset mesh may be something like:

```
Set u to some initial value for all grid points
For istep in 1..nsteps
  Solve for u at non-hole, non-fringe points
  Interpolate to update u at fringe points
```

_Note:_ When designing overset meshes for solvers that follow the above procedure, it is important 
for the grids to overlap enough to allow donors to contain only non-fringe points (otherwise data
will be used in the interpolation stencils that has not yet been updated to the current step).

## Data structures

Overkit has three main data structures: domains, grids, and connectivities.

### Domain (`ovk_domain`)

The domain is the primary Overkit data structure. It contains a set of grids and their associated
connectivities.

### Grid (`ovk_grid`)

The grid represents a single component grid. It holds the grid's metadata (size, periodicity, etc.),
its coordinates, and a field representing the state of each grid point (interior point, boundary
point, hole, etc.).

### Connectivity (`ovk_connectivity`)

The connectivity holds inter-grid connectivity data for a given pair of grids.
For each connection (donor/receiver pair), the following data is stored:

1. _Donor extents:_ the range of points on the donor grid from which data will be interpolated,
   stored as a lower and upper point
2. _Donor interpolation coefficients:_ the interpolation weights for each point in the interpolation
   stencil, stored in tensor product format (i.e., an array of one-dimensional stencils, with the
   first index indicating the point in the stencil and the second index indicating the dimension)
3. _Donor coordinates:_ the mapped coordinates of the receiver point inside the donor stencil (the
   coordinate range depends on the interpolation scheme; e.g., for linear interpolation the
   coordinates range from 0 to 1, with (0,0) representing the lower left vertex and (1,1)
   representing the upper right vertex). _Note:_ These are typically only needed when generating
   interpolation coefficients separately from Overkit
4. _Receiver point:_ the point on the receiver grid to which the interpolated data will be sent.

### Auxiliary types

**Cart (`ovk_cart`):** Type representing a grid's structure, namely its size and periodicity.

**Field (`ovk_field_<int,large_int,real,logical>`):** Type for grid data. Data is stored in the
`values` member.

**Array (`ovk_array_<int,large_int,real,logical>`):** Type for linear arrays. Data is stored in the
`values` member.
  
## Domain creation

Domains are created with the following command:

```fortran
  call ovkCreateDomain(Domain, <numdims>, <numgrids>, StatusLogFile=<statusunit>, &
    ErrorLogFile=<errorunit>)
```

`<numdims>` is the dimension of the problem (2 or 3).

`<numgrids>` is the number of grids in the domain.

`<statusunit>` and `<errorunit>` are files to which to send the status and error output. They can
be set to `OUTPUT_UNIT` and `ERROR_UNIT` from `iso_fortran_env` to output to the command line, or
the arguments can be omitted to suppress output.

_Note:_ There is also a command to clean up once finished:

```fortran
  call ovkDestroyDomain(Domain)
```

## Grid creation

Grids are created with the following command:

```fortran
  call ovkCreateGrid(Domain, <gridid>, <numpoints>, Periodic=<periodic>, &
    PeriodicStorage=<periodicstorage>, PeriodicLength=<periodiclength>, &
    GeometryType=<geometrytype>)
```

`<gridid>` is the index of the grid in the domain.

`<numpoints>` is the number of points in each direction.

`<periodic>` specifies whether the grid is periodic or not (`.true.` or `.false.`) in each direction.
_(optional; default=`.false.`)_.

`<periodicstorage>` can either be **`OVK_PERIODIC_STORAGE_UNIQUE`** or
**`OVK_PERIODIC_STORAGE_DUPLICATED`**. The former indicates that the grid contains only one instance
of each point along any of its periodic dimensions; the latter indicates that the first point is 
duplicated at the end in any periodic dimension. _(optional; default=`OVK_PERIODIC_STORAGE_UNIQUE`)_.
  
`<periodiclength>` specifies the offset that must be applied to the grid's coordinates
when crossing over from one period to the next. The length can be set to 0 for O-shaped grids.
_Note:_ Overkit currently assumes that only the corresponding dimension's coordinate changes when
crossing periods, e.g., if crossing from point (i,N,k) to (i,1,k), the X and Z coordinates will
remain the same, and the Y coordinate will be offset by PeriodicLength(2). _(optional; default=0.0)_.

`<geometrytype>` can be used to provide hints to Overkit about whether the grid's coordinates have
any special properties that can be used to improve performance. Possible values are:

* **`OVK_GEOMETRY_UNIFORM`** - an axis-aligned uniformly-spaced grid
* **`OVK_GEOMETRY_ORIENTED_UNIFORM`** - a rotated uniform grid
* **`OVK_GEOMETRY_RECTILINEAR`** - a non-uniform grid in which all of the cells are
  rectangular and axis-aligned
* **`OVK_GEOMETRY_ORIENTED_RECTILINEAR`** - a rotated rectilinear grid
* **`OVK_GEOMETRY_CURVILINEAR`** - a general non-uniform structured grid.

_(optional; default=`OVK_GEOMETRY_CURVILINEAR`)_.

_Note:_ There is also a command to destroy a grid:

```fortran
  call ovkDestroyGrid(Domain, <gridid>)
```

(but most of the time this isn't necessary, as the domain will clean up all of its contained grids
and connectivities when `ovkDestroyDomain` is called).

### Setting coordinates

Grids' coordinates are set using the following sequence of commands:

```fortran
  ! Request to edit the grid (returns pointer to ovk_grid)
  call ovkEditGrid(Domain, <gridid>, Grid)

  ! Request to edit the grid's coordinates (returns pointers to ovk_field_real)
  ! 2D is similar but with Z omitted
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  call ovkEditGridCoords(Grid, 3, Z)
  
  ! Set the coordinates
  ! 2D is similar but with Z omitted and k equal to 1
  do k = 1, Nz
    do j = 1, Ny
      do i = 1, Nx
        X%values(i,j,k) = <xvalue>
        Y%values(i,j,k) = <yvalue>
        Z%values(i,j,k) = <zvalue>
      end do
    end do
  end do
  
  ! Tell the grid that you're done editing the coordinates
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)
  call ovkReleaseGridCoords(Grid, Z)

  ! Tell the domain that you're done editing the grid
  call ovkReleaseGrid(Domain, Grid)
```

Grid coordinates are initially set to uniformly vary from 0 to N-1 in each direction (which is
sometimes convenient as a starting point for generating grids).

### Setting state

Grids' states are set with the following sequence of commands:

```fortran
  ! Request to edit the grid (returns pointer to ovk_grid)
  call ovkEditGrid(Domain, <gridid>, Grid)

  ! Request to edit the grid's state (returns pointer to ovk_field_int)
  call ovkEditGridState(Grid, State)
  
  ! Set the state
  ! 2D is similar but with k equal to 1
  do k = 1, Nz
    do j = 1, Ny
      do i = 1, Nx
        State%values(i,j,k) = <statevalue>
      end do
    end do
  end do
  
  ! Tell the grid that you're done editing the state
  call ovkReleaseGridState(Grid, State)

  ! Tell the domain that you're done editing the grid
  call ovkReleaseGrid(Domain, Grid)
```
  
`<statevalue>` can be one of the following:

* **`OVK_INTERIOR_POINT`** - (the default) a normal point
* **`OVK_DOMAIN_BOUNDARY_POINT`** - a point that lies along a boundary of the simulation domain
* **`OVK_INTERNAL_BOUNDARY_POINT`** - a point along a grid edge that isn't a domain boundary, but
  where a fringe is undesirable (perhaps due to the presence of some other type of inter-grid
  interface)
* **`OVK_EXTERIOR_POINT`** - a point that is outside of the simulation domain.

_Note 1:_ Coordinates and state can be edited at the same time (in fact it's slightly more efficient
to do this).

_Note 2:_ Overkit expects valid grid coordinates even for points set to `OVK_EXTERIOR_POINT`.

_Note 3:_ The above state values are actually combinations of individual state flags of the form
`OVK_STATE_*` (more on these later).

## Assembly

Once the domain and grids are set up, assembly can be run with the following commands:

```fortran
  AssemblyOptions = ovk_assembly_options_(<numdims>, <numgrids>)
  call ovkAssemble(Domain, AssemblyOptions)
```

However, for anything to actually happen we will need to first modify the assembly options. The
sections below describe the assembly process in Overkit and explain what each of the options in the
**`ovk_assembly_options`** structure do. Commands for modifying assembly options will have one of
the following forms:

```fortran
  ! Per-grid settings
  call ovkSetAssemblyOption<optionname>(AssemblyOptions, <gridid>, <optionvalue>)
  
  ! Per-grid-pair settings
  call ovkSetAssemblyOption<optionname>(AssemblyOptions, <gridid1>, <gridid2>, <optionvalue>)
```

A value of **`OVK_ALL_GRIDS`** can be passed to any of the `<gridid*>` parameters to modify the
option value for multiple grids at once.

### Overlap detection

Overkit first needs to determine how the component grids are situated with respect to each other.
This consists of looking at every grid point and determining within which cell(s) on the other grids
the point resides (if any). Naively, overlap detection has quadratic computational complexity in the
number of grid points (because it involves testing each point against each cell among all pairs of
grids), but this can be improved to linear complexity through the use of spatial partitioning data
structures.

Overkit currently uses a hybrid approach consisting of a kd-tree at coarse levels (to prune out
empty space and to separate regions of vastly different grid resolution) and a spatial hash grid at
finer levels.
  
Relevant options:

* **`Overlappable(m,n)`** - whether grid `m` is considered to potentially overlap with grid `n`
   _(default=`.false.`)_
* **`OverlapTolerance(m,n)`** - tolerance for points on grid `n` that lie slightly outside of grid
  `m` but should be considered to overlap (most often occurs along curved boundaries). Specified as
  a fraction of the boundary cells' sizes (e.g. a value of 0.1 means that a point can be up to 10%
  outside of a nearby cell and still be considered to overlap) _(default=0.0)_
* **`OverlapAccelDepthAdjust(m)`** - adjusts the desired minimum number of cells in each kd-tree
  node on grid `m`. Each increment by 1.0 halves the number (increasing the tree depth), and each
  decrement by 1.0 doubles the number (decreasing the tree depth) _(default=0.0)_
* **`OverlapAccelResolutionAdjust(m)`** - adjusts the size of hash grid bins on grid `m`. Each
  increment by 1.0 halves the bin size (resulting in fewer cells overlapping with each bin,
  i.e., better performance at the cost of more memory) and each decrement by 1.0 doubles the
  bin size _(default=0.0)_.

### Boundary hole cutting

The collection of domain boundary points on all of the grids defines the overall shape of the
domain. Depending on the input grid configuration, some points may lie outside of this region and
will need to be removed. (For example, in a two-grid airfoil configuration with a uniform background
grid and a second grid that wraps around the airfoil geometry, points on the background grid will
lie inside the non-domain region defined by inner edge of the second grid.)

Overkit handles this by (1) projecting a representation of other grids' boundaries onto the grid
being cut, (2) detecting which side is outside of the domain, then (3) flooding this side and
removing the detected points from the domain.

Relevant options:

* **`InferBoundaries(n)`** - assume grid edges in non-overlapping regions are domain boundaries on
  grid `n` _(default=.false.)_
* **`CutBoundaryHoles(m,n)`** - whether the boundary points on grid `m` are allowed to be used in
  cutting grid `n` _(default=.false.)_.

### Receiver selection

Next we need to determine which points should be receiving data. Overkit divides the points detected
in this step into two categories, _**outer fringe**_ points and _**occluded**_ points. Outer fringe
points are points near non-boundary grid edges. Occluded points are overlapped points for which one
of the overlapping grids has been determined to be a better representation of that part of the
domain than the present one.

Relevant options:

* **`FringeSize(n)`** - how many layers of points away from non-boundary edges of grid `n` should be
  counted as fringe _(default=0)_
* **`Occludes(m,n)`** - how to treat points on grid `n` that are overlapped by grid `m`. Possible
  values are **`OVK_OCCLUDES_ALL`**, meaning all such points are occluded, **`OVK_OCCLUDES_NONE`**,
  meaning none are occluded, and **`OVK_OCCLUDES_COARSE`**, meaning only points with a lower
  resolution are occluded _(default=`OVK_OCCLUDES_NONE`)_
* **`EdgePadding(m,n)`** - in some cases it may be desirable to prevent points from receiving data
  from points that are within some distance of the donor grid's edge. This option specifies how
  far from the edge cells should be on the overlapping grid `m` in order to occlude points on the
  overlapped grid `n` _(default=0)_
* **`EdgeSmoothing(n)`** - used in conjunction with `EdgePadding`, helps to remove any sharp edges
  or disconnectedness in the set of occluded points on grid `n`. Higher values increase smoothing
  strength _(default=0)_.

### Overlap minimization

For overset meshes we typically don't need to keep all of the receiver points detected in the
previous step. Instead we can cut out most of the occluded points and leave only the few layers
nearest to the unoccluded part of the grid. Overkit calls these remaining points the _**inner
fringe**_.

The `EdgePadding` option described in the previous section reduces the size of a grid's
occluded region, which has the side effect of increasing the amount of overlap remaining after
overlap minimization. This can be useful in certain situations such as that of the explicit solver
described above.

Relevant options:

* **`MinimizeOverlap(m,n)`** - whether points on grid `n` that are occluded by grid `m` and aren't
  part of the inner fringe should be removed _(default=`.false.`)_
* **`FringeSize(n)`** - (defined above)
* **`EdgePadding(m,n)`** - (defined above).

### Donor selection

By this point we have our final set of receiver points and the only thing left to do is to locate
the corresponding donors. This involves choosing which overlapping grid will donate to each receiver
and then actually generating the interpolation stencils on the donor grids. In determining the
best donor grid for a given receiver, both the local grid resolution and edge distance are used.
If a suitable donor cannot be found, the point will be marked as an orphan and a warning message
will be produced.

Relevant options:

* **`ConnectionType(m,n)`** - how donors on grid `m` will communicate with receivers on grid `n`.
  Possible values are **`OVK_CONNECTION_NONE`**, meaning communication is not allowed, and
  **`OVK_CONNECTION_NEAREST`**, **`OVK_CONNECTION_LINEAR`**, and **`OVK_CONNECTION_CUBIC`** for
  nearest-neighbor, linear, and cubic interpolation respectively _(default=`OVK_CONNECTION_NONE`)_
* **`EdgePadding(m,n)`** - (defined above; used again here as a threshold for deciding whether to
  prioritize edge distance or resolution).

## Mesh remapping

In addition to overset meshes, Overkit can also generate interpolation data for multi-block or
overset mesh remapping. For this application one would create a single domain containing both the
source and destination grids, and set assembly options as follows (omitting options left at the
default value):

* `Overlappable(<srcgrids>,<destgrids>)=.true.`
* `Occludes(<srcgrids>,<destgrids>)=OVK_OCCLUDES_ALL`
* `ConnectionType(<srcgrids>,<destgrids>)=<desired interp scheme>`.

## Interpolation data access

Once assembly is completed, interpolation data can be accessed using the following commands:

```fortran
  ! Get the connectivity data structure from the domain (returns pointer to ovk_connectivity)
  call ovkGetConnectivity(Domain, <gridid1>, <gridid2>, Connectivity)
  
  ! Access the donor extents (returns pointer to integer array of shape (3,2,<numconnections>))
  call ovkGetConnectivityDonorExtents(Connectivity, DonorExtents)
  
  ! Access the donor interpolation coefficients (returns pointer to real(kind=ovk_rk) array of
  ! shape (<stencilsize>,<numdims>,<numconnections>))
  call ovkGetConnectivityDonorInterpCoefs(Connectivity, DonorInterpCoefs)
  
  ! Access the donor coordinates (returns pointer to real(kind=ovk_rk) array of shape
  ! (<numdims>,<numconnections>))
  call ovkGetConnectivityDonorCoords(Connectivity, DonorCoords)
  
  ! Access the receiver points (returns pointer to integer array of shape (3,<numconnections>))
  call ovkGetConnectivityReceiverPoints(Connectivity, ReceiverPoints)
```

## Troubleshooting

Occasionally after running the assembly, you may get a result that doesn't match what you expected.
In these instances it can be helpful to look at the intermediate results of the steps described
above. For this purpose Overkit records some of these results in the grid state and provides a
command `ovkFilterGridState` to access them:

```fortran
  ! Get the grid from the domain
  call ovkGetGrid(Domain, <gridid>, Grid)
  
  ! Returns ovk_field_logical
  call ovkFilterGridState(Grid, <stateflags>, <filtermode>, Mask)
```

`<stateflags>` consists of one or more state flags combined using Fortran's bitwise 'or' (`ior`)

`<filtermode>` tells the command how to treat these flags. **`OVK_ANY`** returns points matching
any of the specified flags and **`OVK_ALL`** returns points matching all of them. **`OVK_NONE`**
and **`OVK_NOT_ALL`** are the respective complements of `OVK_ANY` and `OVK_ALL`.

Currently-supported state flags are:

* **`OVK_STATE_GRID`** - point is part of the domain
* **`OVK_STATE_INTERIOR`** - point is on the interior
* **`OVK_STATE_BOUNDARY`** - point is a boundary
* **`OVK_STATE_EXTERIOR`** - point is outside of the domain
* **`OVK_STATE_DOMAIN_BOUNDARY`** - point is a domain boundary
* **`OVK_STATE_INTERNAL_BOUNDARY`** - point is an internal boundary
* **`OVK_STATE_OVERLAPPED`** - point is overlapped by another grid
* **`OVK_STATE_INFERRED_DOMAIN_BOUNDARY`** - point is a domain boundary that was detected during
  boundary inference
* **`OVK_STATE_BOUNDARY_HOLE`** - point was removed during boundary hole cutting
* **`OVK_STATE_OCCLUDED`** - point is occluded by another grid
* **`OVK_STATE_FRINGE`** - point is part of the fringe
* **`OVK_STATE_OUTER_FRINGE`** - point is part of the outer fringe
* **`OVK_STATE_INNER_FRINGE`** - point is part of the inner fringe
* **`OVK_STATE_OVERLAP_MINIMIZED`** - point was removed during overlap minimization
* **`OVK_STATE_RECEIVER`** - point is a receiver
* **`OVK_STATE_ORPHAN`** - point was a receiver but no suitable donors were found.

Overkit also provides commands to write the grids and states to a file (PLOT3D format) that can be
opened in visualization tools such as TecPlot and ParaView. See the examples' source code
for details.

# Citing

[![DOI](https://zenodo.org/badge/82192918.svg)](https://zenodo.org/badge/latestdoi/82192918)

```tex
@misc{overkitfortran,
  author = {Smith, M.},
  title = {Overkit},
  year = {2019},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/majosm/overkit-fortran}},
  doi = {<get from Zenodo link above>}
}
```

# Acknowledgements

Overkit is one of several software projects produced by the University of Illinois Center for
Exascale Simulation of Plasma-Coupled Combustion (XPACC): <http://xpacc.illinois.edu>.

This material is based in part upon work supported by the Department of Energy, National Nuclear
Security Administration, under Award Number DE-NA0002374.

These examples show the functionality of the pv-atmos package.

- example_flat.py will construct a rectangular box with zonal wind data and a full labeled grid.
- example_sphere.py will construct a spherical geometry with "shells" to mimic of radial grid points.
- example_ocean.py will construct a rectangular ocean basin, including all continents, and show an isosurface of oxygen concentration.
- example_pole.py will project topographic and atmospheric data onto the North pole, and show isosurfaces of geopotential height during the Stratospheric Sudden Warming of Feb 23 1979.

When running the examples, be sure to set the variable pvAtmosPath to the path of the package files basic.py and grids.py and the examples folder.

Note that in the interest of smaller files, all data has been regridded to quite low resolution.

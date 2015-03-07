# pv_atmos 

Python scripting for scientific visualization software
[ParaView](http://www.paraview.org). Historically, pv_atmos has been developed to work with geophysical, and in particular, atmospheric model data (hence the name). However, pv_atmos has evolved into a very general package, and contains
routines for visualizing netCDF data, and the capability to show arbitrary axes and labels in a large variety of geometries (linear and logarithmic axes, spherical geometry).

This package is described in an open access peer-reviewed article:
Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
Journal of Open Research Software 2(1):e4, DOI: [http://dx.doi.org/10.5334/jors.al](http://dx.doi.org/10.5334/jors.al).
Please cite this work if you use this software for your publications.

## Components

Components are briefly described below. For more information on each function, please use 
`help(.)`. A list of all function within the modules is provided here.

### pv_atmos.basic

Provides functionality to read data on a 2D or 3D linear or logarithmic coordinates grid, including time evolution (if present) from a netCDF file. The netCDF should loosely correspond to the [Climate and Forecast (FC) conventions](https://en.wikipedia.org/wiki/Climate_and_Forecast_Metadata_Conventions). The important attribute is the time coordinate: ParaView will be looking for the "units: xxxx since xxxx" attribute to decide which dimension corresponds to time.

All functions are written in CamelCase, and variables in camelCase (sorry, it seems that both versions refer to the same animal).


#### Functions

Higher-order functions that you might want to use regularly are:

- LoadData()
- Cart2Spherical()
- Make3D()
- TransformCoords()
- MakeSelectable()
- DeleteAll()
- HideAll()
- ShowAll()
- CartWind2Sphere()
- Sphere2xyz()
- xyz2Sphere()

Helper functions for the above to work are:

- Cart2Log()
- GridAspectRatio()
- ConvertLogCoordString()
- ExtractBounds()


##### LoadData()

Read a netCDF file, convert linear to logarithmic (e.g. pressure to log-pressure) coordinate if desired, and transform according to prefered aspect ratio.

##### Cart2Spherical()

Transform rectangular geometry into a sphere with given radius.

##### Make3D()

Take a dataset with (a) 2D variable(s), and expand the chosen variable as third dimension. Classic example: Ocean bathymetry on a lon-lat grid. 

##### TransformCoords()

If not already done when loading the data, apply coordinate transformation in Cartesian coordinates, according to specified aspect ratio and logarithmic coordinates.

##### MakeSelectable()

In order to be able to switch any filter's visibility on/off in the GUI's pipeline, call this helper function.

##### DeleteAll(), HideAll(), ShowAll()

Delete, hide, or show all filters present in the pipeline.

##### CartWind2Sphere()

Converts zonal and meridional winds (or any velocity) from m/s into degrees longitude per time step and degrees latitude per time step. It can also convert pressure velocity
from hPa/s into the new vertical coordinate measure per time step.

##### Sphere2xyz(), xyz2Sphere()

Convert a given point in spherical (Cartesian) coordinates into Cartesian (spherical) coordinates, given the transformations applied to the data. Helpful to position labels, camera, etc.

### pv_atmos.grids

Provides the possibility to add axes, grid lines, planes (cuts), and labels. In
case of spherical geometry, one can also add shells, which are spheres
of a radius corresponding to a given vertical level. Planes
and shells contain data information, and can therefore be used for data
analysis as well as grid information.

These routines are not limited to any kind of data, and can be used with
any data, or even without data, to add a custom grid to a visualization. 

#### Functions

Higher-order functions that you might want to use regularly are:

- AddGrid()
- SphericalShells()
- AddGridPlane()
- AddGridLabel()
- AddAxisLabel()
- SphericalLables()
- WaterMark()

Helper functions for the above to work are:

- Lin2Log()
- BoundAspectRatio()

##### AddGrid()

Add a full grid, including grid lines at custom levels of all dimensions. This includes the appropriate lables of the grid lines, and labeling the axes.

##### SphericalShells()

Similar to AddGrid() in Cartesian coordinates, this adds shells around a sphere to serve as grid. These shells are labeled with the appropriate level value, and a water mark can be added to the outermost shell.

##### AddGridPlane()

Add one grid plane along one dimension.

##### AddGridLabel()

Add one label along one dimension.

##### AddAxisLabel()

Label one given axis

##### SphericalLabels()

Label any number of vertical levels in spherical geometry.

##### WaterMark()

Add a water mark to one of the spherical shells. Nice to brand your viz.

##### LonLat2Polar()

Project 2D longitude-latitude or 3D longitude-latitude-vertical data onto polar coordinates around the North or South pole. This has a little 3D twist, in that the projection can be domed in the vertical.


# Releases

See releases with changelogs in the releases panel of the GitHub
distribution. Version 1.0 corresponds to the description of Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView. Journal of Open Research Software 2(1):e4, DOI: [http://dx.doi.org/10.5334/jors.al](http://dx.doi.org/10.5334/jors.al).


# Installation

This is a regular python package, which can be installed with [pip](https://pypi.python.org/pypi/pip).
However, be advised that the functions depend on paraview.simple, which is not available as independent python package. See below how to run pv_atmos once it's installed.


1) Using pip: `pip install pv_atmos`

2) Manual python package install: Download the .zip file from this repository, and unpack it. Run `python setup.py install` for installation.

3) For use in the ParaView python console without installation: No python installation and/or
command shell is needed. Download the .zip file, unpack it where
convenient. Start ParaView, and open the Python Shell contained within
ParaView. If the pv_atmos files are not unpacked in the run directory,
use the `Run Script` button and choose either `basic.py` or `grids.py` (or one after the other of course), and you are ready to use the pv_atmos functions.

# Use

Follow any of the below bullet points to get going with pv_atmos.

* No Python outside ParaView is needed, as ParaView ships with its own distribution. If you don't want to use full python in a console, but simply want to work with the ParaView GUI:

  - Click on `Run Script`, and double-click on `basic.py` and `grids.py`
  - Or: Open `Tools -> Python Shell`, and then:
```
  $ from pv_atmos.basic import *
  $ from pv_atmos.grids import *
```

* Run the version of python shipped with ParaView: This will automatically adjust your python path, and paraview.simple will be recognized:
```
$ /Applications/paraview.app/Contents/bin/pvpython
$ from pv_atmos.basic import *
$ from pv_atmos.grids import *
```
* Set the PYTHONPATH to where paraview.simple resides. On a Mac, this is typically
```
$ export DYLD_FALLBACK_LIBRARY_PATH="/Applications/paraview.app/Contents/Libraries"
$ export LD_LIBRARY_PATH="/Applications/paraview.app/Contents/Libraries"
$ export DYLD_FALLBACK_FRAMEWORK_PATH="/Applications/paraview.app/Contents/Frameworks"
$ export PYTHONPATH="/Applications/paraview.app/Contents/Python:/Applications/paraview.app/Contents/Libraries
$ python
$ from pv_atmos.basic import *
$ from pv_atmos.grids import *
```

# Examples

The `examples` directory contains four example scripts and the data files
uv_daily.nc, ocean_depth.nc, ECMWF_19790223.nc, and ocean_o2.nc. The examples can be run within the python terminal of ParaView, or a general python session, provided paraview.simple is located in the python path.
The example files contain the 3D structure of zonal and meridional wind over three daily time steps, created from GCM output; ocean topography data from GFDL's CM2.1 model; geopotential height anomalies from ERA-Interim reanalysis; and oxygen data from GFDL's ESM2M model, provided by Thomas Froelicher. One script will create a spherical, one a rectangular plot of zonal wind. The pole script will project reanalysis and bathymetry data onto the North pole. The ocean script will create a rectangular ocean basin.
When using the example files, make sure to set `pvAtmosPath` within the example scripts or the command line to the directory containing `basic.py` and `grids.py`:
```
$ pvAtmosPath = '/path/to/pv_atmos/'
```
The examples will import pv_atmos themselves, so no loading of pv_atmos is necessary prior to running the examples.

# Dependencies

Needs the python module `math` for coordinate conversion (Pi and log), and of course `paraview.simple` (which comes with ParaView).

# Remarks

Input netCDF files should generally conform to the Climate and Forecast (CF) metadata convention. A not so rigorous but sufficient test is to load the file manually into ParaView, and select "CF" convention; if ParaView reads in the data, the here included scripts can be used without problems. For instance, the time coordinate does not have to conform to CF conventions, in particular the "calendar" attribute: The important thing to make ParaView accept "time" as time is the attribute "units", which must contain the word "since".

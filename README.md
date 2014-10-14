pv_atmos 
========

Python scripting for scientific visualization software
[ParaView](http://www.paraview.org). In particular, pv_atmos contains
routines for visualizing data from geophysical netCDF data, and the capability to show arbitrary axes and labels (linear and logarithmic axes).

This package is described in an open access peer-reviewed article:
Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
Journal of Open Research Software 2(1):e4, DOI: [http://dx.doi.org/10.5334/jors.al](http://dx.doi.org/10.5334/jors.al).
Please cite this work if you use this software for your publications.

No Python outside ParaView is needed, as ParaView ships with its own
distribution.

If, on the other hand, pv_atmos shall be used as a python package
outside the ParaView python console, make sure paraview.simple is in the
python path (see Installation & Use below).

### atmos_basic
### 
Provides functionality to read data on a 2D or 3D linear or logarithmic coordinates grid, including time evolution
(if present) from a netCDF file. The netCDF should loosely correspond to
the [Climate and Forecast (FC)
conventions](https://en.wikipedia.org/wiki/
Climate_and_Forecast_Metadata_Conventions). The important attribute is
the time coordinate: ParaView will be looking for the "units: xxxx since
xxxx" attribute to decide which dimension corresponds to time.

The most important function is "LoadData", which will read the netCDF
file, and convert linear to logarithmic (e.g. pressure to log-pressure) coordinate if desired. In
addition, "Cart2Spherical" will transform the rectangular geometry into
a sphere with given radius, and "CartWind2Atmos" converts zonal and
meridional winds from m/s into degrees longitude per time step and
degrees latitude per time step. It can also convert pressure velocity
from hPa/s into the new vertical coordinate measure per time step.

### atmos_grids
### 
Provides the possibility to add axes, grid lines, planes, and labels. In
case of spherical geometry, one can also add shells, which are spheres
of a radius corresponding to a given pressure or height level. Planes
and shells contain data information, and can therefore be used for data
analysis as well as grid information.

These routines are not limited to any kind of data, and can be used with
any data, or even without data, to add a custom grid to a visualization. By default, the axes are named 'lon', 'lat', and 'pressure [hPa]', but this can be changed with the variable 'AxisNames' in the function 'AddGrid'.


Releases
--------

See releases with changelogs in the releases panel of the GitHub
distribution.

Installation & Use
------------------

This is a regular python package, which can be installed with [pip](https://pypi.python.org/pypi/pip).
However, be advised that the functions depend on paraview.simple, which is not available as independent python package.
Nevertheless, once installed, pv_atmos can be imported inside the ParaView python console.

1) Using pip: pip install pv_atmos

2) Manual python package install: Download the .zip file from this repository, and unpack it. Run 'python setup.py install' for installation.

With any of these methods, you can now load pv_atmos in the ParaView python console like this:

>>> from pv_atmos.atmos_basic import *
>>> from pv_atmos.atmos_grids import *

or:

>>> import pv_atmos.atmos_basic as ab
>>> import pv_atmos.atmos_grids as ag

3) For use in the ParaView python console without installation: No python installation and/or
command shell is needed. Download the .zip file, unpack it where
convenient. Start ParaView, and open the Python Shell contained within
ParaView. If the pv_atmos files are not unpacked in the run directory,
use the "run script" button and choose atmos_basic.py first,
atmos_grids.py second, and you are ready to use the pv_atmos functions.


Examples
--------

The examples directory contains three example scripts and the data files
uv_daily.nc, ocean_depth.nc, and ocean_o2.nc. The examples can be run within the python terminal of ParaView, or a general python session, provided paraview.simple is located in the python path.
The example file contain the 3D structure of zonal and meridional wind
over three daily time steps, created from GCM output; ocean topography data from GFDL's CM2.1 model, and oxygen data from GFDL's ESM2M model, provided by Thomas Froelicher. One script will
create a spherical, one a rectangular plot of zonal wind. The ocean script will create a rectangular ocean basin.
When using the
example files, make sure to set "pvAtmosPath"
within the example scripts to the directory containing atmos_basic.py
and atmos_grids.py.


Dependencies
------------

Needs the python module math for coordinate conversion (Pi and log).
atmos_grids needs atmos_basic, atmos_basic needs paraview.simple.

Remarks
-------

Input netCDF files should generally conform to the Climate and Forecast
(CF) metadata convention. A not so rigorous but sufficient test is to
load the file manually into ParaView, and select "CF" convention; if
ParaView reads in the data, the here included scripts can be used
without problems. For instance, the time coordinate does not have to
conform to CF conventions, in particular the "calendar" attribute: The
important thing to make ParaView accept "time" as time is the attribute
"units", which must contain the word "since".

pv-atmos
========

Python scripting for scientific visualization software [ParaView](http://www.paraview.org). In particular, pv-atmos contains routines for visualizing data from geophysical netCDF data.

atmos_basic
-----------

Provides functionality to read data on a latitude - longitude and, if desired, pressure or height coordinates grid, including time evolution (if present) from a netCDF file. The netCDF should loosely correspond to the [Climate and Forecast (FC) conventions](https://en.wikipedia.org/wiki/Climate_and_Forecast_Metadata_Conventions). The important attribute is the time coordinate: ParaView will be looking for the "unit: xxxx since xxxx" attribute to decide which dimension corresponds to time.

The most important function is "loadData", which will read the netCDF file, and convert pressure to log-pressure corrdinate if desired. In addition, "Cart2Spherical" will transform the rectangular geometry into a sphere with given radius, and "CartWind2Atmos" converts zonal and meridional winds from m/s into degrees longitude per time step and degrees latitude per time step. It can also convert pressure velocity from hPa/s into the new vertical coordinate measure per time step.

atmos_grids
-----------

Provides the possibility to add axes, grid lines, planes, and labels. In case of spherical geometry, one can also add shells, which are spheres of a radius corresponding to a given pressure or height level. Planes and shells contain data information, and can therefore be used for data analysis as well as grid information.

These routines are not limited to any kind of data, and can be used with any data, or even without data, to add a custom grid to a visualization. 


Examples
--------

The examples directory contains two example scripts and the data file uv_daily.nc, which can be run within the python terminal of ParaView. One will create a spherical, one a rectangular plot of zonal wind.


Dependencies
------------

Needs the python module math for coordinate conversion (Pi and log). atmos_grids needs atmos_basic.

Remarks
-------

Input netCDF files should generally conform to the Climate and Forecast (CF) metadata convention. A not so rigorous but sufficient test is to load the file manually into ParaView, and select "CF" convention; if ParaView reads in the data, the here included scripts can be used without problems.
For instance, the time coordinate does not have to conform to CF conventions, in particular the "calendar" attribute: The important thing to make ParaView accept "time" as time is the attribute "units", which must contain the word "since".
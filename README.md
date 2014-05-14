pv-atmos
========

Python scripting for scientific visualization software ParaView. Applied to atmospheric netCDF data.

Input netCDF files should generally conform to the Climate and Forecast (CF) metadata convention. A not so rigorous but sufficient test is to load the file manually into ParaView, and select "CF" convention; if ParaView reads in the data, the here included scripts can be used without problems.
For instance, the time coordinate does not have to conform to CF conventions, in particular the "calendar" attribute: The important thing to make ParaView accept "time" as time is the attribute "units", which must contain the word "since".
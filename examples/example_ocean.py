#!/usr/bin/python
# Filename: example_ocean.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al

######################################################################################
# set path to the locations of grids.py, basic.py, and the examples folder
# you can define this in the session or here
try:
    pvAtmosPath
except:
    pvAtmosPath = '../'
try: #is pv_atmos installed?
    from pv_atmos.basic import *
    from pv_atmos.grids import *
except:
    execfile(pvAtmosPath + 'basic.py')
    execfile(pvAtmosPath + 'grids.py')

## show me where the files are, relative to pvAtmosPath ##
# path to ocean files
oceanPath = pvAtmosPath + 'examples/'

# file containing bathymetry
topoFile = 'ocean_depth.nc'
topoDims = ['rlon','rlat']
depthVar = 'deptho'

# file containing oxygen data
dataFile = 'ocean_o2.nc'
dataDims = ['xt_ocean','yt_ocean','st_ocean']
# the values we will be interested in
dataName = 'o2'
dataContours = [8e-5]

## how would you like the transformation to work ##
logCoord = [] #no logarithmic coordinates
aspRat   = [1,1,0.01] # divide bathymetry by 100


### get the data ###

# check if files exist
import os.path
if not os.path.isfile(oceanPath+topoFile):
    raise ValueError, oceanPath+topoFile+' does not exist!'
if not os.path.isfile(oceanPath+dataFile):
    raise ValueError, oceanPath+dataFile+' does not exist!'


# bathymetry
(depth_out,depth_coor)=LoadData(oceanPath+topoFile,ncDims=topoDims,logCoords=logCoord )
# get the bounds of the topography. This is important if bathymetry and data files have different origins
topoBds = depth_out.GetDataInformation().GetBounds()

# data
(o2_out,o2_coor)=LoadData(oceanPath+dataFile,ncDims=dataDims,aspectRatios=aspRat,logCoords=logCoord )
# we want to replace the fill values with NaNs here
o2_out.ReplaceFillValueWithNan = 1
# get the bounds of the data file. Again important if different from bathymetry files
dataBds = o2_out.GetDataInformation().GetBounds()
# instead of NaNs, there are -1e10 values in this file that we don't want
o2_thresh = Threshold(o2_coor,ThresholdRange=[0,1])
MakeSelectable()

# the bathymetry file is in cell data, need to convert to point data
c2p=CellDatatoPointData(depth_coor)
MakeSelectable()

# now make the bathymetry 3D
bathy = Make3D(depthVar, expandDir='-z', aspectRatios=aspRat, logCoords=logCoord, src=c2p)
rep=Show(bathy)
rep.ColorArrayName = depthVar
# this works only for paraview > 4.1
try:
    depthVal = bathy.PointData.GetArray(depthVar)
    lkpD = AssignLookupTable(depthVal,'erdc_blue_BW')
    #invert the colors
    valPts=lkpD.RGBPoints[::4]
    lkpD.RGBPoints[::4]=valPts[::-1]
    rep.LookupTable = lkpD
except:
    pass


#### now add data to the ocean #######
o2_cont = Contour(o2_thresh, ContourBy=['POINTS',dataName], Isosurfaces=dataContours)
rep=Show()
rep.ColorArrayName = dataName
try:
    dataVal = o2_cont.PointData.GetArray(dataName)
    lkpW = AssignLookupTable(dataVal,'GnYlRd')
    rep.LookupTable = lkpW
except:
    pass

#### finally, add a grid  ####
AddGrid(xlevels=[-270,-225,-180,-135,-90,-45,0,45], ylevels=[-60,-30,0,30,60], zlevels=range(-5000,1000,1000), bounds=[-280,80,-90,90,-5500,0], ratios=aspRat, logCoord=logCoord, AxisNames=["longitude","latitude","depth [m]"], AxisColor=[0,0,0], AxisWidth=1.0,LabelSize=5.0)


Render()

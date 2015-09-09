#!/usr/bin/python
# Filename: example_pole.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al

######################################################################################
# set path to the locations of pv_atmos/grids.py, pv_atmos/basic.py, and the examples folder
# you can define this in the session or here
try:
    pvAtmosPath
except:
    pvAtmosPath = raw_input("Please provide the path where pv_atmos lives: ")
    pvAtmosPath = pvAtmosPath+'/'
try: #is pv_atmos installed?
    from pv_atmos.basic import *
    from pv_atmos.grids import *
except:
    execfile(pvAtmosPath + 'basic.py')
    execfile(pvAtmosPath + 'grids.py')


## there will be two files: one with Earth's surface, one with atmospheric data

# polar projection parameters:
# we want to project the full Earth
cut = -90
# - this is for illustration, usually you want less -
# and map it to a quarter of a sphere
alpha = 0.25


# first, look at Earth's surface:
# path to example files
filePath = pvAtmosPath + 'examples/'
# file containing bathymetry
topoFile = 'ocean_depth.nc'
topoDims = ['rlon','rlat']
depthVar = 'deptho'
topoLogCoord = [] #no logarithmic coordinates
# correct for missing degrees in longitude - this is cheating
# divide bathymetry by 1e5
topoAspRat   = [360./356.,1,2e-5]

# then, the file containing geopotential height data
dataFile = 'ECMWF_19790223.nc'
dataDims = ['longitude','latitude','level']
# the values we will be interested in
dataName = 'height'
dataContours = [-10e3,-7.5e3,-5e3,7.5e3,10e3,12.5e3,15e3,17.5e3,20e3]
# define aspect ratio of coordinate system: keep lon and lat, and log-p. Base log-p on 1000hPa
dataAspRat = [1,1,0.5]
dataLogCoord = [2]
basis = [1e3]

## now load example data file ##
topoFileName = filePath + topoFile
dataFileName = filePath + dataFile
# check if files exist
import os.path
if not os.path.isfile(topoFileName):
    raise ValueError, topoFileName+' does not exist!'
if not os.path.isfile(dataFileName):
    raise ValueError, dataFileName+' does not exist!'

## Earth

# get bathymetry: this is exactly the same as in example_ocean.py
(depth_out,depth_coor)=LoadData(topoFileName,ncDims=topoDims,logCoords=topoLogCoord, replaceNaN = False )
# the bathymetry file is in cell data, need to convert to point data
c2p=CellDatatoPointData(depth_coor)
MakeSelectable()

# now make the bathymetry 3D
bathy = Make3D(depthVar, expandDir='-z', aspectRatios=topoAspRat, logCoords=topoLogCoord, src=c2p)

# here is where we apply the polar projection
# we only keep the Northern Hemisphere
polBathy = LonLat2Polar(alpha=alpha, src=bathy, cutLat=cut)
# now some colors
repD=Show(polBathy)
repD.ColorArrayName = depthVar
# this works only for paraview > 4.1
try:
    depthVal = bathy.PointData.GetArray(depthVar)
    lkpD = AssignLookupTable(depthVal,'erdc_blue_BW')
    #invert the colors
    valPts=lkpD.RGBPoints[::4]
    lkpD.RGBPoints[::4]=valPts[::-1]
    repD.LookupTable = lkpD
except:
    pass


## Atmosphere
(h_out,h_coor) = LoadData(dataFileName,ncDims=dataDims,aspectRatios=dataAspRat,logCoords=dataLogCoord,basis=basis)
# project this onto the pole
polAtmos = LonLat2Polar(alpha=alpha, src=h_coor, cutLat=cut)

# extract the contours defined above
polConts = Contour(polAtmos, ContourBy=['POINTS',dataName], Isosurfaces=dataContours)
repH=Show(polConts)
repH.ColorArrayName = dataName
# this works only for paraview > 4.1
try:
    dataVal = polAtmos.PointData.GetArray(dataName)
    lkpH = AssignLookupTable(dataVal,'BuRd')
    repH.LookupTable = lkpH
except:
    pass
# change representation type
repH.Representation = 'Points'

# now set the camera
renderView = GetActiveView()
renderView.CameraPosition = [11.991470393200085, -5.743012400242967, 6.854643686759562]
renderView.CameraFocalPoint = [-0.22555362372853732, 0.18852991370207475, 0.41220680834358314]
renderView.CameraViewUp = [-0.4101252214668251, 0.13495605517576154, 0.9019890054142482]
renderView.CameraParallelScale = 3.9563125825408108

## now, let's do a full, flat polar projection of the northern hemisphere
# get layout
viewLayout = GetLayout()
# split cell
viewLayout.SplitHorizontal(0, 0.5)
# Create a new 'Render View'
renderView = CreateView('RenderView')
# place view in the layout
viewLayout.AssignView(2, renderView)

polTopoFlat = LonLat2Polar(alpha=0, cutLat=0, src=bathy)
repDFlat = Show(polTopoFlat)
repDFlat.ColorArrayName = depthVar
try:
    repDFlat.LookupTable = lkpD
except:
    pass

# same thing for the data
polAtmosFlat = LonLat2Polar(alpha=0, src=h_coor, cutLat=0)

# extract the contours defined above
polContsFlat = Contour(polAtmosFlat, ContourBy=['POINTS',dataName], Isosurfaces=dataContours)
repHFlat=Show(polContsFlat)
repHFlat.ColorArrayName = dataName
# this works only for paraview > 4.1
try:
    repHFlat.LookupTable = lkpH
except:
    pass
# change representation type
repHFlat.Representation = 'Points'
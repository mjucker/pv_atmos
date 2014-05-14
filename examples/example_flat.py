#!/usr/bin/python

## import paraview ##
try: paraview.simple
except: from paraview.simple import *
## other needed modules from standard python ##
import math
## load atmos modules ##
# set path to the locations of atmos_grids.py and atmos_basic.py
pvAtmosPath = '../'
try: from atmos_basic import *
except: execfile(pvAtmosPath + 'atmos_basic.py')
try: from atmos_grids import *
except: execfile(pvAtmosPath + 'atmos_grids.py')

# define aspect ratio of coordinate system: keep lon and lat, multiply log-p by 20. Base log-p on 1e3hPa
ratio = [1,1,20]
basis = 1e3

## now load example file ##
fileName = pvAtmosPath + 'examples/uv_daily.nc'
# the file is 3D+time, in pressure coordinates, and we adjust the axis aspect ratio
(output_nc,CorrZ,Coor,AspRat) = loadData(fileName, ['pfull','lat','lon'], 1, ratio)
# we have now read in a 4D file, with log-pressure in the Z-direction

# as an example, apply contour filter on zonal wind to show tropospheric jet streams and stratospheric polar vortex
# then, add contour at 25m/s
Cont = Contour(ContourBy='ucomp',Isosurfaces=[25])
# color it
repU = GetDisplayProperties()
repU.ColorArrayName = 'ucomp'
# get lookup table for coloring ucomp, this only works in >v4.0
try:
    ucomp = output_nc.PointData.GetArray('ucomp')
    lkpU = AssignLookupTable(ucomp,'Cool to Warm')
    repU.LookupTable = lkpU
except:
    pass
# make it transparent
repU.Opacity = 0.7
Show()

# now we want to add arrows showing the horizontal wind vector
(W,normW,clipWS,clipWN) = CartWind2Atmos(src=AspRat,ratios=ratio)
# add arrows
Arrows = Glyph()
Arrows.Scalars = ['POINTS','normW']
Arrows.Vectors = ['POINTS','W']
Arrows.ScaleMode = 'scalar'
# often the glyphs are too large by default. make them smaller
Arrows.SetScaleFactor = Arrows.SetScaleFactor/5
# for animation, might be better to keep the same points at all time steps
Arrows.KeepRandomPoints = 1
repA = Show()
repA.ColorArrayName = 'GlyphVector'
# also, create a colormap lookup table for future use. Again, only >v4.0
try:
    wind = Arrows.PointData.GetArray('GlyphVector')
    lkpW = AssignLookupTable(wind,'X Ray')
    lkpW.ColorSpace = 'Lab'
    repA.LookupTable = lkpW
except:
    pass

## add axes and grid
LabSze = 6.0
Bounds = [0,360,-90,90,1e3,0.01]
# outside box
AddGrid(press=[1e3,0.1],lats=[-90,90],lons=[0,360],bounds=Bounds, ratios=ratio, basis=basis, AxisWidth=2.0,LabelSize=LabSze)
# inside grid lines: add only  line labels, not axis labels
AddGrid(press=[1e2,10,1],lats=[-45,0,45],lons=[90,180,270],bounds=Bounds,ratios=ratio, basis=basis, AxisWidth=1.0,LabelSize=-LabSze)

# add a plane showing wind stength at 100hPa
Plane100=AddPresPlane(100,Bounds,ratio,basis,1,normW)
RenameSource("Plane100",Plane100)
Plane100rep = GetDisplayProperties(Plane100)
Plane100rep.ColorArrayName = 'normW'
try:
    Plane100rep.LookupTable = lkpW
except:
    pass
Plane100rep.Opacity = 0.7
# add a plane showing zonal wind strength at 180 longitude
Plane180=AddLonPlane(180,Bounds,ratio,basis,1,normW)
RenameSource("Plane180E",Plane180)
Plane180rep = GetDisplayProperties(Plane180)
Plane180rep.ColorArrayName = 'ucomp'
try:
    Plane180rep.LookupTable = lkpU
except:
    pass
Plane180rep.Opacity = 0.7
# add a plane showing meridional wind strength at the equator
PlaneEQ=AddLatPlane(0,Bounds,ratio,basis,1,normW)
RenameSource("PlaneEQ",PlaneEQ)
PlaneEQrep = GetDisplayProperties(PlaneEQ)
PlaneEQrep.ColorArrayName = 'vcomp'
try:
    vcomp = output_nc.PointData.GetArray('vcomp')
    lkpV = AssignLookupTable(vcomp,'bone_Matlab')
    PlaneEQrep.LookupTable = lkpV
except:
    pass
PlaneEQrep.Opacity = 0.7

# some camera settings
(Left,Right,Near,Far,Bottom,Top) = BoundAspectRatio(Bounds,ratio,basis)
view=GetActiveView()
view.CameraPosition = [2*Right,0.5*(Near+Far),Top*1.2]
view.CameraFocalPoint = [0.5*(Right+Left),0.5*(Near+Far),0.5*(Bottom+Top)]
view.CameraParallelProjection = 0
view.CameraViewAngle = 30
view.CameraViewUp = [0,0,1]

Render()



#!/usr/bin/python
# Filename: example_flat.py
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
# define aspect ratio of coordinate system: keep lon and lat, multiply log-p by 20. Base log-p on 1e3hPa
ratio = [1,1,20]
logCoord = [2]
basis = [1e3]

## now load example file ##
fileName = pvAtmosPath + 'examples/uv_daily.nc'
# check if files exist
import os.path
if not os.path.isfile(fileName):
    raise ValueError, fileName+' does not exist!'
# the file is 3D+time, in pressure coordinates, and we adjust the axis aspect ratio
(output_nc,Coor) = LoadData(fileName, ['lon','lat','pfull'], ratio, logCoord, basis)
# we have now read in a 4D file, with log-pressure in the Z-direction

# as an example, apply contour filter on zonal wind to show tropospheric jet streams and stratospheric polar vortex
# then, add contour at 25m/s
Cont = Contour(ContourBy='ucomp',Isosurfaces=[25])
# color it
repU = Show()
repU.ColorArrayName = 'ucomp'
# get lookup table for coloring ucomp, this only works in ParaView >= v4.1
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
(W,normW,clipWS,clipWN) = CartWind2Sphere(src=Coor,ratios=ratio)
# add arrows
Arrows = Glyph()
Arrows.Scalars = ['POINTS','normW']
Arrows.Vectors = ['POINTS','W']
Arrows.ScaleMode = 'scalar'
# often the glyphs are too large by default. make them smaller
Arrows.ScaleFactor = 0.5
# make sure seed is the same for reproducibility
Arrows.Seed = 10339
repA = Show()
repA.ColorArrayName = 'GlyphVector'
# also, create a colormap lookup table for future use. Again, only >= v4.1
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
# outside box, rename axes
AddGrid(xlevels=[0,360],ylevels=[-90,90],zlevels=[1e3,0.1],bounds=Bounds, ratios=ratio, logCoord=logCoord, basis=basis, AxisWidth=2.0,LabelSize=LabSze, AxisNames=["longitude","latitude","pressure [hPa]"])
# inside grid lines: add only line labels, not axis labels
AddGrid(xlevels=[90,180,270],ylevels=[-45,0,45],zlevels=[1e2,10,1],bounds=Bounds,ratios=ratio, logCoord=logCoord, basis=basis, AxisWidth=1.0,LabelSize=-LabSze)

# add a plane showing wind stength at 100hPa
Plane100=AddGridPlane(2, 100,Bounds,ratio,logCoord,basis,1,normW)
RenameSource("Plane100",Plane100)
Plane100rep = GetDisplayProperties(Plane100)
Plane100rep.ColorArrayName = 'normW'
try:
    Plane100rep.LookupTable = lkpW
except:
    pass
Plane100rep.Opacity = 0.7
# add a plane showing zonal wind strength at 180 longitude
Plane180=AddGridPlane(0, 180,Bounds,ratio,logCoord,basis,1,normW)
RenameSource("Plane180E",Plane180)
Plane180rep = GetDisplayProperties(Plane180)
Plane180rep.ColorArrayName = 'ucomp'
try:
    Plane180rep.LookupTable = lkpU
except:
    pass
Plane180rep.Opacity = 0.7
# add a plane showing meridional wind strength at the equator
PlaneEQ=AddGridPlane(1, 0,Bounds,ratio,logCoord,basis,1,normW)
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
(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax) = BoundAspectRatio(Bounds,ratio,logCoord,basis)
view=GetActiveView()
view.CameraPosition = [2*Xmin,0.5*(Ymin+Ymax),Zmax*1.2]
view.CameraFocalPoint = [0.5*(Xmin+Xmax),0.5*(Ymin+Ymax),0.5*(Zmin+Zmax)]
view.CameraParallelProjection = 0
view.CameraViewAngle = 30
view.CameraViewUp = [0,0,1]

Render()



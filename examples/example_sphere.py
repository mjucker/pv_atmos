#!/usr/bin/python
# Filename: example_sphere.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al

######################################################################################
# set path to the locations of atmos_grids.py, atmos_basic.py, and the examples folder
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

# define aspect ratio of coordinate system: keep lon and lat, and log-p. Base log-p on 1e3hPa
ratio = [1,1,1]
logCoord = [2]
basis = [1e3]
radius = 1

## now load example file ##
fileName = pvAtmosPath + 'examples/uv_daily.nc'
# check if files exist
import os.path
if not os.path.isfile(fileName):
    raise ValueError, fileName+' does not exist!'
# the file is 3D+time, in pressure coordinates, and we adjust the axis aspect ratio
(output_nc,Coor) = LoadData(fileName, ['lon','lat','pfull'], ratio, logCoord, basis)
# we have now read in a 4D file, with log-pressure in the Z-direction

# convert box into spherical coordinates
Globe = Cart2Spherical(radius,Coor)
MakeSelectable(Globe)

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


## add grid in form of 'shells' in the atmosphere. also add a watermark
Shells = SphericalShells(radius,ratio,logCoord,basis,Coor,[10,1],1,waterMark='bob 2014')

# also add a shell at 100hPa, colored by zonal wind
addShell = [100]
Plane100=SphericalShells(radius,ratio,logCoord,basis,Coor,addShell,1)
Srep = Show(Plane100[0])
Srep.ColorArrayName = 'ucomp'
try:
    Srep.LookupTable = lkpU
except:
    pass
Srep.Opacity = 0.5
# add a label for this plane
Label100 = SphericalLabels(radius,ratio,logCoord,basis,addShell)

# finally, add a sphere marking the surface
Surface = Sphere(Radius=radius)
RenameSource("Surface",Surface)
Surface.ThetaResolution = 128
Surface.PhiResolution   = 64
Show()


# some example for camera position
view=GetActiveView()
view.CameraPosition = [-15,-10,5.3]
view.CameraFocalPoint = [0,0,0.5]
view.CameraViewUp = [0,0,1]

Render()



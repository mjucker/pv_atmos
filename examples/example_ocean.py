
def transformTopo(src=GetActiveSource(),moveXFunction=''):
    depth = Calculator(src)
    depth.Function = 'iHat*(coordsX'+moveXFunction+') + jHat*coordsY - kHat*abs('+str(depthVar)+')'
    depth.CoordinateResults = 1
    MakeSelectable()
    
    aspect = Calculator(depth)
    aspect.Function = 'iHat*coordsX*'+str(aspRat[0])+' + jHat*coordsY*'+str(aspRat[1])+' + kHat*coordsZ*'+str(aspRat[2])
    aspect.CoordinateResults = 1
    rep = Show(aspect)
    rep.ColorArrayName = depthVar
    # assign a colormap lookup table. Only ParaView > v4.0
    try:
        depthVal = aspect.PointData.GetArray(depthVar)
        lkpD = AssignLookupTable(depthVal,'erdc_blue_BW')
        #invert the colors
        valPts=lkpD.RGBPoints[::4]
        lkpD.RGBPoints[::4]=valPts[::-1]
        rep.LookupTable = lkpD
    except:
        pass

#########

try:
    from atmos_basic import *
    from atmos_grids import *
except:
    pvAtmosPath='../'
    execfile(pvAtmosPath + 'atmos_basic.py')
    execfile(pvAtmosPath + 'atmos_grids.py')

## show me where the files are ##
oceanPath='./'

topoFile = 'ocean_depth.nc'
topoDims = ['rlon','rlat']
depthVar = 'deptho'

dataFile = 'ocean_o2.nc'
dataDims = ['xt_ocean','yt_ocean','st_ocean']
# the values we will be interested in
dataName = 'o2'
dataContours = [1e-4]

## how would you like the transformation to work ##
logCoord = [] #no logarithmic coordinates
aspRat   = [1,1,0.01]


### get the data ###

# topography
(depth_out,depth_coor)=LoadData(oceanPath+topoFile,ncDims=topoDims,logCoords=logCoord )
# get the bounds of the topography
topoBds = depth_out.GetDataInformation().GetBounds()

# data
(o2_out,o2_coor)=LoadData(oceanPath+dataFile,ncDims=dataDims,aspectRatios=aspRat,logCoords=logCoord )
# we want to replace the fill values with NaNs here
o2_out.ReplaceFillValueWithNan = 1
# get the bounds of the data file
dataBds = o2_out.GetDataInformation().GetBounds()
# instead of NaNs, there are -1e10 values in this file that we don't want
o2_thresh = Threshold(o2_coor,ThresholdRange=[0,1])

#the topography file is in cell data, need to convert to point data
c2p=CellDatatoPointData(depth_coor)
MakeSelectable()

### see if we have to move the topography file to align it with the data file
### this is not necessary with the provided files, but might be with other files
swapTopo = False
if topoBds[0] != dataBds[0] or topoBds[1] != dataBds[1]:
    if topoBds[0] < dataBds[0]:
        moveClip = Clip(c2p, ClipType="Plane")
        moveClip.ClipType.Origin=[dataBds[0],0,0]
        moveClip.ClipType.Normal=[-1,0,0]
        MakeSelectable()
        RenameSource('moveRight',moveClip)
        moveXFunction = '+ 360'
        stayClip = Clip(c2p, ClipType="Plane")
        stayClip.ClipType.Origin=[dataBds[0],0,0]
        stayClip.ClipType.Normal=[+1,0,0]
        MakeSelectable()
        RenameSource('stayClip',stayClip)
    elif topoBds[1] > dataBds[1]:
        moveClip = Clip(c2p, ClipType="Plane")
        moveClip.ClipType.Origin=[dataBds[1],0,0]
        moveClip.ClipType.Normal=[1,0,0]
        MakeSelectable()
        RenameSource('moveLeft',moveClip)
        moveXFunction = '- 360'
        stayClip = Clip(c2p, ClipType="Plane")
        stayClip.ClipType.Origin=[dataBds[0],0,0]
        stayClip.ClipType.Normal=[-1,0,0]
        MakeSelectable()
        RenameSource('stayClip',stayClip)
    swapTopo = True


### work on the topography #####

if swapTopo :
    transformTopo(moveClip,moveXFunction)
    transformTopo(stayClip,'')
else:
    transformTopo(c2p,'')


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


Render()

#!/usr/bin/python
# Filename: atmos_basic.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite CITATION HERE
#
# Python interface for ParaView (www.paraview.org). Reads netCDF file on a latitude - longitude and, if desired, pressure or height coordinates grid, including time evolution (if present). netCDF file needs to correspond to Climate and Forecast (FC) conventions (https://en.wikipedia.org/wiki/Climate_and_Forecast_Metadata_Conventions).
# Also Provides functions to modify Cartesian coordinates and wind components.

##### needed modules: paraview.simple, math #########################

# some global constants
strPi = str(math.pi)[0:7]

##### define auxiliary functions ##################################
# make sure Z-coordinate (pressure or height) has only positive values
def CorrectZCoord(src=GetActiveSource()):
    calc=Calculator(src)
    calc.Function = 'iHat*coordsX + jHat*coordsY + kHat*abs(coordsZ)'
    calc.CoordinateResults = 1
    return calc

# define pressure to z coordinate conversion
def ConvertPressureString(pString, ratio=1.0, basis=1e3):
    if ratio == 1.0:
        expression = 'log10(abs(' + pString + ')/' + str(basis) + ')'
    else:
        expression = 'log10(abs(' + pString + ')/' + str(basis) + ')*' + str(ratio)
    return expression

def Pressure2Z(plevel, ratio=1.0, basis=1e3):
    level = -math.log10(plevel/basis)*ratio
    return level

# convert pressure coordinates to Cartesian coordinates
def Pressure2Cart(src=GetActiveSource(), ratio=1.0, basis=1e3):
    pFun = ConvertPressureString('coordsZ',ratio,basis)
    calc=Calculator(src)
    calc.Function = 'iHat*coordsX + jHat*coordsY - kHat*'+pFun
    calc.CoordinateResults = 1
    return calc

# convert Cartesian to spherical coordinates. The sphere will have a minimum radius of 'radius'
def Cart2Spherical(radius=1.0, src=GetActiveSource()):
    calc=Calculator(src)
    strRad = str(radius)
    try:
        calc.Function = 'iHat*('+strRad+'+coordsZ)*cos(coordsY*'+strPi+'/180)*cos(coordsX*'+strPi+'/180) + jHat*('+strRad+'+coordsZ)*cos(coordsY*'+strPi+'/180)*sin(coordsX*'+strPi+'/180) + kHat*('+strRad+'+coordsZ)*sin(coordsY*'+strPi+'/180)'
    except:
        calc.Function = 'iHat*'+strRad+'*cos(coordsY*'+strPi+'/180)*cos(coordsX*'+strPi+'/180) + jHat*'+strRad+'*cos(coordsY*'+strPi+'/180)*sin(coordsX*'+strPi+'/180) + kHat*'+strRad+'*sin(coordsY*'+strPi+'/180)'
    calc.CoordinateResults = 1
    RenameSource('Cart2Spherical',calc)
    return calc

# adjust aspect ratio of Cartesian grid: multiplies ratios x coordinates
def GridAspectRatio(ratios, src=GetActiveSource()):
    calc=Calculator(src)
    try:
        calc.Function = 'iHat*'+str(ratios[0])+'*coordsX + jHat*'+str(ratios[1])+'*coordsY + kHat*'+str(ratios[2])+'*coordsZ'
    except:
        calc.Function = 'iHat*'+str(ratios[0])+'*coordsX + jHat*'+str(ratios[1])+'*coordsY'
    calc.CoordinateResults = 1
    return calc

# adjust aspect ratio of bounding box vector
def BoundAspectRatio(bounds, ratios, basis=1e3):
    Left   = bounds[0]*ratios[0]
    Right  = bounds[1]*ratios[0]
    Near   = bounds[2]*ratios[1]
    Far    = bounds[3]*ratios[1]
    if len(bounds) == 6 :
        Bottom = Pressure2Z(bounds[-2],ratios[2],basis)
        Top    = Pressure2Z(bounds[-1],ratios[2],basis)
        return Left,Right,Near,Far,Bottom,Top
    else:
        return Left,Right,Near,Far

# make filter selectable in pipeline browser, but don't show it
def MakeSelectable(src=GetActiveSource()):
    rep=Show(src)
    rep.Visibility=0


######### read in data, redefine pressure coordinates and change aspect ratio ###############
def loadData( fileName, outputDimensions=['pfull','lat','lon'], presCoords=1, aspectRatios=[1,1,1] ):
    # outputDimensions must be in same sequence as in netCDF file, except time (e.g. ['pfull','lat','lon'] )
    output_nc = NetCDFReader( FileName=[fileName] )

    if len(outputDimensions)>0 :
        outDims = '('+ outputDimensions[0]
        for dim in range(1,len(outputDimensions)):
            outDims = outDims + ', ' + outputDimensions[dim]
        outDims += ')'
        output_nc.Dimensions = outDims
    
    output_nc.SphericalCoordinates = 0
    output_nc.OutputType = 'Unstructured'
    output_nc.ReplaceFillValueWithNan = 0
    MakeSelectable()
    RenameSource(fileName,output_nc)
    
    CorrZ = CorrectZCoord(output_nc)
    MakeSelectable(CorrZ)
    RenameSource('CorrZ',CorrZ)
    
    if presCoords>0 :
        Coor = Pressure2Cart(CorrZ)
        RenameSource('LogP',Coor)
        MakeSelectable(Coor)
    
    try:
        if presCoords>0:
            AspRat = GridAspectRatio(aspectRatios, Coor)
        else:
            AspRat = GridAspectRatio(aspectRatios, CoorZ)
        RenameSource('AspectRatio',AspRat)
        MakeSelectable(AspRat)
    except:
        pass
    
    return output_nc,CorrZ,Coor,AspRat

######## some other usefull tools #################################################

# convert wind components from m/s to lat/timeStep, lon/timeStep, z/timeStep, and store it as vector W. Works with both pressure and height velocity, as long as vertAsp = [initial vertical range]/[present vertical range] is given
def CartWind2Atmos(src=GetActiveSource(), zonalComponentName='ucomp', meridionalComponentName='vcomp', secondsPerTimeStep=86400, verticalComponentName='none', vertAsp=1):
    W=Calculator(src)
    if verticalComponentName != 'none' :
        W.Function = '(' + \
        'iHat*'+zonalComponentName+'/(6.28*6.4e6*cos(coordsY*'+strPi+'/180))*360 +' + \
        'jHat*'+meridionalComponentName+'/('+strPi+'*6.4e6)*180 +' + \
        'kHat*'+verticalComponentName+'/'+vertAsp + \
        ')*'+str(secondsPerTimeStep) 
    else:
        W.Function =  '(' + \
        'iHat*'+zonalComponentName+'/(6.28*6.4e6*cos(coordsY*'+strPi+'/180))*360 +' + \
        'jHat*'+meridionalComponentName+'/('+strPi+'*6.4e6)*180' + \
        ')*'+str(secondsPerTimeStep)  
    W.ResultArrayName = 'W'
    RenameSource('CartWind2Atmos',W)
    MakeSelectable(W)
    # add the magnitdue of the wind vector, i.e wind strength. nice to have for color, threshold, glyph filters later on
    norm = Calculator(W)
    norm.Function = 'mag(W)'
    norm.ResultArrayName = 'normW'
    RenameSource('normW',norm)
    MakeSelectable(norm)
    # conversion invloves a division by zero over the poles. to avoid large numbers, cut away anything higher than 80 degrees
    clipS = Clip(norm)
    clipS.ClipType = 'Plane'
    clipS.ClipType.Normal = [0.0, 1.0, 0.0]
    clipS.ClipType.Origin  = [0.0, -80.0, 0.0]
    RenameSource('clipS',clipS)
    MakeSelectable(clipS)
    clipN = Clip(clipS)
    clipN.ClipType = 'Plane'
    clipN.ClipType.Normal = [0.0,-1.0, 0.0]
    clipN.ClipType.Origin  = [0.0, 80.0, 0.0]
    RenameSource('clipN',clipN)
    MakeSelectable(clipN)
    return W,norm,clipS,clipN




    
    

#!/usr/bin/python
# Filename: atmos_basic.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite CITATION HERE
#
# Python interface for ParaView (www.paraview.org). Reads netCDF file on a latitude - longitude and, if desired, pressure or height coordinates grid, including time evolution (if present). netCDF file needs to correspond to Climate and Forecast (FC) conventions (https://en.wikipedia.org/wiki/Climate_and_Forecast_Metadata_Conventions).
# Also Provides functions to modify Cartesian coordinates and wind components.

##### needed modules: paraview.simple, math #########################
from paraview.simple import *
from math import pi,log10

# some global constants
strPi = str(pi)[0:7]

##### define auxiliary functions ##################################
def CorrectZCoord(src=GetActiveSource()):
    """Make sure Z coordinate points in positive direction.

    Adds a Calculator filter to the pipeline, which takes the absolute value of coordsZ:
    src -- filter in pipeline to attach Calculator
    """
    calc=Calculator(src)
    calc.Function = 'iHat*coordsX + jHat*coordsY + kHat*abs(coordsZ)'
    calc.CoordinateResults = 1
    return calc

# define pressure to z coordinate conversion
def ConvertPressureString(pString, ratio=1.0, basis=1e3):
    """Convert Z coordinate conversion into a string for Calculator filter.

    Output is the string to be used inside the Calculator filter:
    pString -- value in hPa
    ratio   -- multiplicative factor for vertical coordinate
    basis   -- basis (surface) pressure to normalize
    """
    if ratio == 1.0:
        expression = 'log10(abs(' + pString + ')/' + str(basis) + ')'
    else:
        expression = 'log10(abs(' + pString + ')/' + str(basis) + ')*' + str(ratio)
    return expression

def Pressure2Z(plevel, ratio=1.0, basis=1e3):
    """Convert pressure to log-pressure (height).

    plevel -- the pressure level to convert
    ratio  -- multiplicative factor for vertical coordinate
    basis  -- basis (surface) pressure to normalize
    """
    level = -log10(plevel/basis)*ratio
    return level
 
def Pressure2Cart(src=GetActiveSource(), ratio=1.0, basis=1e3):
    """Convert pressure coordinates to Cartesian coordinates.

    Adds a Calculator filter to the pipeline
    src   -- filter in pipeline to attach to
    ratio -- multiplicative factor for vertical coordinate
    basis -- basis (surface) pressure to normalize
    """
    pFun = ConvertPressureString('coordsZ',ratio,basis)
    calc=Calculator(src)
    calc.Function = 'iHat*coordsX + jHat*coordsY - kHat*'+pFun
    calc.CoordinateResults = 1
    return calc

def Cart2Spherical(radius=1.0, src=GetActiveSource()):
    """Convert Cartesian to spherical coordinates. 

    Assumes X coordinate is longitude, Y coordinate latitude, Z coordinate vertical.
    Adds Calculator filter to the pipeline.
    radius -- radius of the sphere, where coordZ = basis
    src    -- filter in pipeline to attach to
    """
    calc=Calculator(src)
    strRad = str(radius)
    try:
        calc.Function = 'iHat*('+strRad+'+coordsZ)*cos(coordsY*'+strPi+'/180)*cos(coordsX*'+strPi+'/180) + jHat*('+strRad+'+coordsZ)*cos(coordsY*'+strPi+'/180)*sin(coordsX*'+strPi+'/180) + kHat*('+strRad+'+coordsZ)*sin(coordsY*'+strPi+'/180)'
    except:
        calc.Function = 'iHat*'+strRad+'*cos(coordsY*'+strPi+'/180)*cos(coordsX*'+strPi+'/180) + jHat*'+strRad+'*cos(coordsY*'+strPi+'/180)*sin(coordsX*'+strPi+'/180) + kHat*'+strRad+'*sin(coordsY*'+strPi+'/180)'
    calc.CoordinateResults = 1
    RenameSource('Cart2Spherical',calc)
    return calc

# 
def GridAspectRatio(ratios, src=GetActiveSource()):
    """Adjust aspect ratio of Cartesian grid: multiplies ratios x coordinates.

    Adds Calculator filter to the pipeline.
    ratios -- 2- or 3-vector with multiplicative factors for each spatial coordinate
    """
    calc=Calculator(src)
    try:
        calc.Function = 'iHat*'+str(ratios[0])+'*coordsX + jHat*'+str(ratios[1])+'*coordsY + kHat*'+str(ratios[2])+'*coordsZ'
    except:
        calc.Function = 'iHat*'+str(ratios[0])+'*coordsX + jHat*'+str(ratios[1])+'*coordsY'
    calc.CoordinateResults = 1
    return calc

# adjust aspect ratio of bounding box vector
def BoundAspectRatio(bounds, ratios, basis=1e3):
    """Adjust aspect ratio of bounding box (axes).

    Inputs are:
    bounds -- Physical bounds of 2D or 3D axes [Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]
    ratios -- Corrections to actually plotted axes
    basis  -- basis (surface) pressure to normalize
    Outputs are:
    Left,Right,Near,Far,Bottom,Top -- [Xmin,Xmax,Ymin,Ymax,Zmin,Zmax] of axes
    """
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

def MakeSelectable(src=GetActiveSource()):
    """Make filter selectable in pipeline browser, but don't show it."""
    rep=Show(src)
    rep.Visibility=0


######### read in data, redefine pressure coordinates and change aspect ratio ###############
def loadData( fileName, outputDimensions=['pfull','lat','lon'], presCoords=1, aspectRatios=[1,1,1], basis=1e3 ):
    """Load netCDF file, convert coordinates into useful aspect ratio.

    Adds file output_nc, Calculator CorrZ, Calculator LogP, and Calculator AspRat to the pipeline
    fileName         -- full path and file name of data to be read
    outputDimensions -- names of the dimensions within the netCDF file. Time should be excluded. Ordering matters!
    presCoords       -- whether (1) or not (0) Z coordinate should be logarithmic
    aspectRatios     -- how to scale coordinates [xscale,yscale,zscale]. Z coordinate is scaled after applying log10 if presCoords=1
    basis            -- basis (surface) pressure to normalize
    """ 
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
        Coor = Pressure2Cart(CorrZ,1,basis)
        RenameSource('LogP',Coor)
        MakeSelectable(Coor)
    else:
        Coor = []
    
    if presCoords>0:
        AspRat = GridAspectRatio(aspectRatios, Coor)
    else:
        AspRat = GridAspectRatio(aspectRatios, CorrZ)
    RenameSource('AspectRatio',AspRat)
    MakeSelectable(AspRat)
    
    return output_nc,CorrZ,Coor,AspRat

######## some other usefull tools #################################################

# 
def CartWind2Atmos(src=GetActiveSource(), zonalComponentName='ucomp', meridionalComponentName='vcomp', secondsPerTimeStep=86400, verticalComponentName='none', ratios=[1,1,1]):
    """Convert wind components from m/s to lat/timeStep, lon/timeStep, z/timeStep, and store it as vector W. 

    Works with both pressure and height velocity, as long as vertAsp = [initial vertical range]/[present vertical range] is given. 
    src                     -- filter in pipeline to attach to
    zonalComponentName      -- name of zonal wind component in pipeline
    meridionalComponentName -- name of meridional wind component in pipeline
    secondsPerTimeStep      -- duration of time step in seconds: 86400 for daily
    verticalComponentName   -- name of vertical component, or 'none'
    ratios                  -- Corrections to actually plotted axes
    Adds two Calculators to the pipeline:
    W     -- wind vector calculation
    normW -- magnitude of wind vector
    Adds two slices to the pipeline to remove division by zero close to poles:
    clipS -- remove south pole
    clipN -- remove north pole
    """
    W=Calculator(src)
    if verticalComponentName != 'none' :
        W.Function = '(' + \
        'iHat*'+zonalComponentName+'/(6.28*6.4e6*cos(coordsY/'+str(ratios[1])+'*'+strPi+'/180))*360 +' + \
        'jHat*'+meridionalComponentName+'/('+strPi+'*6.4e6)*180 +' + \
        'kHat*'+verticalComponentName+'/'+vertAsp + \
        ')*'+str(secondsPerTimeStep) 
    else:
        W.Function =  '(' + \
        'iHat*'+zonalComponentName+'/(6.28*6.4e6*cos(coordsY/'+str(ratios[1])+'*'+strPi+'/180))*360 +' + \
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
    clipS.ClipType.Origin  = [0.0, -80.0*ratios[1], 0.0]
    RenameSource('clipS',clipS)
    MakeSelectable(clipS)
    clipN = Clip(clipS)
    clipN.ClipType = 'Plane'
    clipN.ClipType.Normal = [0.0,-1.0, 0.0]
    clipN.ClipType.Origin  = [0.0, 80.0*ratios[1], 0.0]
    RenameSource('clipN',clipN)
    MakeSelectable(clipN)
    return W,norm,clipS,clipN

#
def DeleteAll():
    """Delete all objects in the pipeline browser."""
    for src in GetSources().values():
        Delete(src)
#
def HideAll():
    """Make all objects in pipeline browser invisible."""
    for src in GetSources().values():
    	Hide(src)
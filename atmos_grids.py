#!/usr/bin/python
# Filename: atmos_grids.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al
#
# Python interface for ParaView (www.paraview.org). Provides means to add axes in Cartesian or spherical coordinates. Needs atmos_basic.py

##### needed modules: atmos_basic ####################
from paraview.simple import *
##------------ General functions for adding/extracting and labeling planes ---------------------

# add a plane perpendicular to the Z direction
def AddZPlane(z, bounds=[0.0,360.0,-90.0,90.0], ratios=[1,1,1], logCoord=[2], basis=[1e3], data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    """Extract or add horizonal plane, i.e. perpendicular to Z direction.
    
    if data ==1, extracts a horizontal plane from data, else adds an empty plane for grid only
    
    z       -- physical value at which plane should be extracted/created
    bounds  -- extremities of grid box [Xmin,Xmax,Ymin,Ymax[,Zmin,Zmax]]
    ratios  -- aspect ratios of axes
    logCoord-- index/indices of coordinates to be treated logarithmically
    basis   -- basis to normalize logarithmic coordinate(s). If len==1, applied to all logCoord, otherwise must be same length as logCoord
    data    -- whether (1) or not (0) plane extracts data
    src     -- filter in pipeline to attach to
    AxisColor -- color of lines in RGB
    AxisWidth -- width of lines
    Adds one plane to the pipeline, either containing data or not
    """
    if 2 in logCoord :
    	zlevel = Lin2Log(z,ratios[2],basis[logCoord.index(2)])
    else:
    	zlevel = z*ratios[2]
    
    if len(bounds) == 4 :
        (Xmin,Xmax,Ymin,Ymax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    else:
        (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    if data == 1 : #extract plane from data
        Plane1 = Slice(src)
        Plane1.SliceType.Origin = [Xmin, Ymin, zlevel]
        Plane1.SliceType.Normal = [0.0, 0.0, 1.0]
        Plane1Rep = Show(Plane1)
        Plane1Rep.Representation = 'Surface'
    else: #create new plane otherwise
        Plane1 = Plane()
        Plane1.Origin = [Xmin  , Ymin, zlevel]
        Plane1.Point1 = [Xmax  , Ymin, zlevel] #Origin -> Point1 defines x-axis of plane
        Plane1.Point2 = [Xmin  , Ymax, zlevel] #Origin -> Point2 defines y-axis of plane
        Plane1Rep = Show(Plane1)
        Plane1Rep.Representation = 'Outline'
        Plane1Rep.AmbientColor = AxisColor
        Plane1Rep.LineWidth = AxisWidth
    return Plane1

def AddPresPlane(z, bounds=[0.0,360.0,-90.0,90.0], ratios=[1,1,1], logCoord=[2], basis=[1e3], data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    """AddPresPlane() is deprecated. Please use AddZPlane() instead"""
    warnings.warn("AddPresPlane() is deprecated. Please use AddZPlane() instead",DeprecationWarning)
    plane = AddZPlane(z, bounds, ratios, logCoord, basis, data, src, AxisColor, AxisWidth)
    return plane

# add a plane perpendicular to the Y direction
def AddYPlane(y,  bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], logCoord=[2], basis=[1e3], data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    """Add plane perpendicular to Y direction.
    
    if data ==1, extracts a plane from data, else adds an empty plane for grid only
    
    y       -- physical value at which plane should be extracted/created
    bounds  -- extremities of grid box [Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]
    ratios  -- aspect ratios of axes
    logCoord-- index/indices of coordinates to be treated logarithmically
    basis   -- basis to normalize logarithmic coordinate(s). If len==1, applied to all logCoord, otherwise must be same length as logCoord
    data    -- whether (1) or not (0) plane extracts data
    src     -- filter in pipeline to attach to
    AxisColor -- color of lines in RGB
    AxisWidth -- width of lines
    Adds one plane to the pipeline, either containing data or not
        """
    if 1 in logCoord :
    	yLoc = Lin2Log(y,ratios[1],basis[logCoord.index(1)])
    else:
    	yLoc = y*ratios[1]
    
    (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    if data == 1 : #extract plane from data
        Plane1 = Slice(src)
        Plane1.SliceType.Origin = [Xmin, yLoc, Zmin ]
        Plane1.SliceType.Normal = [0.0, 1.0, 0.0]
        Plane1Rep = Show(Plane1)
        Plane1Rep.Representation = 'Surface'
    else: #create new plane otherwise
        Plane1 = Plane()
        Plane1.Origin = [Xmin  , yLoc, Zmin ]
        Plane1.Point1 = [Xmax  , yLoc, Zmin ] #Origin -> Point1 defines x-axis of plane
        Plane1.Point2 = [Xmin  , yLoc, Zmax   ] #Origin -> Point2 defines y-axis of plane
        Plane1Rep = Show(Plane1)
        Plane1Rep.Representation = 'Outline'
        Plane1Rep.AmbientColor = AxisColor
        Plane1Rep.LineWidth = AxisWidth
    return Plane1

def AddLatPlane(lat,  bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], logCoord=[2], ratios=[1,1,1], basis=[1e3], data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    """AddLatPlane() is deprecated. Please use AddYPlane() instead."""
    warnings.warn("AddLatPlane() is deprecated. Please use AddYPlane() instead",DeprecationWarning)
    plane = AddYPlane(lat,  bounds, ratios, logCoord, basis, data, src, AxisColor, AxisWidth)
    return plane

# add a plane perpendicular to the X direction
def AddXPlane(x, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], logCoord=[2], basis=[1e3], data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    """Add plane perpendicular to the X direction.
    
    if data ==1, extracts a plane from data, else adds an empty plane for grid only
    
    x       -- physical value at which plane should be extracted/created
    bounds  -- extremities of grid box [Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]
    ratios  -- aspect ratios of axes
    logCoord-- index/indices of coordinates to be treated logarithmically
    basis   -- basis to normalize logarithmic coordinate(s). If len==1, applied to all logCoord, otherwise must be same length as logCoord
    data    -- whether (1) or not (0) plane extracts data
    src     -- filter in pipeline to attach to
    AxisColor -- color of lines in RGB
    AxisWidth -- width of lines
    Adds one plane to the pipeline, either containing data or not
    """
    if 0 in logCoord :
        xLoc = Lin2Log(x,ratios[0],basis[logCoord.index(0)])
    else:
        xLoc = x*ratios[0]

    (Xmin,Xmax ,Ymin,Ymax,Zmin ,Zmax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    if data == 1 : #extract plane from data
        Plane1 = Slice(src)
        Plane1.SliceType.Origin = [xLoc, Ymin, Zmin ]
        Plane1.SliceType.Normal = [1.0, 0.0, 0.0]
        Plane1Rep = Show(Plane1)
        Plane1Rep.Representation = 'Surface'
    else: #create new plane otherwise
        Plane1 = Plane()
        Plane1.Origin = [xLoc, Ymax , Zmin ]
        Plane1.Point1 = [xLoc, Ymin , Zmin ] #Origin -> Point1 defines x-axis of plane
        Plane1.Point2 = [xLoc, Ymax , Zmax   ] #Origin -> Point2 defines y-axis of plane
        Plane1Rep = Show(Plane1)
        Plane1Rep.Representation = 'Outline'
        Plane1Rep.AmbientColor = AxisColor
        Plane1Rep.LineWidth = AxisWidth
    return Plane1

def AddLonPlane(lon, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], basis=1e3, data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    """AddLonPlane() is deprecated. Please use AddXPlane() instead."""
    warnings.warn("AddLonPlane() is deprecated. Please use AddXPlane() instead",DeprecationWarning)
    plane = AddXPlane(lon, bounds, ratios, [2], [basis], data, src, AxisColor, AxisWidth)
    return plane

# add a label at a given Z coordinate, outside the domain
def AddZLabel(zlevel, LabelSize=5.0, bounds=[0.0,360,-90.0,90.0], ratios=[1,1,1], logCoord=[2], basis=[1e3], AxisColor=[0,0,0]):
    """Adds a label at a given Z coordinate, outside the domain (e.g. for axes labeling).

    Adds the label at a horizontal position relative to bounds, and at the correct height. 
    The label text is the same as the value of zlevel.
    zlevel    -- Z position and name of label
    LabelSize -- Size of the label text
    bounds    -- bounds of axes: label is positioned just outside these bounds [Xmin,Xmax,Ymin,Ymax[,Zmin,Zmax]]
    ratios    -- aspect ratios of axes
    logCoord  -- index/indices of coordinates to be treated logarithmically
    basis     -- basis to normalize logarithmic coordinate(s). If len==1, applied to all logCoord, otherwise must be same length as logCoord
    AxisColor -- color of lines in RGB
    Adds the Label itself and 4 transforms (one at each corner of the axes) to the pipeline
    """
    if 2 in logCoord :
        Z = Lin2Log(zlevel,ratios[2],basis[logCoord.index(2)])
    else:
        Z = zlevel*ratios[2]

    if len(bounds) == 4 :
        (Xmin,Xmax ,Ymin,Ymax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    else:
        (Xmin,Xmax ,Ymin,Ymax,Zmin ,Zmax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    LabelScale = [abs(LabelSize), abs(LabelSize), abs(LabelSize)]
    Label = a3DText(Text=str(zlevel))
    
    percentOff = 0.02
    Trans1 = Transform()
    Trans1.Transform.Translate = [ Xmax , Ymax+percentOff*(Ymax-Ymin), Z ]
    Trans1.Transform.Rotate = [ 90.0,90.0, 0.0 ]
    Trans1.Transform.Scale = LabelScale
    Trans1.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface' #hidden from behind
    Trans2 = Transform(Label)
    Trans2.Transform.Translate = [ Xmax +percentOff*(Xmax -Xmin), Ymin, Z ]
    Trans2.Transform.Rotate = [ 90.0, 0.0, 0.0 ]
    Trans2.Transform.Scale = LabelScale
    Trans2.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    Trans3 = Transform(Label)
    Trans3.Transform.Translate = [Xmin-percentOff*(Xmax -Xmin), Ymax, Z ]
    Trans3.Transform.Rotate = [ 90.0,180.0, 0.0 ]
    Trans3.Transform.Scale = LabelScale
    Trans3.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    Trans4 = Transform(Label)
    Trans4.Transform.Translate = [ Xmin, Ymin-percentOff*(Ymax-Ymin), Z ]
    Trans4.Transform.Rotate = [ 90.0,-90.0, 0.0 ]
    Trans4.Transform.Scale = LabelScale
    Trans4.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    return Label,Trans1,Trans2,Trans3,Trans4


def AddPresLabel(plevel, LabelSize=5.0, bounds=[0.0,360,-90.0,90.0], ratios=[1,1,1], basis=1e3, AxisColor=[0,0,0]):
    """AddPresLabel() is deprecated. Please use AddZLabel() instead."""
    warnings.warn("AddPresLabel() is deprecated. Please use AddZLabel() instead.",DeprecationWarning)
    (label,trans1,trans2,trans3,trans4)=AddZLabel(plevel, LabelSize, bounds, ratios, [2], [basis], AxisColor)
    return label,trans1,trans2,trans3,trans4

# add label at a given Y coordinate outside the domain
def AddYLabel(ylevel, LabelSize=5.0, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], logCoord=[2], basis=[1e3], AxisColor=[0,0,0]):
    """Adds a label at a given Y position.

    Adds the label at a horizontal and vertical position relative to bounds. 
    
    The label text is the same as the value of ylevel.
    ylevel    -- Y position and name of label
    LabelSize -- Size of the label text
    bounds    -- bounds of axes: label is positioned just outside these bounds [Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]
    ratios    -- aspect ratios of axes
    logCoord  -- index/indices of coordinates to be treated logarithmically
    basis     -- basis to normalize logarithmic coordinate(s). If len==1, applied to all logCoord, otherwise must be same length as logCoord
    AxisColor -- color of lines in RGB
    Adds the Label itself and 2 transforms (one on each side of the axes) to the pipeline
        """
    if 1 in logCoord :
        yLoc = Lin2Log(ylevel,ratios[1],basis[logCoord.index(1)])
    else:
        yLoc = ylevel*ratios[1]

    (Xmin,Xmax ,Ymin,Ymax,Zmin ,Zmax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    Z = Zmin  - LabelSize*1.5
    LabelScale = [LabelSize, LabelSize, LabelSize]
    Label = a3DText(Text=str(ylevel))
    if ylevel < 0.0:
        yOffset = LabelSize*1.5
    elif ylevel > 0.0:
        yOffset = LabelSize
    else:
        yOffset = LabelSize*0.5
    Trans1 = Transform()
    Trans1.Transform.Translate = [ Xmax , yLoc-yOffset, Z ]
    Trans1.Transform.Rotate = [ 90.0, 90.0, 0.0 ]
    Trans1.Transform.Scale = LabelScale
    Trans1.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface' #hidden from behind
    Trans2 = Transform(Label)
    Trans2.Transform.Translate = [ Xmin, yLoc+yOffset, Z ]
    Trans2.Transform.Rotate = [ 90.0,-90.0, 0.0 ]
    Trans2.Transform.Scale = LabelScale
    Trans2.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    return Label,Trans1,Trans2

def AddLatLabel(lat, LabelSize=5.0, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], basis=1e3, AxisColor=[0,0,0]):
    """AddLatLabel() is deprecated. Please use AddYLabel()."""
    warnings.warn("AddLatLabel() is deprecated. Please use AddYLabel() instead.",DeprecationWarning)
    (label,trans1,trans2) = AddYLabel(lat,LabelSize,bounds,ratios,[2],[basis],AxisColor)
    return label,trans1,trans2

# add label at a given X coordinate outside the domain
def AddXLabel(xlevel, LabelSize=5.0, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], logCoord=[2], basis=[1e3], AxisColor=[0,0,0]):
    """Adds a label at a given longitude.

    Adds the label at a horizontal and vertical position relative to bounds. 
    
    The label text is the same as the value of xlevel.
    xlevel    -- X position and name of label
    LabelSize -- Size of the label text
    bounds    -- bounds of axes: label is positioned just outside these bounds [Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]
    ratios    -- aspect ratios of axes
    logCoord  -- index/indices of coordinates to be treated logarithmically
    basis     -- basis to normalize logarithmic coordinate(s). If len==1, applied to all logCoord, otherwise must be same length as logCoord
    AxisColor -- color of lines in RGB
    Adds the Label itself and 2 transforms (one on each side of the axes) to the pipeline
        """
    if 0 in logCoord :
        xLoc = Lin2Log(xlevel,ratios[0],basis[logCoord.index(0)])
    else:
        xLoc = xlevel*ratios[0]

    (Xmin,Xmax ,Ymin,Ymax,Zmin ,Zmax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    Z = Zmin  - LabelSize*1.5
    LabelScale = [LabelSize, LabelSize, LabelSize]
    Label = a3DText(Text=str(xlevel))

    if xlevel < 100.0:
        xOffset = LabelSize*1.0
    else:
        xOffset = LabelSize*1.5
    Trans1 = Transform()
    Trans1.Transform.Translate = [xLoc-xOffset, Ymin, Z ]
    Trans1.Transform.Rotate = [ 90.0, 0.0, 0.0 ]
    Trans1.Transform.Scale = LabelScale
    Trans1.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface' #hidden from behind
    Trans2 = Transform(Label)
    Trans2.Transform.Translate = [xLoc+xOffset, Ymax, Z ]
    Trans2.Transform.Rotate = [ 90.0,180.0, 0.0 ]
    Trans2.Transform.Scale = LabelScale
    Trans2.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    return Label,Trans1,Trans2

def AddLonLabel(lon, LabelSize=5.0, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], basis=1e3, AxisColor=[0,0,0]):
    """AddLonLabel() is deprecated. Please use AddXlabel() instead."""
    warnings.warn("AddLonLabel() is deprecated. Please use AddXLabel() instead.",DeprecationWarning)
    (label,trans1,trans2) = AddXLabel(lon,LabelSize,bounds,ratios,[2],[basis],AxisColor)
    return label,trans1,trans2

# move generic label with given position and rotation
def AddAxisLabel(Labl,Trans,Rot, AxisColor=[0,0,0], LabelSize=5.0,):
    """Move a generic label.

    Modifies Labl according to the transformations Trans, Rot, LabelSize:
    Labl  -- input object, usually a3DText
    Trans -- translation vector [dx,dy,dz]
    Rot   -- rotation vector [rx,ry,rz]
    LabelSize -- Size of the label text
    AxisColor -- text color
    Adds a Transform filter to the pipeline
    """
    LabelScale = [LabelSize, LabelSize, LabelSize]
    TransLab = Transform(Labl)
    TransLab.Transform.Translate = Trans
    TransLab.Transform.Rotate = Rot
    TransLab.Transform.Scale = LabelScale
    TransLab.SMProxy.InvokeEvent('UserEvent','HideWidget')
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    return TransLab

# Add vertical labels in spherical coordinates
def SphericalLabels(radius=1, ratios=[1,1,1], logCoord=[2], basis=[1e3], shellValues=[100,10,1], labelPosition=[170, 10], labelSize=1.0):
    """Label pressure surface(s) that has(have) been extracted in spherical geometry (shell).

    radius   -- radius of sphere with r(basis)=radius
    ratios   -- aspect ratios of axes before converting to sphere
    logCoord -- index/indices of coordinates to be treated logarithmically
    basis    -- basis to normalize logarithmic coordinate(s). If len==1, applied to all logCoord, otherwise must be same length as logCoord
    shellValues -- vector of values to be labelled, in original units (before ratios conversion)
    labelPosition -- the position in [longitude,latitude] on the sphere
    labelSize -- multiplicative factor for size of the labels
    Adds the label text, a Transform and a Calculator filter to the pipeline
        """
    if 2 in logCoord :
        if len(basis) > 0:
            bas = basis[logCoord.index(2)]
        else:
            bas = basis[0]
    for p in range(len(shellValues)):
        ps = shellValues[p]
        if 2 in logCoord :
            if ps >= bas :
                fac = 1.01
            else:
                fac = 0.99
            labelRadius = radius + Lin2Log(fac*ps,ratios[2],bas)
        else:
            labelRadius = radius + 1.01*ps*ratios[2]
        txt=a3DText()
        txt.Text = str(abs(ps))
        MakeSelectable(txt)
        RenameSource('Text'+str(ps)+'[Z]',txt)
        Transform1=Transform(txt)
        Transform1.Transform="Transform"
        Transform1.Transform.Rotate    =[0.0,0.0,90.0]
        if abs(ps) >= 1:
            Transform1.Transform.Scale     =[3.0*labelSize,3.0*labelSize, 1.0*labelSize]
        else:
            Transform1.Transform.Scale     =[2.0*labelSize,2.0*labelSize, 1.0*labelSize]
        Transform1.Transform.Translate =[labelPosition[0], labelPosition[1], labelRadius]
        MakeSelectable(Transform1)
        Text2Sphere = Cart2Spherical(0,Transform1)
        Text_disp=Show()
        Text_disp.DiffuseColor = [0.0, 0.0, 0.0]
        Text_disp.Opacity=0.5


def AtmosLabels(radius=1, ratios=[1,1,1], basis=1e3, shellValues=[100,10,1], labelPosition=[170, 10]):
    """AtmosLabels() is deprecated. Please use SphericalLabels() instead."""
    SphericalLabels(radius,ratios,[2],[basis],shellValues,labelPosition)

# add a water mark to denote authorship of plot
def WaterMark(waterMark, markRadius=1, markPosition=[250, 10], markSize=1.0):
    """Add a water mark on the sphere as signature.

    waterMark  -- text to project onto the sphere
    markRadius -- radius of water mark, already converted from physical units
    markPosition -- position of water mark in [longitude, latitude] on the sphere
    markSize   -- size of the water mark label
    Adds a a3DText, Transform, and a Calculator filter to the pipeline
    """
    txt=a3DText()
    txt.Text = waterMark
    MakeSelectable()
    RenameSource('WaterMark',txt)
    Transform2=Transform()
    Transform2.Transform="Transform"
    Transform2.Transform.Scale     =[2.0*markSize,2.0*markSize, 1.0*markSize]
    Transform2.Transform.Translate =[markPosition[0], markPosition[1], markRadius]
    MakeSelectable(Transform2)
    Mark2Sphere = Cart2Spherical(0,Transform2)
    Text_disp=Show()
    Text_disp.DiffuseColor = [0.0, 0.0, 0.0]
    Text_disp.Opacity=0.1

######## combine some of the above to create a suite of atmospheric shells. #########################
########  this replaces the grid in spherical geometry                      #########################
def SphericalShells(radius=1, ratios=[1,1,1], logCoord=[2], basis=[1e3], src=GetActiveSource(), shellValues=[100,10,1], labels=1, labelPosition=[170, 10], waterMark='none', markPosition=[250, 10], labelSize=1.0):
    """Add spherical shells as grid to spherical geometry, or to visualize specific pressure levels.

    Adds as many shells as there are shellValues.
    
    Combines AddZPlane, SphericalLabels, and WaterMark to create a grid in spherical coordinates.
    radius        -- radius of sphere with r(basis)=radius
    ratios        -- aspect ratios of axes before converting to sphere
    logCoord      -- index/indices of coordinates to be treated logarithmically
    basis         -- basis to normalize logarithmic coordinate(s). If len==1, applied to all logCoord, otherwise must be same length as logCoord
    src           -- filter in pipeline to attach to
    shellValues   -- vector of values to be labelled, in original units (hPa,km)
    labels        -- add labels (>0) or not (0)
    labelPosition -- the position in [longitude,latitude] on the sphere
    waterMark     -- string with text for water mark, or 'none'
    markPosition  -- position of water mark in [longitude, latitude] on the sphere
    labelSize     -- multiplicative factor for size of labels and water mark
    Adds a Plane, two Calculators, a a3DText, and a Transform to the pipeline.
    If a water mark is included, adds an additional a3DText, Transform, and a Calculator filter to the pipeline.
    Returns a list of all pressure planes for further treatment.
    """
    Planes=[]
    for ps in shellValues:
        TropoSlice = AddZPlane(ps, ratios=ratios, logCoord=logCoord, basis=basis, data=1, src=src)
        MakeSelectable(TropoSlice)
        RenameSource(str(ps)+'[Z]',TropoSlice)
        Cart2Sphere = Cart2Spherical(radius,TropoSlice)
        TropoSlice_disp=Show()
        TropoSlice_disp.Opacity=0.1
        Planes.append(Cart2Sphere)

        if labels>0:
            SphericalLabels(radius, ratios, logCoord, basis, [ps], labelPosition, labelSize)
    
    # add watermark
    if waterMark != 'none':
        if 2 in logCoord :
            if len(basis) > 0:
                bas = basis[logCoord.index(2)]
            else:
                bas = basis[0]
            labelRadius = Lin2Log(min(shellValues),ratios[2],bas)
        else:
            labelRadius = radius + shellValues[-1]*ratios[2]
        WaterMark(waterMark, labelRadius, markPosition, labelSize)
    return Planes

def AtmosShells(radius=1, ratios=[1,1,1], basis=1e3, src=GetActiveSource(), shellValues=[100,10,1], labels=1, labelPosition=[170, 10], waterMark='none', markPosition=[250, 10]):
    """AtmosShells() is deprecated. Please use SphericalShells() instead."""
    SphericalShells(radius,ratios,[2],[basis],src,shellValues,labels,labelPosition,waterMark,markPosition)

###### add full set of grids and lables in rectangular geometry ############################
def AddGrid(xlevels=[0,90,180,270], ylevels=[-60,-30,0,30,60], zlevels=[100,10,1,0.1], bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], logCoord=[2], basis=[1e3], AxisNames=["lon","lat","pressure [hPa]"], AxisColor=[0,0,0], AxisWidth=2.0,LabelSize=5.0):
    """Add a full grid with grid lines and labels.

Adds as many X,Y,Z grid lines as needed.
    xlevels   -- vector with X grid positions
    ylevels   -- vector with Y grid positions
    zlevels   -- vector with Z grid levels
    bounds    -- grid bounds 
    ratios    -- axes ratios
    basis     -- basis (surface) pressure
    AxisNames -- names of x,y,z axis
    AxisColor -- color of lines in RGB
    AxisWidth -- width of lines
    LabelSize -- Size of the label text
    Remark: This adds a lot of objects and filters to the pipeline, and should probably only be used once the visualization itself is finished. This function can be called even if there is no data loaded, or with data that has nothing to do with geophysical fluids.
        """
    (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax) = BoundAspectRatio(bounds, ratios, logCoord, basis)
    absLsz = abs(LabelSize)
    #Z grid
    if len(zlevels) > 0:
        for z in zlevels:
            PlaneTmp = AddZPlane(z, bounds[0:4], ratios, logCoord, basis, AxisColor=AxisColor, AxisWidth=AxisWidth)
            RenameSource("Plane"+str(z),PlaneTmp)
            if abs(LabelSize) > 0.0 :
                (label,transa,transb,transc,transd) = AddZLabel(z, absLsz, bounds[0:4], ratios, logCoord, basis, AxisColor)
                RenameSource("Label"+str(z),label)
                Show(transa)
                Show(transb)
                Show(transc)
                Show(transd)
        
        #Z axis label
        if LabelSize > 0.0 :
            BoxH = Zmax - Zmin 
            LabelTmp = a3DText(Text=AxisNames[2])
            RenameSource("ZLabel",LabelTmp)
            LabelPushFac = len(str(max(max(zlevels),abs(min(zlevels)))))+2
            if max(zlevels) < abs(min(zlevels)):
                LabelPushFac += 2
            Transx = [
                      Xmax , Xmax +LabelPushFac*absLsz, Xmin, Xmin-LabelPushFac*absLsz ]
            Transy = [Ymax+LabelPushFac*absLsz, Ymin, Ymin-LabelPushFac*absLsz, Ymax ]
            Rotx = [ 180.0,   0.0, 180.0, 180.0 ]
            Roty = [  90.0, -90.0,  90.0,  90.0 ]
            Rotz = [   0.0,  90.0, 180.0,  90.0 ]
            for i in range(len(Transx)):
                Trans = [ Transx[i], Transy[i], Zmin + BoxH*0.5 -4.0*absLsz ]
                Rot = [ Rotx[i], Roty[i], Rotz[i] ]
                TransPres = AddAxisLabel(LabelTmp, Trans, Rot, AxisColor, absLsz)
                
    # for other coordinate labels
    Z = Zmin  - absLsz*3.0
    #Y grid
    if len(ylevels) > 0:
        for y in ylevels:
            PlaneTmp = AddYPlane(y, bounds, ratios, logCoord, basis, AxisColor=AxisColor, AxisWidth=AxisWidth)
            RenameSource("Plane"+str(y),PlaneTmp)
            if abs(LabelSize) > 0.0 :
                (label,transa,transb) = AddYLabel(y, absLsz, bounds, ratios, logCoord, basis, AxisColor)
                RenameSource("Label"+str(y),label)
                Show(transa)
                Show(transb)
        #Y axis label
        if LabelSize > 0.0 :
            LabelTmp = a3DText(Text=AxisNames[1])
            RenameSource("YLabel",LabelTmp)
            midY = 0.5*(Ymax+Ymin)
            Trans = [Xmax , midY-2.5*absLsz, Z]
            Rot = [ 90.0, 90.0, 0.0 ]
            TransLat = AddAxisLabel(LabelTmp, Trans, Rot, AxisColor, absLsz)
            
            Trans = [Xmin, midY+2.5*absLsz, Z]
            Rot = [ 90.0,-90.0, 0.0 ]
            TransLat = AddAxisLabel(LabelTmp, Trans, Rot, AxisColor, absLsz)

    #X grid
    if len(xlevels) > 0 :
        for x in xlevels:
            PlaneTmp = AddXPlane(x, bounds, ratios, logCoord, basis, AxisColor=AxisColor, AxisWidth=AxisWidth)
            RenameSource("Plane"+str(x)+"[X]",PlaneTmp)
            if abs(LabelSize) > 0.0 :
                (label,transa,transb) = AddXLabel(x, absLsz, bounds, ratios, logCoord, basis, AxisColor)
                Show(transa)
                Show(transb)

        #X axis label
        if LabelSize > 0.0 :
            LabelTmp = a3DText(Text=AxisNames[0])
            RenameSource("XLabel",LabelTmp)
            Trans = [0.5*(Xmin+Xmax )-3.0*absLsz, Ymin, Z]
            Rot = [ 90.0,   0.0, 0.0 ]
            TransLon = AddAxisLabel(LabelTmp, Trans, Rot, AxisColor, absLsz)
            
            Trans = [0.5*(Xmin+Xmax)+3.0*absLsz, Ymax, Z]
            Rot = [ 90.0, 180.0, 0.0 ]
            TransLon = AddAxisLabel(LabelTmp, Trans, Rot, AxisColor, absLsz)






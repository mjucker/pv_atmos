#!/usr/bin/python
# Filename: atmos_grids.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite CITATION HERE
#
# Python interface for ParaView (www.paraview.org). Provides means to add axes in Cartesian or spherical coordinates. Needs atmos_basic.py

##### needed modules: paraview.simple, math, atmos_basic ####################


##------------ General functions for adding/extracting and labeling planes ---------------------

#extract or add horizonal plane, i.e. plane at given pressure
def AddPresPlane(pres, bounds=[0.0,360.0,-90.0,90.0], ratios=[1,1,1], basis=1e3, data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    plevel = Pressure2Z(pres,ratios[2],basis)
    if len(bounds) == 4 :
        (Left,Right,Near,Far) = BoundAspectRatio(bounds,ratios,basis)
    else:
        (Left,Right,Near,Far,Bottom,Top) = BoundAspectRatio(bounds,ratios,basis)
    if data == 1 : #extract plane from data
        Plane1 = Slice(src)
        Plane1.SliceType.Origin = [Left, Near, plevel]
        Plane1.SliceType.Normal = [0.0, 0.0, 1.0]
        Plane1Rep = GetDisplayProperties(Plane1)
        Plane1Rep.Representation = 'Surface'
    else: #create new plane otherwise
        #level = Pressure2Z(plevel,ratio, basis)
        Plane1 = Plane()
        Plane1.Origin = [Left  , Near, plevel]
        Plane1.Point1 = [Right , Near, plevel] #Origin -> Point1 defines x-axis of plane
        Plane1.Point2 = [Left  , Far , plevel] #Origin -> Point2 defines y-axis of plane
        Plane1Rep = GetDisplayProperties(Plane1)
        Plane1Rep.Representation = 'Outline'
        Plane1Rep.AmbientColor = AxisColor
        Plane1Rep.LineWidth = AxisWidth
    return Plane1

# plane parallel to longitude, at given latitude
def AddLatPlane(lat,  bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], basis=1e3, data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    (Left,Right,Near,Far,Bottom,Top) = BoundAspectRatio(bounds,ratios,basis)
    if data == 1 : #extract plane from data
        Plane1 = Slice(src)
        Plane1.SliceType.Origin = [Left, lat, Bottom]
        Plane1.SliceType.Normal = [0.0, 1.0, 0.0]
        Plane1Rep = GetDisplayProperties(Plane1)
        Plane1Rep.Representation = 'Surface'
    else: #create new plane otherwise
        Plane1 = Plane()
        Plane1.Origin = [Left  , lat, Bottom]
        Plane1.Point1 = [Right , lat, Bottom] #Origin -> Point1 defines x-axis of plane
        Plane1.Point2 = [Left  , lat, Top   ] #Origin -> Point2 defines y-axis of plane
        Plane1Rep = GetDisplayProperties(Plane1)
        Plane1Rep.Representation = 'Outline'
        Plane1Rep.AmbientColor = AxisColor
        Plane1Rep.LineWidth = AxisWidth
    return Plane1

# plane parallel to latitude, at given longitude
def AddLonPlane(lon, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], basis=1e3, data=0, src=GetActiveSource(), AxisColor=[0,0,0], AxisWidth=1.0):
    (Left,Right,Near,Far,Bottom,Top) = BoundAspectRatio(bounds,ratios,basis)
    if data == 1 : #extract plane from data
        Plane1 = Slice(src)
        Plane1.SliceType.Origin = [lon, Near, Bottom]
        Plane1.SliceType.Normal = [1.0, 0.0, 0.0]
        Plane1Rep = GetDisplayProperties(Plane1)
        Plane1Rep.Representation = 'Surface'
    else: #create new plane otherwise
        Plane1 = Plane()
        Plane1.Origin = [lon, Far , Bottom]
        Plane1.Point1 = [lon, Near, Bottom] #Origin -> Point1 defines x-axis of plane
        Plane1.Point2 = [lon, Far , Top   ] #Origin -> Point2 defines y-axis of plane
        Plane1Rep = GetDisplayProperties(Plane1)
        Plane1Rep.Representation = 'Outline'
        Plane1Rep.AmbientColor = AxisColor
        Plane1Rep.LineWidth = AxisWidth
    return Plane1# add label of pressure surface in plane geometry
def AddPresLabel(plevel, LabelSize=5.0, bounds=[0.0,360,-90.0,90.0], ratios=[1,1,1], basis=1e3, AxisColor=[0,0,0]):
    Z = Pressure2Z(plevel,ratios[2],basis)
    if len(bounds) == 4 :
        (Left,Right,Near,Far) = BoundAspectRatio(bounds,ratios)
    else:
        (Left,Right,Near,Far,Bottom,Top) = BoundAspectRatio(bounds,ratios,basis)
    LabelScale = [abs(LabelSize), abs(LabelSize), abs(LabelSize)]
    Label = a3DText(Text=str(plevel))
    
    Trans1 = Transform()
    Trans1.Transform.Translate = [ Right, Far+1, Z ]
    Trans1.Transform.Rotate = [ 90.0,90.0, 0.0 ]
    Trans1.Transform.Scale = LabelScale
    Trans1.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface' #hidden from begind
    Trans2 = Transform(Label)
    Trans2.Transform.Translate = [ Right+1.0, Near, Z ]
    Trans2.Transform.Rotate = [ 90.0, 0.0, 0.0 ]
    Trans2.Transform.Scale = LabelScale
    Trans2.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    Trans3 = Transform(Label)
    Trans3.Transform.Translate = [Left-1.0, Far, Z ]
    Trans3.Transform.Rotate = [ 90.0,180.0, 0.0 ]
    Trans3.Transform.Scale = LabelScale
    Trans3.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    Trans4 = Transform(Label)
    Trans4.Transform.Translate = [ Right, Near-1.0, Z ]
    Trans4.Transform.Rotate = [ 90.0,-90.0, 0.0 ]
    Trans4.Transform.Scale = LabelScale
    Trans4.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    return Label,Trans1,Trans2,Trans3,Trans4

# add label of latitude surface in plane geometry
def AddLatLabel(lat, LabelSize=5.0, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], basis=1e3, AxisColor=[0,0,0]):
    (Left,Right,Near,Far,Bottom,Top) = BoundAspectRatio(bounds,ratios,basis)
    Z = Bottom - LabelSize*1.5
    LabelScale = [LabelSize, LabelSize, LabelSize]
    Label = a3DText(Text=str(lat))
    if lat < 0.0:
        latOffset = LabelSize*1.5
    elif lat > 0.0:
        latOffset = LabelSize
    else:
        latOffset = LabelSize*0.5
    Trans1 = Transform()
    Trans1.Transform.Translate = [ Right, lat-latOffset, Z ]
    Trans1.Transform.Rotate = [ 90.0, 90.0, 0.0 ]
    Trans1.Transform.Scale = LabelScale
    Trans1.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface' #hidden from begind
    Trans2 = Transform(Label)
    Trans2.Transform.Translate = [ Left, lat+latOffset, Z ]
    Trans2.Transform.Rotate = [ 90.0,-90.0, 0.0 ]
    Trans2.Transform.Scale = LabelScale
    Trans2.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    return Label,Trans1,Trans2

# add label for longitude surface in plane geometry
def AddLonLabel(lon, LabelSize=5.0, bounds=[0.0,360.0,-90.0,90.0,1e3,0.1], ratios=[1,1,1], basis=1e3, AxisColor=[0,0,0]):
    (Left,Right,Near,Far,Bottom,Top) = BoundAspectRatio(bounds,ratios,basis)
    Z = Bottom - LabelSize*1.5
    LabelScale = [LabelSize, LabelSize, LabelSize]
    Label = a3DText(Text=str(lon))
    if lon < 100.0:
        lonOffset = LabelSize*1.0
    else:
        lonOffset = LabelSize*1.5
    Trans1 = Transform()
    Trans1.Transform.Translate = [lon-lonOffset, Near, Z ]
    Trans1.Transform.Rotate = [ 90.0, 0.0, 0.0 ]
    Trans1.Transform.Scale = LabelScale
    Trans1.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface' #hidden from begind
    Trans2 = Transform(Label)
    Trans2.Transform.Translate = [lon+lonOffset, Far, Z ]
    Trans2.Transform.Rotate = [ 90.0,180.0, 0.0 ]
    Trans2.Transform.Scale = LabelScale
    Trans2.SMProxy.InvokeEvent('UserEvent','HideWidget') #don't show box
    Rep = GetDisplayProperties()
    Rep.Representation = 'Surface'
    Rep.DiffuseColor = AxisColor
    Rep.BackfaceRepresentation = 'Cull Backface'
    return Label,Trans1,Trans2

# add generic axis label
def AddAxisLabel(Labl,Trans,Rot, LabelSize=5.0, AxisColor=[0,0,0]):
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


# label pressure surface that has been extracted in spherical geometry
def AtmosLabels(radius=1, ratios=[1,1,1], basis=1e3, shellValues=[100,10,1], labelPosition=[170, 10]):
    for p in range(len(shellValues)):
        ps = shellValues[p]
        labelRadius = radius+Pressure2Z(.99*ps,ratios[2],basis)
        txt=a3DText()
        txt.Text = str(abs(ps))
        MakeSelectable(txt)
        RenameSource('Text'+str(ps)+'hPa',txt)
        Transform1=Transform(txt)
        Transform1.Transform="Transform"
        Transform1.Transform.Rotate    =[0.0,0.0,90.0]
        if abs(ps) >= 1:
            Transform1.Transform.Scale     =[3.0,3.0, 1.0]
        else:
            Transform1.Transform.Scale     =[2.0,2.0, 1.0]
        Transform1.Transform.Translate =[labelPosition[0], labelPosition[1], labelRadius]
        MakeSelectable(Transform1)
        Text2Sphere = Cart2Spherical(0,Transform1)
        Text_disp=Show()
        Text_disp.DiffuseColor = [0.0, 0.0, 0.0]
        Text_disp.Opacity=0.5

# add a water mark to denote authorship of plot
def WaterMark(waterMark, markRadius=1, markPosition=[250, 10]):
    txt=a3DText()
    txt.Text = waterMark
    MakeSelectable()
    RenameSource('WaterMark',txt)
    Transform2=Transform()
    Transform2.Transform="Transform"
    Transform2.Transform.Scale     =[2.0,2.0, 1.0]
    Transform2.Transform.Translate =[markPosition[0], markPosition[1], markRadius]
    MakeSelectable(Transform2)
    Mark2Sphere = Cart2Spherical(0,Transform2)
    Text_disp=Show()
    Text_disp.DiffuseColor = [0.0, 0.0, 0.0]
    Text_disp.Opacity=0.5

######## combine some of the above to create a suite of atmospheric shells. #########################
########  this replaces the grid in spherical geometry                      #########################
def AtmosShells(radius=1, ratios=[1,1,1], basis=1e3, src=GetActiveSource(), shellValues=[100,10,1], labels=1, labelPosition=[170, 10], waterMark='none', markPosition=[250, 10]):
    Planes=[]
    for ps in shellValues:
        TropoSlice = AddPresPlane(ps, ratios=ratios, basis=basis, data=1, src=src)
        MakeSelectable(TropoSlice)
        RenameSource(str(ps)+'hPa',TropoSlice)
        Cart2Sphere = Cart2Spherical(radius,TropoSlice)
        TropoSlice_disp=Show()
        TropoSlice_disp.Opacity=0.1
        Planes.append(Cart2Sphere)

        if labels>0:
            AtmosLabels(radius,ratios,basis,[ps], labelPosition)
    
    # add watermark
    if waterMark != 'none':
        labelRadius = Pressure2Z(.99*min(shellValues),ratios[2],basis)
        WaterMark(waterMark, labelRadius, markPosition)
    return Planes

###### add full set of grids and lables in rectangular geometry ############################

def AddGrid(press=[100,10,1,0.1], lats=[-60,-30,0,30,60], lons=[0,90,180,270], bounds=[0.0,360.0,-90.0,90.0,1000.0,0.1], ratios=[1,1,1], basis=1e3, AxisColor=[0,0,0], AxisWidth=2.0,LabelSize=5.0):
    
    (Left,Right,Near,Far,Bottom,Top) = BoundAspectRatio(bounds,ratios,basis)
    absLsz = abs(LabelSize)
    #pressure grid
    if len(press) > 0:
        for pres in press:
            PlaneTmp = AddPresPlane(pres, bounds[0:4], ratios, basis, AxisColor=AxisColor, AxisWidth=AxisWidth)
            RenameSource("Plane"+str(pres),PlaneTmp)
            if abs(LabelSize) > 0.0 :
                (label,transa,transb,transc,transd) = AddPresLabel(pres, absLsz, bounds[0:4], ratios, basis, AxisColor)
                RenameSource("Label"+str(pres),label)
                Show(transa)
                Show(transb)
                Show(transc)
                Show(transd)
        
        #pressure label
        if LabelSize > 0.0 :
            BoxH = Top - Bottom
            LabelTmp = a3DText(Text='pressure')
            RenameSource("PresLabel",LabelTmp)
            LabelPushFac = 3.7
            Transx = [
                      Right, Right+LabelPushFac*absLsz, Left, Left-LabelPushFac*absLsz ]
            Transy = [Far+LabelPushFac*absLsz, Near, Near-LabelPushFac*absLsz, Far ]
            Rotx = [ 180.0,   0.0, 180.0, 180.0 ]
            Roty = [  90.0, -90.0,  90.0,  90.0 ]
            Rotz = [   0.0,  90.0, 180.0,  90.0 ]
            for i in range(len(Transx)):
                Trans = [ Transx[i], Transy[i], BoxH*0.5-3.0*absLsz ]
                Rot = [ Rotx[i], Roty[i], Rotz[i] ]
                TransPres = AddAxisLabel(LabelTmp, Trans, Rot)
                
    # for coordinate labels
    Z = Bottom - absLsz*3.0
    #latitude grid
    if len(lats) > 0:
        for lat in lats:
            PlaneTmp = AddLatPlane(lat, bounds, ratios, basis, AxisColor=AxisColor, AxisWidth=AxisWidth)
            RenameSource("Plane"+str(lat),PlaneTmp)
            if abs(LabelSize) > 0.0 :
                (label,transa,transb) = AddLatLabel(lat, absLsz, bounds, ratios, basis, AxisColor)
                RenameSource("Label"+str(lat),label)
                Show(transa)
                Show(transb)
        #Latitude label
        if LabelSize > 0.0 :
            LabelTmp = a3DText(Text='latitude')
            RenameSource("LatLabel",LabelTmp)
            midLat = 0.5*(Far-Near)
            Trans = [Right, midLat-2.5*absLsz, Z]
            Rot = [ 90.0, 90.0, 0.0 ]
            TransLat = AddAxisLabel(LabelTmp, Trans, Rot)
            
            Trans = [Left, midLat+2.5*absLsz, Z]
            Rot = [ 90.0,-90.0, 0.0 ]
            TransLat = AddAxisLabel(LabelTmp, Trans, Rot)

    #longitude grid
    if len(lons) > 0 :
        for lon in lons:
            PlaneTmp = AddLonPlane(lon, bounds, ratios, basis, AxisColor=AxisColor, AxisWidth=AxisWidth)
            RenameSource("Plane"+str(lon)+"E",PlaneTmp)
            if abs(LabelSize) > 0.0 :
                (label,transa,transb) = AddLonLabel(lon, absLsz, bounds, ratios, basis, AxisColor)
                Show(transa)
                Show(transb)

        #Longitude label
        if LabelSize > 0.0 :
            LabelTmp = a3DText(Text='longitude')
            RenameSource("LonLabel",LabelTmp)
            Trans = [0.5*(Left+Right)-3.0*absLsz, Near, Z]
            Rot = [ 90.0,   0.0, 0.0 ]
            TransLon = AddAxisLabel(LabelTmp, Trans, Rot)
            
            Trans = [0.5*(Left+Right)+3.0*absLsz, Far, Z]
            Rot = [ 90.0, 180.0, 0.0 ]
            TransLon = AddAxisLabel(LabelTmp, Trans, Rot)






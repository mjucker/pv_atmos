import unittest
import random
import math
import numpy as np

#adjust the path to pv-atmos here:
pvAtmosPath='../'
try: from atmos_basic import *
except: execfile(pvAtmosPath + 'atmos_basic.py')
try: from atmos_grids import *
except: execfile(pvAtmosPath + 'atmos_grids.py')


class TestSequenceFunctions(unittest.TestCase):
    def setUp(self):
        self.p = random.uniform(1e-5,1e3) # choose pressure surface randomly
        self.rat = [random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)] # choose box aspect ration randomly
        self.pCoord = round(random.uniform(0,1)) # choose whether or not data uses pressure coordinates
        self.np = round(random.uniform(1,5)) # number of pressure levels for grid
        self.nx = round(random.uniform(1,5)) # number of lon gridlines for grid
        self.ny = round(random.uniform(1,5)) # number of lat gridlines for grid
        mx1 = random.uniform(-90,90)
        mx2 = random.uniform(-90,90)
        self.mx = [min(mx1,mx2),max(mx1,mx2)] # range of latitude for grid
        my1 = random.uniform(0,360)
        my2 = random.uniform(0,360)
        self.my = [min(my1,my2),max(my1,my2)]  # range of longitude for grid
        self.top = random.uniform(1e-5,self.p) # top of the box
        self.rad = random.uniform(1e-5,10) # radius of spherical geometry

    def test_basic(self):
        # load data
        fileName = pvAtmosPath+'examples/uv_daily.nc'
        (output_nc,CorrZ,Coor,AspRat) = loadData(fileName, ['pfull','lat','lon'], self.pCoord, self.rat, math.ceil(self.p))
        self.assertEqual(output_nc.OutputType,'Unstructured')
        (W,norm,clipS,clipN)=CartWind2Atmos(src=AspRat,ratios=self.rat)
        self.assertEqual(clipS.ClipType.Origin[1],-80.0*self.rat[1])
        self.assertEqual(clipN.ClipType.Origin[1], 80.0*self.rat[1])

    def test_allGrids(self):
        AddGrid(press=np.arange(self.p,self.top,-self.p/self.np), lats=np.arange(self.mx[0],self.mx[1],(self.mx[1]-self.mx[0])/self.nx), lons=np.arange(self.my[0],self.my[1],(self.my[1]-self.my[0])/self.ny), bounds=[self.my[0],self.my[1],self.mx[0],self.mx[1],self.p,self.top], ratios=self.rat, basis=math.ceil(self.p), AxisColor=[0,0,0], AxisWidth=2.0,LabelSize=5.0)
        print 'Check pipeline for grid'

    def test_shells(self):
        # load data
        fileName = pvAtmosPath+'examples/uv_daily.nc'
        (output_nc,CorrZ,Coor,AspRat) = loadData(fileName, ['pfull','lat','lon'], 1, self.rat, math.ceil(self.p))
        AtmosShells(radius=self.rad,ratios=self.rat,basis=math.ceil(self.p),src=AspRat,shellValues=np.arange(self.p,self.top,-self.p/self.np),labels=self.pCoord,labelPosition=[self.my[0],self.mx[0]],waterMark='WaterMark',markPosition=[self.my[1],self.mx[1]])
        print 'Check pipeline for spherical shells'


if __name__ == '__main__':
    unittest.main()


suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)

runner = unittest.TextTestRunner()
runner.run(suite)

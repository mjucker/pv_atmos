 #!/usr/bin/python
# Filename: tools.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al
#
# This file provides helper functions that can be useful as pre-viz-processing of files and data

############################################################################################

## create a generic netCDF file that can be read with pv_atmos
#
def WriteGenericNc(x, y, z, t, data, varName, outFileName='ncGeneric.nc'):
    """Write netCDF file that is compatible with pv_atmos.
        
        Takes one variable, and stores it inside a generic netCDF file.
        If the file already exists, adds the variable to the file.
        
        INPUTS:
        x,y,z,t     - spatial and time coordinates. If not all used, set to []
        data        - data to be written. Should be one variable.
        varName     - name of the variable data in the new netCDF file
        outFileName - name of (new) netCDF file. If exists already, variable is added.
        NOTE: if outFileName already exists, the dimensions in the file have to
        correspond to the dimensions given as inputs.
        """
    import netCDF4 as nc
    import os.path as op
    
    dims = ()
    if op.isfile(outFileName):
        mode = 'r+'
    else:
        mode = 'w'
    outFile = nc.Dataset(outFileName, mode, format='NETCDF3_64BIT')
    if len(x) > 0:
        xD = outFile.createDimension('x',len(x))
        xV = outFile.createVariable('x','f4', ('x',))
        xV[:] = x
        xV.setncattr('long_name','X axis')
        xV.setncattr('cartesian_axis','X')
        dims = ('x',) + dims
    if len(y) > 0:
        yD = outFile.createDimension('y',len(y))
        yV = outFile.createVariable('y','f4', ('y',))
        yV[:] = y
        yV.setncattr('long_name','Y axis')
        yV.setncattr('cartesian_axis','Y')
        dims = ('y',) + dims
    if len(z) > 0:
        zD = outFile.createDimension('z',len(z))
        zV = outFile.createVariable('z','f4', ('z',))
        zV[:] = z
        zV.setncattr('long_name','Z axis')
        zV.setncattr('cartesian_axis','Z')
        dims = ('z',) + dims
    if len(t) > 0:
        tD = outFile.createDimension('time',len(t))
        tV = outFile.createVariable('time','f4', ('time',))
        tV[:] = t
        tV.setncattr('long_name','time')
        tV.setncattr('cartesian_axis','T')
        tV.setncattr('units','time since 0000-00-00 00:00:00')
        dims = ('time',) + dims
    vOut = outFile.createVariable(varName,'f4',dims)
    vOut[:] = data
    outFile.close()
    print 'Done, created file'+outFileName



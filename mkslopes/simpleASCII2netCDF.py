#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np
import sys

try:
    sa_file = sys.argv[1]
except IndexError:
    sys.exit(f"File and variable names needed: '{sys.argv[0]} FILE VARIABLE'")

print(sa_file)

try:
    vname = sys.argv[2]
except IndexError:
    sys.exit(f"File and variable names needed: '{sys.argv[0]} FILE VARIABLE'")

# Load ParFlow simple ASCII data
with open(sa_file, 'r') as f:
    nx, ny, nz = map(int, f.readline().split())
var = np.loadtxt(sa_file, skiprows=1)
myvar = var.reshape((nx, ny, nz), order='C')

# Write to netCDF
nc_file = f"{vname}.nc"
with Dataset(nc_file, 'w', format='NETCDF4') as nc:
    nc.createDimension('lon', nx)
    nc.createDimension('lat', ny)
    nc.createDimension('time', nz)
    nc.description = 'Converted from ParFlow simple ASCII to netCDF'
    var = nc.createVariable(f'{vname}', 'f4', ('time', 'lat', 'lon'))
    var[:] = myvar

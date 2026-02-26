#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np

nc = Dataset('HSURFBurnedAndMod.nc', 'r')
var = nc.variables['HSURFBurnedAndMod'][:]
nc.close()
var_fortran = np.array(var, order='F')
with open('HSURFBurnedAndMod.sa', 'w') as f:
     f.write(f"{var.shape[1]} {var.shape[1]} {var.shape[0]}\n1\n")
     for val in var_fortran.flatten():
         f.write(f"{val}\n")

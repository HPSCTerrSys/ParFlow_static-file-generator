#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np
import sys

try:
    name = sys.argv[1]
except IndexError:
    sys.exit(f"Argument needed: '{sys.argv[0]} myvar' with myvar a variable in myvar.nc.")

nc = Dataset(f"{name}.nc", 'r')
var = nc.variables[f"{name}"][:]
nc.close()
var_fortran = np.array(var, order='F')
with open(f"{name}.sa", 'w') as f:
     f.write(f"{var.shape[1]} {var.shape[1]} {var.shape[0]}\n1\n")
     for val in var_fortran.flatten():
         f.write(f"{val}\n")

#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np
import sys

try:
    fname = sys.argv[1]
except IndexError:
    sys.exit(f"File and variable names needed: '{sys.argv[0]} FILE VARIABLE'")

try:
    vname = sys.argv[2]
except IndexError:
    sys.exit(f"File and variable names needed: '{sys.argv[0]} FILE VARIABLE'")

nc = Dataset(f"{fname}.nc", 'r')
var = nc.variables[f"{vname}"][:]
nc.close()
var_fortran = np.array(var, order='F')
with open(f"{fname}.sa", 'w') as f:
     f.write(f"{var.shape[1]} {var.shape[1]} {var.shape[0]}\n1\n")
     for val in var_fortran.flatten():
         f.write(f"{val}\n")

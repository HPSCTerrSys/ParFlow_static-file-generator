#!/usr/bin/env python3
import netCDF4 as nc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import sloth.IO 

LLSMFileName = '../mklandmask/EUR-11_TSMP_FZJ-IBG3_444x432_LAND-LAKE-SEA-MASK.nc'
LLSMVarName  = 'LLSM'

with nc.Dataset(LLSMFileName, "r") as nc_file:
    LLSM = nc_file.variables[LLSMVarName][0,...]

ny, nx = LLSM.shape
nz = 15
dy = dx = 12500.0
dz = 2.0

# write PFBMask:
# 0 where LLSM is land or leak
# 1 where LLSM is ocean
PFBMask = np.where((LLSM<1), 0, 1)
print(f'PFBMask.shape: {PFBMask.shape}')

plt.imshow(PFBMask, origin='lower', interpolation='none')
plt.colorbar()
plt.savefig('./PFBMask.pdf')

## Broadcast PFBMask from (y,x) to (z,y,x)
#PFBMask = np.broadcast_to(PFBMask, (nz, ny, nx))
#print(f'PFBMask.shape: {PFBMask.shape}')

# save as .pfb
PFBMaskFileName = './PfbMask4SolidFile.pfb'
#sloth.IO.create_pfb(PFBMaskFileName, PFBMask[...], delta=(dz,dy,dx))
sloth.IO.create_pfb(PFBMaskFileName, PFBMask[np.newaxis,...], delta=(dz,dy,dx))

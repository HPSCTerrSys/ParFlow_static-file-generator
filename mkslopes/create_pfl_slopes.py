#!/usr/bin/env python3
''' 
author: Marco VAN HULTEN
e-mail: Marco.van.Hulten@uni-bonn.de
version: 2026-05-29

Description:
create_pfl_slopes.py creates slopes from a Digital Elevation Model (DEM).

Usage:
    ./create_pfl_slopes.py

INPUT FILES:
- DEM (topography, HSURF)
- Shapefiles holding the correct river-network (multiple files posible)
- land-lake-see mask
- target Grid (lats, lons)

OUTPUT FILES:
- x and y slopes
'''
import sys
from pathlib import Path
import numpy as np
import netCDF4 as nc
import priority_flow as prf
import matplotlib.pyplot as plt

# Read into array
ncin = nc.Dataset("HSURF.nc", "r")
data = ncin.variables['HSURF'][:]
np.shape(data)

# Write to a npy data file
np.save('data/DEM_tsmp.npy', np.array(data))
ncin.close()

# Load TSMP data
exec( open("tsmp_data_loader.py").read() )
DEM_tsmp = load_dem()
DEM=DEM_tsmp[::-1]  # so DEM is our direction-corrected DEM!
np.shape(DEM)
exec( open("plotting.py").read() )

nc_watershed_mask = nc.Dataset("../mklandmask/EUR-11_TSMP_FZJ-IBG3_444x432_LAND-LAKE-SEA-MASK.nc")
watershed_mask = nc_watershed_mask.variables['LLSM'][0,:,:]
watershed_mask[watershed_mask==2] = 1.
watershed_mask2 = watershed_mask[::-1]
_plot_inputs(watershed_mask2)

# Flow direction
init = prf.init_queue(DEM)
trav_hs = prf.d4_traverse_b( DEM, init["queue"].copy(),
                             init["marked"].copy(),
                             basins=init["basins"].copy(),
                             epsilon=0,
                             n_chunk=10,
)
dem_diff = trav_hs["dem"] - DEM
dem_diff[dem_diff==0] - np.nan
targets = init["marked"].copy()
targets[targets==0] = np.nan
_plot_step1(trav_hs, dem_diff, targets, watershed_mask2)

np.save("data/flow_direction.npy", trav_hs["direction"])

# Create the slopes
# TODO: use prf.slope_calc_upwind() for downwinding (needed in TSMP2)
slopes_std = prf.slope_calc_standard(
    dem=DEM_tsmp.copy()[::-1],
    direction=trav_hs["direction"],
    mask=watershed_mask2.copy(),
    minslope=1e-5,
    maxslope=1,
    dx=12000, dy=12000,
    secondary_th=-1,
)
slopex = slopes_std["slopex"]
slopey = slopes_std["slopey"]
_plot_slopes(slopex, slopey)

# Write the slope files to ParFlow pfb format
from parflow.tools.io import read_pfb, write_pfb
write_pfb("slopex.pfb", slopex)
write_pfb("slopey.pfb", slopey)
write_pfb("flow_direction.pfb", trav_hs["direction"])

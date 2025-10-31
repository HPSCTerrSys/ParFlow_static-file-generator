#!/usr/bin/env python3
import numpy as np
import csv
import os
import netCDF4 as nc
from datetime import datetime
import sloth.mapper

###############################################################################
# THE USER HAS TO CHANGE BELOW ################################################
###############################################################################
EXTPRA_IN  = '${RAWDATA_DIR}/domain202210280000.nc'
EXTPRA_OUT = './EUR-11_TSMP_FZJ-IBG3_464x452_EXTPAR.nc'
# pixel to adjust, as e.g. Bosporus break through
adjPxIN    = f'./Breakthrough-Bosporus.csv'

###############################################################################
# THE USER HAS TO CHANGE ABOVE ################################################
###############################################################################

# remove output file if exist, to avoid interference
try:
    os.remove(EXTPRA_OUT)
except OSError:
    pass
os.system(f'cp -v {EXTPRA_IN} {EXTPRA_OUT}')

# pixel to adjust
# Read in vertices of regions to adjust
with open(adjPxIN, "r") as f:
    ToAdjust = {}
    reader = csv.reader(f, delimiter=",")
    next(reader, None)  # skip the headers
    for line in reader:
        regionID = line[0]
        # check that regionID is a key in ToAdjust dict and if not create
        if regionID not in ToAdjust.keys():
            ToAdjust[regionID] = {'id': [], 'lon':[], 'lat':[], 
                                  'lonIdx':[], 'latIdx':[]}
        ToAdjust[regionID]['id'].append(line[1])
        ToAdjust[regionID]['lon'].append(float(line[2]))
        ToAdjust[regionID]['lat'].append(float(line[3]))

with nc.Dataset(EXTPRA_OUT, 'a') as nc_file:
    lonIN   = nc_file.variables['lon'][...]
    latIN   = nc_file.variables['lat'][...]
    FR_LAND = nc_file.variables['FR_LAND'][...]

    # Initialize mapper from SLOTH
    Mapper = sloth.mapper.mapper(SimLons=lonIN, SimLats=latIN)
    for regionID in ToAdjust.keys():
        Mapper.ObsLons = np.array(ToAdjust[regionID]['lon'])
        Mapper.ObsLats = np.array(ToAdjust[regionID]['lat'])
        Mapper.ObsIDs  = np.array(ToAdjust[regionID]['id'])
        # Map points to sim coords / get array index of points in sim coords
        Mapper.MapRaw()
        print(f'Found for regionID {regionID}:')
        print(f'Mapped x index: {Mapper.MapXIdx_raw}')
        print(f'Mapped y index: {Mapper.MapYIdx_raw}')
        ToAdjust[regionID]['lonIdx'] = Mapper.MapXIdx_raw.tolist()
        ToAdjust[regionID]['latIdx'] = Mapper.MapYIdx_raw.tolist()
    
    # Adjust pixel
    for regionID in ToAdjust:
        print(f'handling regionID {regionID}')
        CoordIdx = list(zip(ToAdjust[regionID]['latIdx'], ToAdjust[regionID]['lonIdx']))
        LATref, LONref = CoordIdx[0]
        for Plat, Plon in CoordIdx[1:]:
            print(f'handling: {Plat} | {Plon}')
            # set FR_LAND to 0 (no land - only water)
            FR_LAND[Plat, Plon] = 0

    # overwrite FR_LAND in output file
    nc_file["FR_LAND"][...] = FR_LAND[...]


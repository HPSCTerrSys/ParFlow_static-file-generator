#!/usr/bin/env python3
''' 
author: Niklas WAGNER
e-mail: n.wagner@fz-juelich.de
version: 2021-03-19

It could be needed, that some pixel got 'adjustet' to correctly place river
streams. The thing is, that the GRASS algo in later sva* script is trying
(or forcing) to route water out of the domain. If now rivers do flow e.g.
to the blacksea, and the blacksea is (for coarse resolved topos) not
connected to the mediteraniensea, than GRASS will route water through another
way e.g. through a big river in eastern europe. Saying the flow dir of this river
is getting reverted in order to route water out of the blacksea (even ocean
pixel are masked within ParFlow). To solve this one can 'break through' the
land pixel between black and mediteraniensea, that GRAS is routing water this
way. One can think of other needed 'adjustments' which gets defined below.
'''
import numpy as np
import netCDF4 as nc
import sys
import os
import csv
import fiona
from shapely.geometry import shape
import sloth.mapper

###############################################################################
#### START OF: stuff you need to adjust!
###############################################################################
# Topo to 'adjust'
with nc.Dataset(f'./HSURFBurned.nc', 'r') as nc_file:
    HSURFin = nc_file.variables['HSURFBurned'][...]
    lon2D   = nc_file.variables['lon'][...]
    lat2D   = nc_file.variables['lat'][...]
# file to save adjusted topo to
nc_file      = nc.Dataset(f'./HSURFBurnedAndMod.nc', 'w')
latDim       = nc_file.createDimension('lat',HSURFin.shape[0])
lonDim       = nc_file.createDimension('lon',HSURFin.shape[1])
out_HSURF    = nc_file.createVariable(f'HSURFBurnedAndMod','f8',('lat','lon'),
                                           zlib=True)
# pixel to adjust
# Read in vertices of regions to adjust
with open("./BurnRiversAndCanyons_Vertices_EUR11.csv", "r") as f:
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

# Initialize mapper from SLOTH
#Mapper = sloth.mapper.mapper(SimLons=lon2D, SimLats=lat2D,
#                            ObsLons=pointLons, ObsLats=pointLats,
#                            ObsIDs=pointIDs)
Mapper = sloth.mapper.mapper(SimLons=lon2D, SimLats=lat2D)
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
# simple pass one list per 'region' to adjust with a leading ref-pixel.
# lets say you need to break through land pixel between blacksea and 
# mediteraniensea. 'breaking' through is done by setting topo height of the 
# pixels to break through to a ref hight - the sea hight.
#ToAdjust = [
#        [(138,339), (137,339), (136,339)],
#        [(246,221),(246,220),(247,220),(247,219),(248,219),(248,218),(249,218),(249,217),(250,217),(250,216)]
#        ]
###############################################################################
#### END OF: stuff you need to adjust!
###############################################################################

tmp_out_HSURF = HSURFin.copy()
for regionID in ToAdjust:
    print(f'handling regionID {regionID}')
    CoordIdx = list(zip(ToAdjust[regionID]['latIdx'], ToAdjust[regionID]['lonIdx']))
    LATref, LONref = CoordIdx[0]
    for Plat, Plon in CoordIdx[1:]:
        print(f'handling: {Plat} | {Plon}')
        tmp_out_HSURF[Plat, Plon] = tmp_out_HSURF[LATref,LONref]

out_HSURF[...] = tmp_out_HSURF[...]
nc_file.close()


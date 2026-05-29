#!/usr/bin/env python3
''' 
author: Niklas WAGNER
e-mail: n.wagner@fz-juelich.de
version: 2021-03-18

Description:
burnShape2Topo.py is aimed to burn the 'correct' river-corridors on top of a
given topography.
This is needed because coarse topography resolutions poorly reflect the correct 
river streams and is done in two steps:
  1) mapping the ‘correct’ river-corridors provided by official river-network 
     datasets (as e.g. HydroSHEDS) from shape-files onto the topography raster.  
  2) lowering the topography for pixels where a river stream is detected as 
     well as surrounding pixel, to force the model to drain towards the correct 
     river-stream.

Usage:
INPUT FILES:
- HSURF (height of surface = topo) file to burn river on
- Shapefiles holding the correct river-network (multiple files posible)
- land-lake-see mask
- target Grid (lats, lons)

OUTPUT FILES:
- HSURF used to burn rivers on
- HSURF with burned rivers
'''
import numpy as np
import netCDF4 as nc
import sys
import os
import fiona
import datetime

from shapely.geometry import shape
import sloth.IO

###############################################################################
#### START OF: stuff you need to adjust!
###############################################################################
# Target grid
# - this file should contian the correct lat/lon values of the target grid
# - def the file as well the name of lat and lon variables inside
rawdata_dir = os.environ['RAWDATA_DIR']
with nc.Dataset(f'{rawdata_dir}/EUR-11_TSMP_FZJ-IBG3_CLMPFLDomain_444x432_gridfile.nc','r') as nc_trgGridFile:
    lons = nc_trgGridFile.variables['lon'][...]
    lats = nc_trgGridFile.variables['lat'][...]
# - `griddesFileName` is the CDO griddes files from the .nc file above
griddesFileName = f'{rawdata_dir}/EUR-11_TSMP_FZJ-IBG3_CLMPFLDomain_444x432_griddes.txt'

# River network shape file to burn to DEM
shapefiles = [
      f'HydroRIVERS_v10_eu_shp',
      f'HydroRIVERS_v10_af_shp',
      ]

# Origin topo
# - Def how many pixels to cut from the boundary (NBCCUT)
# - - This is e.g. needed, if you take the topo or lon/lat values from COSMO, but
# - this files should contain the original topot you want to burn rivers in
# - define also the name of the topo var
# - if this file holdes you lan/lot values also, simple define the file twice
NBCCUT = 10
nc_HSURF = nc.Dataset(f'../mklandmask/EUR-11_TSMP_FZJ-IBG3_464x452_EXTPAR.nc')
HSURF    = nc_HSURF.variables['HSURF'][NBCCUT:-NBCCUT,NBCCUT:-NBCCUT]
print(f'HSURF.shape: {HSURF.shape}')
# Lans-Lake-Sea mask (LLSM)
with nc.Dataset(f'../mklandmask/EUR-11_TSMP_FZJ-IBG3_444x432_LAND-LAKE-SEA-MASK.nc', 'r') as nc_indic:
    LLSM = nc_indic.variables['LLSM'][0,...]
    print(f'LLSM.shape: {LLSM.shape}')
mask_lsm = np.where(LLSM<2, True, False)
###############################################################################
#### END OF: stuff you need to adjust!
###############################################################################

def get_majorRivers(rec):
    ''' filter records of shp for major rivers only.

    ORD_FLOW is a indicator within HydroSHEDS data, to distinguish rivers
    into logarithmic size classes based on averaged flowrates. 
    This may varrii for different datasets, but should be easily adjustable.

    Important is, that this function returns the needed records of the shp
    of interest only.

    INPUT:  records of a shp opend with the py-lib fiona
    RETURN: filtered records.
    '''
    propertie = 'ORD_FLOW'
    validValues = [1,2,3,4]
    return rec['properties'].get(f'{propertie}') in validValues 

def spher_dist_v1(lon1, lat1, lon2, lat2, Rearth=6373):
    """ calculate the spherical / haversine distance

    Source: https://www.kompf.de/gps/distcalc.html
    This function is supposed to proper handle different shaped coords
    latX and lonX is supposed to be passed in rad

    return 2D ndarray
    """
    term1 = np.sin(lat1) * np.sin(lat2)
    term2 = np.cos(lat1) * np.cos(lat2)
    term3 = np.cos(lon2 - lon1)
    return Rearth * np.arccos(term1+term2*term3)

def smothFct_v1(x, strength):
    ''' calculating the smooth coefficient to burn rivers.

    This function is a simple helper function for BurnRiversSmoth().
    The idea of BurnRiversSmoth() is, to push down the topography fo 
    rivers location and sorounding pixels. For the actual river-corridor 
    this effect of pushing down should be strongest and decrease with each
    pixel further away from river-corridor, to 'fade the effect out'.

    INPUT
      x:        count of pixels away from actuall river-coridor
      strength: a factor controling how strong to push topo down
    RETURN
      factor of reducing / pushing down topo in % (0-1) of total hight  
    '''
    return (1-(strength/x))

def smothFct_v2(denominator, nominator=2.):
    ''' calculating the smooth coefficient to burn rivers.

    This function is a simple helper function for BurnRiversSmoth().
    The idea of BurnRiversSmoth() is, to push down the topography fo 
    rivers location and sorounding pixels. For the actual river-corridor 
    this effect of pushing down should be strongest and decrease with each
    pixel further away from river-corridor, to 'fade the effect out'.

    INPUT
      x:         count of pixels away from actuall river-coridor
      nominator: a factor controling how strong this fades out with x
    RETURN
      factor of reducing / pushing down topo in % (0-1) of fix burn-depth  
    '''
    return (float(nominator)/denominator)

def extendNPix(data, N=1):
    ''' simple function to increase river width

    By shiftig a binary array (a mask indicating river or not) 1-pixel north, 
    south, west, and east and sum up all five arrays, I do recieve an binary 
    array where each line (river) is 2 pixel bigger.

    INPUT
      data: the binary array (rivers)
      N=1 : amount of shifting the array - should be one!
    RETURN
      binary array (river mask) with extendes lines (rivers) 
    ''' 
    noth  = np.roll(data, N, axis=0)
    east  = np.roll(data, N, axis=1)
    south = np.roll(data, -N, axis=0)
    west  = np.roll(data, -N, axis=1)
    data  += noth + east + south + west
    data[data != 0] = 1
    return data

def BurnRiversSmoth(topo, rivers, strength=0.2, radius=10):
    ''' burn rivers based on pixels actual topo height

    INPUT
      topo:     ndarray in target raster with topo
      rivers:   ndarray with binary river mask in same shape as topo (1=river; 0=no river)
      strength: initial strength of burning rivers (0.1= push down river-coridors by 10% initially)
      radius:   number of pixel around rivers to also push down (decreased intensity)
    RETURN
      topo with 'burned' rivers in same shape as topo input
    '''
    # copy input array to not change the original
    rivers_smoth = rivers.copy()
    out_topo = topo.copy()
    
    # actuall burning
    # repeat the actuall 'burning' (reduzing of topo) for 'radius'
    # with each itteration (1 pixel further away from river-corridor)
    # the intensity of 'burning' is reduced by smothFct().
    for i in range(1,radius+1):
        out_topo[rivers_smoth==1] *= smothFct_v1(x=i, strength=strength)
        print(f'smothFct_v1: {smothFct_v1(x=i, strength=strength)}')
        rivers_smoth = extendNPix(rivers_smoth, N=1)
    return out_topo

def BurnRiversSmoth_v2(topo, rivers, baseBurnDepth=10, radius=10):
    ''' burn rivers by fix 'brun-depth'

    INPUT
      topo:           ndarray in target raster with topo
      rivers:         ndarray with binary river mask in same shape as topo (1=river; 0=no river)
      baseBurnDepth:  initial depth of burning rivers (10= push down river-coridors by 10[L] (L depends on topo units))
      radius:         number of pixel around rivers to also push down (decreased intensity)
    RETURN
      topo with 'burned' rivers in same shape as topo input
    '''
    # copy input array to not change the original
    rivers_smoth = rivers.copy()
    out_topo = topo.copy()
    
    # actuall burning
    # repeat the actuall 'burning' (reduzing of topo) for 'radius'
    # with each itteration (1 pixel further away from river-corridor)
    # the intensity of 'burning' is reduced by smothFct().
    for i in range(1,radius+1):
        out_topo[rivers_smoth==1] -= baseBurnDepth * smothFct_v2(denominator=i, nominator=2)
        print(f'smothFct_v2: {smothFct_v2(denominator=i, nominator=2)}')
        rivers_smoth = extendNPix(rivers_smoth, N=1)
    return out_topo

def getShapefileOnGrid(shp_file, RiverMaskIn, target_lats, target_lons):
    ''' function to get shapfile entries on target grid / raster

    This is a brutforce function to bring river location from shp-file
    on target raster. Therefore each entry / record of the shapefile is 
    itterated and marked as 'river' on a binary array (0,1) of the same shape
    as to topo. The related pixel index is found via minimum spheric distance 
    between each record and the target lan/lot grids.
    'Riverpixel' may be marked multiple times and got overwritten if.
    Depending on the nuber of rivers this may take a while!

    INPUT
      shp_file:    full path to shp-file (or dir/zip holding the shp-file with same name)
      RiverMask:   the array to be filled (have to have same shape as target_lats and target_lons)
      target_lats: array of same shape as topo holding lat values for each pixel
      target_lons: array of same shape as topo holding lon values for each pixel
    RETURN
      filled RiverMask
    '''
    # init RiverMask in shape of target grid    
    RiverMask = RiverMaskIn
    # need min/max lon/lat values to exclude shp-records outside the domain
    minLon = np.nanmin(target_lons)
    maxLon = np.nanmax(target_lons)
    minLat = np.nanmin(target_lats)
    maxLat = np.nanmax(target_lats)

    # open shp-file
    with fiona.open(f'{shp_file}') as src:
        # keep a filtered set of rivers only (see get_majorRivers() for details)
        hits = filter(get_majorRivers, src)

        counter = 0
        # itterate ofer all records (stream segments in our case)
        for ssmt in hits:
            # coordinates are stored in (lon, lat)
            for Plon, Plat in ssmt['geometry']['coordinates']:
                # check if ssmt is in domain
                # if all or msot river are inside the domain, this slows stuff
                # done, but keep things faster is many rivers are outside the domain.
                if not ( (minLon <= Plon <= maxLon) and (minLat <= Plat <= maxLat) ):
                    continue
                # calculate spheric distance to all points in target grid
                dist = spher_dist_v1(np.deg2rad(lons),
                                     np.deg2rad(lats),
                                     np.deg2rad(Plon),
                                     np.deg2rad(Plat))
                # get index with min distance
                mapped_idx = np.unravel_index(np.argmin(dist, axis=None),dist.shape)
                # mark found point as river
                RiverMask[mapped_idx[0], mapped_idx[1]] = 1
                # some print-steatments to let the user know the program is running...
                print(f'handled point {counter:05d}: ({Plat}, {Plon})')
                print(f'mapped_idx: ({mapped_idx[0]},{mapped_idx[1]})')
                counter += 1
    return RiverMask




# as bringing shp file to target raster may take some time, the river-mask
# is dumped to 'RiverMask.npy' after calculation. So if you run the script
# again, it is frist tryed to open this dumped file to save some calculatio
# time. If this file is not found -> recalculate
RiverMaskFileName = './RiverMask.nc'
try:
    with nc.Dataset(RiverMaskFileName, 'r') as nc_file:
        RiverMask = nc_file.variables['RiverMask'][...]
    print(f'RiverMask.shape: {RiverMask.shape}')
except FileNotFoundError:
    RiverMask = np.zeros_like(lons.filled())
    #shapefiles are defined at the top
    for shapefile in shapefiles:
        RiverMask = getShapefileOnGrid(shp_file=shapefile, RiverMaskIn=RiverMask, target_lats=lats, target_lons=lons)

    descriptor = [
            f'RiverMask -- add some extra info here!',
            ]
    netCDFFileName = sloth.IO.createNetCDF(RiverMaskFileName, domain=griddesFileName,
            author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
            institution='FZJ/IBG-3', history=f'Created: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}',
            description=' '.join(descriptor),
            source='')

    with nc.Dataset(netCDFFileName, 'a') as nc_file:
        nc_RiverMakse = nc_file.createVariable('RiverMask', 'i4', ('rlat', 'rlon'),
                zlib=True)
        nc_RiverMakse.standard_name = "RiverMask"
        nc_RiverMakse.long_name = "RiverMask"
        nc_RiverMakse.units = "--"
        nc_RiverMakse.coordinates = "lon lat"
        nc_RiverMakse.grid_mapping = "rotated_pole"
        print(f'nc_RiverMakse.shape: {nc_RiverMakse.shape}')
        print(f'RiverMask.shape: {RiverMask.shape}')
        nc_RiverMakse[...] = RiverMask[...]

# open HSURFBurned.nc to write new topo with burned river into
descriptor = [
        f'HSURFBurned -- original HSURF with burned major rivers',
        f'Add more description here',
        ]
netCDFFileName_HSURFBurned = sloth.IO.createNetCDF('./HSURFBurned.nc', domain=griddesFileName,
        calcLatLon=True,
        author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
        institution='FZJ/IBG-3', history=f'Created: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}',
        description=' '.join(descriptor),
        source='')
netCDFFileName_HSURF = sloth.IO.createNetCDF('./HSURF.nc', domain=griddesFileName,
        author='Niklas WAGNER', contact='n.wagner@fz-juelich.de',
        institution='FZJ/IBG-3', history=f'Created: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}',
        description='The origin HSURF used to burn in rivers',
        source='')

# working copys
tmp_HSURF = HSURF.copy()
tmp_HSURFBurned = HSURF.copy()
# Ofset-adjustment to avoid negativ height (lowering the topo by % is hard for negativ values)
tmp_offset = 999
tmp_HSURFBurned[mask_lsm==False] += tmp_offset

# Simplest way to burn rivers: set rivers to fix height.
# However this does not work well, and I guess this is because of 
# how the GRAS algo in sva*-script works.
#burnHeight = 0.0
#tmp_HSURFBurned[Rivers==1] = burnHeight
# the bit mor complex way to 'burn' rivers tot topo. For details see BurnRiversSmoth() 
# and BurnRiversSmoth_v2() above
tmp_HSURFBurned = BurnRiversSmoth(topo=tmp_HSURFBurned, rivers=RiverMask, strength=0.3, radius=3)
#tmp_HSURFBurned = BurnRiversSmoth_v2(topo=tmp_HSURFBurned, rivers=RiverMask, baseBurnDepth=10, radius=10)

# redo offset-adjustment
tmp_HSURFBurned[mask_lsm==False] -= tmp_offset
# make sure ocean is lowest point!
tmp_HSURFBurned[mask_lsm==True] = -tmp_offset - 50

# write results to netCDF object
with nc.Dataset(netCDFFileName_HSURFBurned, 'a') as nc_file:
    nc_HSURFBurned = nc_file.createVariable('HSURFBurned', 'f8', ('rlat', 'rlon'),
            zlib=True)
    nc_HSURFBurned.standard_name = "HSURFBurned"
    nc_HSURFBurned.long_name = "HSURFBurned"
    nc_HSURFBurned.units = "m"
    nc_HSURFBurned.coordinates = "lon lat"
    nc_HSURFBurned.grid_mapping = "rotated_pole"
    nc_HSURFBurned[...] = tmp_HSURFBurned
with nc.Dataset(netCDFFileName_HSURF, 'a') as nc_file:
    nc_HSURF = nc_file.createVariable('HSURF', 'f8', ('rlat', 'rlon'),
            zlib=True)
    nc_HSURF.standard_name = "HSURF"
    nc_HSURF.long_name = "HSURF"
    nc_HSURF.units = "m"
    nc_HSURF.coordinates = "lon lat"
    nc_HSURF.grid_mapping = "rotated_pole"
    nc_HSURF[...] = HSURF[...]

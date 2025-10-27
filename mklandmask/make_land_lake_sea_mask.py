#!/usr/bin/env python3
#
# This script builds a land-lake-sea mask (LLSM) on the target grid, with 
# land=2, lake=1, sea=0. 
# The calculation is based on INT2LM / COSMO static file, in particular
# the variables `FR_LAND` and `HSURF`, to ensure 100% matching of the LLSM
# between all component models.
#
# EXECUTION
# - Update below marked variables according the used domain
# - Execution : python make_land_lake_sea_mask.py
#
# INPUT FILES
# - Target reference grid (output), netDCF
# - EXTPAR / INT2LM static file (input), netCDF
#
# OUTPUT FILE
# - LLSM on target grid, netCDF
#
# STEPS
# 1 - Read input data
# 3 - Make land-lake-sea mask
# 4 - Optimize land-lake-sea mask, i.e. change 'lake' values to 'sea' values 
#     along the coasts and in the estuaries
# 5 - Save output
#
# AUTHOR: Alexandre BELLEFLAMME, Niklas WAGNER
# E-MAIL: a.belleflamme@fz-juelich.de, n.wagner@fz-juelich.de
# INSTITUTION : Forschungszentrum Juelich GmbH - IBG-3 Agrosphere 
#               (Integrated Modelling Group)
# VERSION : 2022-10-28
###############################################################################

import numpy as np
import csv
from netCDF4 import Dataset
from datetime import datetime
import sloth.mapper
import matplotlib as mpl
import matplotlib.pyplot as plt

###############################################################################
# THE USER HAS TO CHANGE BELOW ################################################
###############################################################################
ncGridIN  = f'../01_Grids/EUR-11_TSMP_FZJ-IBG3_CLMPFLDomain_444x432_gridfile.nc'  # grid file
ncLsmIN   = f'./EUR-11_TSMP_FZJ-IBG3_464x452_EXTPAR.nc'              # land sea mask file
ncOUT     = f'./EUR-11_TSMP_FZJ-IBG3_444x432_LAND-LAKE-SEA-MASK.nc'  # output file
nboundEXTPAR = 10 # the number of pixels added from griddes to EXTPAR

ncBaseInformation = {
        'author': 'Niklas WAGNER, Klaus GOERGEN, Alexandre BELLEFLAMME',
        'contact': 'n.wagner@fz-juelich.de, k.goergen@fz-juelich.de, a.belleflamme@fz-juelich.de',
        'institution': 'FZJ/IBG-3',
        'history': f'Created {datetime.today().strftime("%d/%m/%y")}',
        'source': ''
        }
###############################################################################
# THE USER HAS TO CHANGE ABOVE ################################################
###############################################################################

# Read output grid
ncFileGridIN = Dataset(ncGridIN)
rlonIN = ncFileGridIN.variables['rlon'][:]
rlatIN = ncFileGridIN.variables['rlat'][:]
lonIN = ncFileGridIN.variables['lon'][:,:]
latIN = ncFileGridIN.variables['lat'][:,:]

nx = rlonIN.shape[0] # output grid dimension along x-axis (lon)
ny = rlatIN.shape[0] # output grid dimension along y-axis (lat)

# Read input land fraction and elevation
ncFileLsmIN = Dataset(ncLsmIN)
lsmIN = ncFileLsmIN.variables['FR_LAND'][nboundEXTPAR:-nboundEXTPAR, nboundEXTPAR:-nboundEXTPAR]
demIN = ncFileLsmIN.variables['HSURF'][nboundEXTPAR:-nboundEXTPAR, nboundEXTPAR:-nboundEXTPAR]

# Make and optimize land-lake-sea-mask, with land = 2, lake = 1, sea = 0
msk = np.full((ny,nx),2.)
msk = np.where((demIN==0.) & (lsmIN<0.5), 0., msk)
msk = np.where((demIN!=0.) & (lsmIN<0.5), 1., msk)


# optimisation loop 1 "from 0" forward
print('loop1')
for j in range(1,ny-1):
    for i in range(1,nx-1):
        if msk[j,i] == 0. :
            if msk[j+1,i] == 1. :
                n = j+1
                while n < ny and msk[n,i] == 1. :
                    msk[n,i] = 0.
                    n += 1
            if msk[j-1,i] == 1. :
                n = j-1
                while n >= 0 and msk[n,i] == 1. :
                    msk[n,i] = 0.
                    n -= 1
            if msk[j,i+1] == 1. :
                m = i+1
                while m < nx and msk[j,m] == 1. :
                    msk[j,m] = 0.
                    m += 1
            if msk[j,i-1] == 1. :
                m = i-1
                while m >= 0 and msk[j,m] == 1. :
                    msk[j,m] = 0.
                    m -= 1

# optimisation loop 1 "from 0" backward
print('loop2')
for j in range(ny-2,0,-1):
    for i in range(nx-2,0,-1):
        if msk[j,i] == 0. :
            if msk[j+1,i] == 1. :
                n = j+1
                while n < ny and msk[n,i] == 1. :
                    msk[n,i] = 0.
                    n += 1
            if msk[j-1,i] == 1. :
                n = j-1
                while n >= 0 and msk[n,i] == 1. :
                    msk[n,i] = 0.
                    n -= 1
            if msk[j,i+1] == 1. :
                m = i+1
                while m < nx and msk[j,m] == 1. :
                    msk[j,m] = 0.
                    m += 1
            if msk[j,i-1] == 1. :
                m = i-1
                while m >= 0 and msk[j,m] == 1. :
                    msk[j,m] = 0.
                    m -= 1

# optimisation loop 2 "from 1" forward
print('loop3')
for j in range(1,ny-1):
    for i in range(1,nx-1):
        if msk[j,i] == 1. :
            if msk[j-1,i] == 0. :
                n = j
                while n < ny and msk[n,i] == 1. :
                    msk[n,i] = 0.
                    n += 1
            if msk[j+1,i] == 0. :
                n = j
                while n >= 0 and msk[n,i] == 1. :
                    msk[n,i] = 0.
                    n -= 1
            if msk[j,i-1] == 0. :
                m = i
                while m < nx and msk[j,m] == 1. :
                    msk[j,m] = 0.
                    m += 1
            if msk[j,i+1] == 0. :
                m = i
                while m >= 0 and msk[j,m] == 1. :
                    msk[j,m] = 0.
                    m -= 1

# optimisation loop 2 "from 1" backward
print('loop4')
for j in range(ny-2,0,-1):
    for i in range(nx-2,0,-1):
        if msk[j,i] == 1. :
            if msk[j-1,i] == 0. :
                n = j
                while n < ny and msk[n,i] == 1. :
                    msk[n,i] = 0.
                    n += 1
            if msk[j+1,i] == 0. :
                n = j
                while n >= 0 and msk[n,i] == 1. :
                    msk[n,i] = 0.
                    n -= 1
            if msk[j,i-1] == 0. :
                m = i
                while m < nx and msk[j,m] == 1. :
                    msk[j,m] = 0.
                    m += 1
            if msk[j,i+1] == 0. :
                m = i
                while m >= 0 and msk[j,m] == 1. :
                    msk[j,m] = 0.
                    m -= 1

# Output netCDF file
print('write data')
ncout = Dataset(ncOUT, 'w', format = 'NETCDF4')

ncout.createDimension('rlon', nx)
ncout.createDimension('rlat', ny)
ncout.createDimension('time', None)

rlon = ncout.createVariable('rlon', 'f4', 'rlon')
rlat = ncout.createVariable('rlat', 'f4', 'rlat')
time = ncout.createVariable('time', 'i4', 'time')
glon = ncout.createVariable('lon', 'f4', ('rlat','rlon'))
glat = ncout.createVariable('lat', 'f4', ('rlat','rlon'))
llsm = ncout.createVariable('LLSM', 'f4', ('time','rlat','rlon'))
rpol = ncout.createVariable('rotated_pole', 'c')

rlon[:] = rlonIN
rlat[:] = rlatIN
glon[:,:] = lonIN
glat[:,:] = latIN
llsm[0,:,:] = msk

today = datetime.today()
time_num = today.toordinal()
time[0] = time_num

rlon.standard_name = 'grid_longitude'
rlat.standard_name = 'grid_latitude'
time.standard_name = 'Time'
glon.standard_name = 'longitude'
glat.standard_name = 'latitude'
llsm.standard_name = 'LLSM'

rlon.long_name = 'longitude in rotated pole grid'
rlat.long_name = 'latitude in rotated pole grid'
time.long_name = 'Time'
glon.long_name = 'geographical longitude'
glat.long_name = 'geographical latitude'
llsm.long_name = 'Land-lake-sea mask'

rlon.units = 'degrees'
rlat.units = 'degrees'
time.units = 'days since Jan 01, 0001'
glon.units = 'degrees_east'
glat.units = 'degrees_north'
llsm.units = '-'

llsm.grid_mapping = 'rotated_pole'
llsm.coordinates = 'lon lat'

rpol.grid_mapping_name = 'rotated_latitude_longitude'
rpol.grid_north_pole_longitude = -162.
rpol.grid_north_pole_latitude = 39.25
rpol.north_pole_grid_longitude = 0.

# Add basic information
ncout.author      = f'{ncBaseInformation["author"]}'
ncout.contact     = f'{ncBaseInformation["contact"]}'
ncout.institution = f'{ncBaseInformation["institution"]}'
ncout.history     = f'{ncBaseInformation["history"]}'
ncout.source      = f'{ncBaseInformation["source"]}'
ncout.description = 'Land=2 Lake=1 Sea=0 Mask'

ncout.close()

plt.imshow(msk, origin='lower', interpolation='none')
plt.savefig('LLSM.pdf')

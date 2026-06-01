#!/usr/bin/env python3
import sys
import numpy as np
#from netCDF4 import Dataset
from pysheds.grid import Grid

try:
    topo_tiff = sys.argv[1]
except IndexError:
    sys.exit(f"File name needed: '{sys.argv[0]} FILE'")

# Read elevation raster
grid = Grid.from_raster(topo_tiff)
dem = grid.read_raster(topo_tiff)

## Fill pits and depressions and resolve flats in DEM
pit_filled_dem = grid.fill_pits(dem)
flooded_dem = grid.fill_depressions(pit_filled_dem)
inflated_dem = grid.resolve_flats(flooded_dem)

# Write pit-filled to GeoTIFF
grid.to_raster(data=inflated_dem, file_name='pit-filled_dem.tiff')


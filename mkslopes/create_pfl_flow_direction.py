#!/usr/bin/env python3
# Calculation of flow direction grid from HSURFBurnedAndMod-Sphere.tiff
# based on https://mattbartos.com/pysheds/

from pysheds.grid import Grid

# Read elevation raster
grid = Grid.from_raster('HSURF-Sphere.tiff')
dem = grid.read_raster('HSURF-Sphere.tiff')

## Fill pits and depressions and resolve flats in DEM
#pit_filled_dem = grid.fill_pits(dem)
#flooded_dem = grid.fill_depressions(pit_filled_dem)
#inflated_dem = grid.resolve_flats(flooded_dem)

# Determine D8 flow directions from DEM
dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
#fdir = grid.flowdir(inflated_dem, dirmap=dirmap)
fdir = grid.flowdir(dem, dirmap=dirmap)

# Calculate the slopes
slopes = grid.cell_slopes(dem, fdir)

# Write flow direction and slope fields to GeoTIFF
grid.to_raster(data=fdir, file_name='flow_direction.tiff')
grid.to_raster(data=slopes, file_name='slopes.tiff')

# Determine accumulation and write to GeoTIFF
#acc = grid.accumulation(fdir, dirmap=dirmap)
#grid.to_raster(data=acc, file_name='accumulation.tiff')

# Determine river network
# N.B. 'branches' is not a Raster but a FeatureCollection containing tuples!
#branches = grid.extract_river_network(fdir, acc > 10, dirmap=dirmap)

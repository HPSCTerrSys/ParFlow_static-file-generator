#!/usr/bin/env python3
# Calculation of flow direction grid from HSURFBurnedAndMod-Sphere.tiff
# based on https://mattbartos.com/pysheds/

from pysheds.grid import Grid

# Read elevation raster
grid_orig = Grid.from_raster('HSURF-Sphere.tiff')
grid_burnt = Grid.from_raster('HSURFBurnedAndMod-Sphere.tiff')
dem_orig = grid_orig.read_raster('HSURF-Sphere.tiff')
dem_burnt = grid_burnt.read_raster('HSURFBurnedAndMod-Sphere.tiff')

#FIXME/TODO: NOT SURE IF THIS IS NEEDED FOR EITHER BURNT OR ORIG
## Fill pits and depressions and resolve flats in DEM
#pit_filled_dem = grid.fill_pits(dem)
#flooded_dem = grid.fill_depressions(pit_filled_dem)
#inflated_dem = grid.resolve_flats(flooded_dem)

# Determine D8 flow directions from DEM
dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
#fdir = grid.flowdir(inflated_dem, dirmap=dirmap)
fdir_orig = grid_orig.flowdir(dem_orig, dirmap=dirmap)
fdir_burnt = grid_burnt.flowdir(dem_burnt, dirmap=dirmap)

# Calculate the slopes
slopes = grid_burnt.cell_slopes(dem_orig, fdir_orig)
slopes = grid_burnt.cell_slopes(dem_orig, fdir_burnt)
#FIXME: In both cases most of values become zero

# Write flow direction and slope fields to GeoTIFF
grid_burnt.to_raster(data=fdir_burnt, file_name='flow_direction.tiff')
grid_orig.to_raster(data=slopes, file_name='slopes.tiff')

# Determine accumulation and write to GeoTIFF
#acc = grid.accumulation(fdir, dirmap=dirmap)
#grid.to_raster(data=acc, file_name='accumulation.tiff')

# Determine river network
# N.B. 'branches' is not a Raster but a FeatureCollection containing tuples!
#branches = grid.extract_river_network(fdir, acc > 10, dirmap=dirmap)

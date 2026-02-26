#!/usr/bin/env R
# Create D4 ParFlow slope files

topo_nc <- paste0(args[1])

# On JSC HPC 'remotes' is already included in the R module
if (!requireNamespace("remotes", quietly=TRUE)) install.packages('remotes')
library(remotes)
install_github("lecondon/PriorityFlow", subdir="Rpkg")

# Install packages that are dependencies but not described as such by PriorityFlow
install.packages('fields')

# Use the PriorityFlow package to calculate slopes
# TODO: complete and correct arguments!
library('PriorityFlow')
PriorityFlow::SlopeCalcUP(
            dem = topo_nc,
            direction = , # e.g. flow_direction.{tiff,nc} from pysheds routine
            dx = 0.11, dy = 0.11,
            mask = "../mklandmask/EUR-11_TSMP_FZJ-IBG3_444x432_LAND-LAKE-SEA-MASK.nc",
            borders = ,
            rivermask = # e.g. RiverMask.nc created by burnShape2Topo.py
)

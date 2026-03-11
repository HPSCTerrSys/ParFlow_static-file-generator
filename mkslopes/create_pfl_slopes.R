#!/usr/bin/env R
# Create D4 ParFlow slope files

if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    topo_nc <- args
    mask_nc <- args
}

# Install PriorityFlow
# On JSC HPC 'remotes' is already included in the R module
if (!requireNamespace("remotes", quietly=TRUE)) install.packages('remotes')
library(remotes)
install_github("lecondon/PriorityFlow", subdir="Rpkg", quiet=TRUE)

# Install packages that are dependencies but not described as such by PriorityFlow
if (!requireNamespace("fields", quietly=TRUE)) install.packages('fields', repos="https://ftp.fau.de/cran/")

# Use the PriorityFlow package to calculate slopes
library('PriorityFlow')

# netCDF support; part of R-bundle-CRAN/2025.11
library('terra')

llsm <- rast(mask_nc)
llsm[llsm==2] <- 1  # treat lakes as land
lsm_df <- data.frame(llsm)
lsm_matrix = data.matrix(lsm_df)
lsm_shaped = matrix(lsm_matrix, nrow=444)
hsurf <- rast(topo_nc)
hsurf_df <- data.frame(hsurf)
hsurf_matrix <- data.matrix(hsurf_df)
hsurf_shaped <- matrix(hsurf_matrix, nrow=444)
zero_matrix <- array( 0, dim(hsurf_shaped) )

# Calculate slopes; ParFlow needs SlopeCalcUP()
slope = PriorityFlow::SlopeCalStan( dem=hsurf_shaped, direction=zero_matrix,
                                    dx=12500, dy=12500, mask=lsm_shaped )
#image(slope$slopex)
#dev.off()

slopex_rast <- rast( slope$slopex )
slopey_rast <- rast( slope$slopey )
slope_dataset <- sds(slopex_rast, slopey_rast)
varnames(slope_dataset) <- c("slopex", "slopey")
writeCDF(slope_dataset, filename="slopes-out.nc")

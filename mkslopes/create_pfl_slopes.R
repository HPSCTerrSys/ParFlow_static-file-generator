#!/usr/bin/env R
# Create D4 ParFlow slope files

# The 'remotes' package is to grab stuff from GitHub &co, but at JSC it is in R
#install.packages('remotes')
library(remotes)
install_github("lecondon/PriorityFlow", subdir="Rpkg")

install.packages('fields')

library('PriorityFlow')

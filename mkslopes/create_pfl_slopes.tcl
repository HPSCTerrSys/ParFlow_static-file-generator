#!/usr/bin/env tclsh
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

pfload "HSURFBurnedAndMod.sa"
pfslopexD4 dataset0
pfsave dataset1 -sa "slopexD4.sa"
pfslopeyD4 dataset0
pfsave dataset2 -sa "slopeyD4.sa"

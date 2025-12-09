# ParFlow static-file generator

This repository shows how to generate land masks, slopes, solids and texture indicator for [ParFlow](https://github.com/parflow/parflow) simulations.
The workflow is taylored to generating files for running ParFlow in the [TSMP2](github.com/HPSCTerrSys/TSMP2) framework (coupled to eCLM), but it may be useful beyond that.
The generator incorporates steps from [the TSMP1 static file generator](https://gitlab.jsc.fz-juelich.de/detect/detect_z03_z04/constant_fields/TSMP_EUR-11), but it is updated, restructured and parts are refactored.
Instead of including raw input data in this repository, the repo is instead kept small.
The input data needed for this generator can be found [here](https://icg4geo.icg.kfa-juelich.de/ExternalReposPublic/tsmp2-static-files/grids_parflow_cordex-eur-11u).

If you are running this generator on a [JSC](https://www.fz-juelich.de/en/ias/jsc) machine, sourcing the provided environment file

```
source jsc.2025.gnu.psmpi
```

makes the necessary utilities and libraries available.
Otherwise you have to make sure that the respective software is installed or made available on your system.

To create all static files needed to run ParFlow, you first need to download static input files:

```
RAWDATA_DIR="grids_parflow_cordex-eur-11u"
git clone https://icg4geo.icg.kfa-juelich.de/ExternalReposPublic/tsmp2-static-files/${RAWDATA_DIR}.git
export RAWDATA_DIR="$(realpath $RAWDATA_DIR)"
```

Then create a land mask (`mklandmask/`), then the slopes (`mkslopes/`), the mask solids (`mksolids/`) and the texture indicator (`mktextureindicator/`) as described in the next sections.

## Creation of the land-lake-sea-mask

The purpose of this process is to create a land mask that distinguishes between land and water pixels.
This land mask ensures consistency and allows all model components in TSMP2 to see the same land and water areas.

To create the land mask, we use a variable called `FR_LAND` from `EXTPAR` files.
A threshold of 0.5 is applied, meaning that if a pixel is covered by less than 50% land, it is considered a pure water pixel.
Conversely, if a pixel is covered by more or equal 50% land, it is considered a pure land pixel.

Additionally, we differentiate between lake and sea pixels based on the topography.
Water pixels with a topographic height of zero or less are classified as sea pixels, while those with a topographic height greater than zero are classified as lake pixels.
However, this classification process introduces some artifacts along the coastline.
To address these artifacts, we run an optimization loop that checks if lake pixels are neighbors of sea pixels.
If they are, we treat them as sea pixels as well, ensuring a more accurate representation of the coastline.

It is important to note that COSMO and ICON also use the `FR_LAND` variable as a land-sea mask with also a threshold of 0.5.
This ensures that all components of TSMP share the same land-sea mask.

Furthermore, special treatment is required due to the coarse horizontal resolution of our target grid.
Some topographic formations, like the Bosporus connecting the Black Sea and the Mediterranean Sea, cannot be adequately resolved.
To address this, we modify the original `FR_LAND` variable by setting the value to 0, based on a predefined shape-files containing the coordinates of the Bosporus breakthrough, effectively designating those areas as total water pixels.

Above steps are performed by two scripts located in this directory:

```
cd mklandmask/
python3 Breakthrough-Bosporus.py
python3 make_land_lake_sea_mask.py
```

## Creation of the flow direction and slopes

Using DEMs straight forward to calculate slopes for ParFlow could lead to smaller or bigger problems in river corridor placement.
In particular for coarse spatial resolutions this is easy to imagine as tight canyons are smoothed out.
For example a 12km resolution as for the EU11 domain does not see the breakthrough valley [`irongate` for donau river](https://de.wikipedia.org/wiki/Eisernes_Tor) leading to the result, that the donau is flowing around the related mountain range.
A further example are the Netherlands, where major parts of the land area are below sea level and rivers do not follow the ‘natural’ river-corridor but are forced to follow artificial canals.
To fix those and other issues / problems the real river-corridors are ‘burned’ to the DEM within this approach.

Burning the correct river corridors is achieved in three steps:

1) `burnShape2Topo.py`  
The correct river positions are mapped to the target grid (in this case hydroSHEDS data were used) and then river pixels (and neighboring ones) are pushed down within the DEM.

2) `modTopo.py`  
Some pixels do need extra treatment, as for example the Elbe river in this setup.
Those need manual adjustment and are corrected ‘pixel-by-pixel’.

3) `sva_static_pfl.ncl`  
GRASS algorithmus is taken to calculate flow direction and main river streams based on the burned DEM.
To keep the correct slope values, those are calculated based on the original DEM, but the flow direction, represented by the sign of the slope value, is taken from flow direction calculated by GRASS algorithm.

This way the slope values are in line with the origin DEM, but flow direction is according to the correct river-corridors.

#### Usage
First extract and adjust the HydroSHEDS river data:

```
cd ../mkslopes
wget https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip
wget https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_af_shp.zip
unzip -u HydroRIVERS_v10_eu_shp.zip
unzip -u HydroRIVERS_v10_af_shp.zip
./burnShape2Topo.py
./modTopo.py
```

We use [pysheds](https://mattbartos.com/pysheds/) to calculate slopes and flow directions.
The following script will install and run pysheds.
pysheds uses XArray as an internal format and reads and writes to disk in GeoTIFF format.
We wrap around that, converting between netCDF and GeoTIFF.

```
./create_pfl_slopes
```

## Creation of the mask solids mask

## Creation of the texture indicator

----

This code is under development and its output files are not validated.
Please use [the prepared ParFlow static files](https://icg4geo.icg.kfa-juelich.de/ExternalReposPublic/tsmp2-static-files/extpar_parflow_cordex-eur-11u) if you just want to run ParFlow on a 12-km CORDEX domain.

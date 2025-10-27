# ParFlow static-file generator

This repository shows how to generate land masks, slopes, solids and texture indicator for [ParFlow](https://github.com/parflow/parflow) simulations.
The workflow is taylored to generating files for running ParFlow in the [TSMP2](github.com/HPSCTerrSys/TSMP2) framework (coupled to eCLM), but it may be useful beyond that.
The generator incorporates steps from [the TSMP1 static file generator](https://gitlab.jsc.fz-juelich.de/detect/detect_z03_z04/constant_fields/TSMP_EUR-11), but it is updated, restructured and rewritten.

If you are running this generator on a [JSC](https://www.fz-juelich.de/en/ias/jsc) machine, sourcing the provided environment file

```
source jsc.2025.gnu.psmpi
```

makes the necessary utilities and libraries available.
Otherwise you have to make sure that the respective software is installed or made available on your system.

To create all static files needed to run ParFlow, you first need to create a land mask (`mklandmask/`), then the slopes (`mkslopes/`), the mask solids (`mksolids/`) and the texture indicator (`mktextureindicator/`).

## Creation of the land-lake-sea-mask

## Creation of the flow direction and slopes

## Creation of the mask solids mask

## Creation of the texture indicator


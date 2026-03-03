# Breach.jl
Detect levees and other water-retaining barriers in DEMs by identifying depressions and their watersheds, as described in Pronk, M. et. al. (2026) Automated Levee Detection in Digital Elevation Models.

## Reproduce
Check `example/data/README.md` to download the required example 15 m DEM data for Zeeland in the Netherlands, including reference data.

Please [install Julia ](https://julialang.org/downloads/), and run the following in this directory:

*installs all required packages*

```julia --project=example -e "using Pkg; Pkg.instantiate()"```

*runs the example*

```julia --project=example example/breach26.jl```

This should result in output including validation statistics and a `example/detected_levees.gpkg` geopackage that you can visualize in (Q)GIS with several layers.

![QGIS](<example.png>)

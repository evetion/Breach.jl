using Breach
using GeoArrays
using Rasters
using GeoDataFrames
using Downloads

dir = @__DIR__
datadir = joinpath(@__DIR__, "data")

demfn = joinpath(datadir, "nl.tif")
isfile(demfn) || Downloads.download("https://github.com/evetion/Breach.jl/releases/download/v0.0.1/nl.tif", demfn)

# Read in a DEM and replacing nodata
dem = GeoArrays.read(demfn)
dem = coalesce(dem, Inf)
dem[isinf.(dem)] .= minimum(dem)

# Calculate Dataframe with objects. The levee is the part of each watershed outline connected to the spill point
# with a positive TPI (maximum over radii of 30, 90 and 150 m). In mountains, the levee can also be limited to a
# height above the spill height, e.g. `Breach.breach(dem; max_height=10)`.
data = Breach.breach(dem)

# Save them to disk for inspection in (Q)GIS
Breach.savegpkg("$dir/detected_levees", data, dem)

# Read in validation data for NL
valfn = joinpath(datadir, "validation_levees.gpkg")
isfile(valfn) || Downloads.download("https://github.com/evetion/Breach.jl/releases/download/v0.0.1/validation_levees.gpkg", valfn)
df = GeoDataFrames.read(valfn)

# Rasterize the validation data, and the predicted data
r = Raster(demfn; lazy=true)
actual = rasterize(df, fill=x -> true, missingval=false, to=r, boundary=:touches)
pred = rasterize(data.watershed, fill=x -> true, missingval=false, to=r, boundary=:touches)

# Compute statistics, for the watershed outlines and for the levee parts of them
@info "Statistics (watersheds)" Breach.stats(pred, actual)
levees = rasterize(collect(skipmissing(data.levee)), fill=x -> true, missingval=false, to=r, boundary=:touches)
@info "Statistics (levees)" Breach.stats(levees, actual)

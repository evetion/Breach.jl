using Breach
using GeoArrays
using Rasters
using GeoDataFrames

dir = @__DIR__
datadir = joinpath(@__DIR__, "data")

# Read in a DEM and replacing nodata
dem = GeoArrays.read(joinpath(datadir, "nl.tif"))
dem = coalesce(dem, Inf)
dem[isinf.(dem)] .= minimum(dem)

# Calculate Dataframe with objects
data = Breach.breach(dem)

# Save them to disk for inspection in (Q)GIS
Breach.savegpkg("$dir/detected_levees", data, dem)

# Read in validation data for NL
fn = joinpath(datadir, "validation_levees.gpkg")
df = GeoDataFrames.read(fn)

# Rasterize the validation data, and the predicted data
r = Raster(joinpath(datadir, "nl.tif"); lazy=true)
actual = rasterize(df, fill=x -> true, missingval=false, to=r, boundary=:touches)
pred = rasterize(data.watershed, fill=x -> true, missingval=false, to=r, boundary=:touches)

# Compute statistics
@info "Statistics" Breach.stats(pred, actual)

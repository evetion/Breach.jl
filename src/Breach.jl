module Breach

using DataStructures
using ProgressMeter
import LocalFilters
import Contour
using GeoArrays
import GeoInterface
import StaticArraysCore
import Geomorphometry
using StaticArrays
using OffsetArrays
using DataFrames
using GeoDataFrames
using Statistics
import GeometryOps as GO

function edges(A::AbstractMatrix)
    CI = CartesianIndices(A)
    edges(CI)
end

function edges(CI::CartesianIndices)
    indices = Vector{CartesianIndex}()
    append!(indices, first(eachrow(CI)))
    append!(indices, last(eachrow(CI)))
    append!(indices, first(eachcol(CI)))
    append!(indices, last(eachcol(CI)))
    unique(indices)
end

const ha = 10_000.0  # m²

const nbs =
    CartesianIndex.([(-1, -1), (-1, 1), (1, -1), (1, 1), (-1, 0), (0, -1), (0, 1), (1, 0)])


function breach(dem; area_threshold=1ha)
    depressions, labels, fdem, dir, outline, acc = Breach.sbreach(dem)

    dem_cellarea = abs(prod(Geomorphometry.cellsize(dem)))
    dem_celllength = sqrt(dem_cellarea)

    filter!(x -> x.second.area > dem_cellarea, depressions)
    filter!(x -> x.second.volume > 0, depressions)

    data = DataFrame(merge(v, (; label=k)) for (k, v) in depressions)

    data[!, :max_depth] = data.spillheight .- data.min  # 
    data[!, :avg_depth] = data.volume ./ data.area
    data[!, :slope] = max.(data.slopein, data.slopeout)

    filter!(x -> isfinite(x.area), data)
    filter!(x -> x.area .> (area_threshold), data)

    data[!, :point] = GeoArrays.coords.(Ref(dem), Tuple.(data.spillpoint), Ref(GeoArrays.Center()))
    data[!, :geometry] = Breach._contours.(Ref(dem), data.spillpoint, data.depressionpoint, data.spillheight .- eps.(data.spillheight), data.ex)
    data[!, :entrypoint] = data.spillpoint .+ getindex.(Ref(dir), data.spillpoint)

    filter!(x -> dir[x.spillpoint] .!= CartesianIndex(0, 0), data)

    watersheds, levees, elevations, _mask, extents = Breach.watershed_contour(dir, data.spillpoint, data.entrypoint, dem, data.spillheight, data.label)
    data[!, :watershed] = watersheds
    data[!, :roundness] = (GeoInterface.npoint.(data.watershed) * dem_celllength) ./ data.area
    data[!, :levee] = limit_line.(data.watershed, elevations, data.spillheight, 10)

    return data
end

mad(x) = median(abs.(x .- median(x)))

function sbreach(dem::AbstractMatrix, queued=falses(size(dem)))
    fdem = deepcopy(dem)
    open = PriorityQueue{CartesianIndex{2},eltype(dem)}()
    pit = PriorityQueue{CartesianIndex{2},eltype(dem)}()

    labels = similar(dem, Int32)
    # -1 is watershed, 0 is inital, and everything else is a label
    fill!(labels, 0)
    label = 1

    cellsize = Geomorphometry.cellsize(dem)
    cellarea = abs(cellsize[1] * cellsize[2])

    w, h = cellsize
    d = sqrt(w^2 + h^2)
    w = abs(w)
    h = abs(h)
    dist = OffsetArrays.centered(SMatrix{3,3,Float64}([d h d; w 0 w; d h d]))

    acc = similar(dem, Float32)
    acc .= cellarea

    order = ones(Int64, length(dem))
    L = LinearIndices(dem)
    dir = fill(CartesianIndex{2}(0, 0), size(dem))

    breach = similar(dem, Bool)
    fill!(breach, false)

    outline = similar(dem, Bool)
    fill!(breach, false)

    breach = similar(dem, Bool)
    fill!(breach, false)

    i = 0
    # Administration
    nt = NamedTuple{(:min, :area, :volume, :spillheight, :spillpoint, :depressionpoint, :ex, :startex, :parent, :level, :acc, :slopein, :slopeout),Tuple{eltype(dem),eltype(dem),eltype(dem),eltype(dem),CartesianIndex{2},CartesianIndex{2},Tuple{CartesianIndex{2},CartesianIndex{2}},CartesianIndex{2},Int,Int,Float32,eltype(dem),eltype(dem)}}
    data = OrderedDict{Int,nt}()

    R = CartesianIndices(dem)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)

    # Start at the edges of the DEM
    @inbounds for cell in edges(R)
        enqueue!(open, cell, dem[cell])
        queued[cell] = true  # queued
    end
    @inbounds while !isempty(open) || !isempty(pit)
        # Either process the next lowest cell from open
        # or continue filling a depression from pit
        ispit = false
        i += 1
        cell = if !isempty(pit)
            ispit = true
            cell = dequeue!(pit)
        else
            cell = dequeue!(open)
        end

        order[i] = L[cell]

        # For all cells, we loop over all neighbours
        # that have not been visited yet
        # Each cell is found twice, once from the queue, once as neighbour
        # for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
        for nb in nbs
            ncell = cell + nb
            ncell in R || continue

            (queued[ncell] || ncell == cell) && continue
            queued[ncell] = true

            dir[ncell] = cell - ncell

            if fdem[ncell] < fdem[cell]
                # Depression behaviour
                if labels[cell] == 0 && !breach[cell]
                    # Found a new depression, cell is the breach
                    data[label] = (;
                        min=Inf,
                        area=0,
                        volume=0,
                        spillheight=dem[cell],
                        spillpoint=cell,
                        depressionpoint=ncell,
                        ex=(min(cell, ncell), max(cell, ncell)),
                        startex=ncell,
                        parent=0,
                        level=0,
                        acc=0,
                        slopein=(dem[cell] - dem[ncell]) / dist[dir[ncell]],
                        slopeout=(dem[cell] - dem[cell+dir[cell]]) / dist[dir[ncell]],
                    )
                    breach[cell] = true # mark for other neighbors
                    labels[ncell] = label  # breach
                    label += 1
                elseif breach[cell]
                    # Expanding the depression from breach
                    labels[ncell] = label - 1
                elseif labels[cell] > 0
                    # Expand current depression, but check for nesting
                    if (dem[cell] > data[labels[cell]].min) && (dem[ncell] < dem[cell]) && !breach[cell]
                        # Current cell ascends (data has not yet been updated), 
                        # but neighbor is descending again. We found a nested depression.
                        data[label] = (;
                            min=Inf,
                            area=0,
                            volume=0,
                            spillheight=dem[cell],
                            spillpoint=cell,
                            depressionpoint=ncell,
                            ex=(min(cell, ncell), max(cell, ncell)),
                            startex=ncell,
                            parent=labels[cell],
                            level=data[labels[cell]].level + 1,
                            acc=0,
                            slopein=(dem[cell] - dem[ncell]) / dist[dir[ncell]],
                            slopeout=(dem[cell] - dem[cell+dir[cell]]) / dist[dir[ncell]],
                        )
                        breach[cell] = true
                        labels[ncell] = label  # breach
                        label += 1
                    elseif dem[ncell] < dem[cell]
                        # Just descending
                        labels[ncell] = labels[cell]
                    elseif dem[ncell] >= data[labels[cell]].spillheight
                        # We're exiting a nested depression
                        labels[ncell] = data[labels[cell]].parent
                    else
                        # Normal depression ascension
                        labels[ncell] = labels[cell]
                    end
                else
                    error("Logic error")
                end
                # Fill depression to spill height
                fdem[ncell] = fdem[cell]
                # DataStructures.enqueue!(pit, ncell)
                enqueue!(pit, ncell, dem[ncell])
            else
                # Normal climbing behaviour
                if labels[cell] != 0
                    outline[ncell] = true
                end
                enqueue!(open, ncell, dem[ncell])
            end
        end
        # ispit || continue
        clabel = labels[cell]
        clabel <= 0 && continue  # Skip watershed cells
        celldata = data[clabel]
        # :min, :area, :volume, :spillheight, :spillpoint, :depressionpoint, :ex
        data[clabel] = (;
            min=min(celldata.min, dem[cell]),
            area=celldata.area + cellarea,
            volume=celldata.volume + (celldata.spillheight - dem[cell]) * cellarea,
            spillheight=celldata.spillheight,
            spillpoint=celldata.spillpoint,
            depressionpoint=celldata.depressionpoint,
            ex=_mergeex(celldata.ex, (cell, cell)),
            startex=isless(cell, celldata.startex) ? cell : celldata.startex,
            parent=celldata.parent,
            level=celldata.level,
            acc=celldata.acc,
            slopein=celldata.slopein,
            slopeout=celldata.slopeout
        )
    end

    for (key, celldata) in data
        data[key] = merge(celldata, (; acc=acc[celldata.spillpoint]))
    end

    # Update parents with data from children (e.g. fix volume)
    pairs = sort!(collect(data), by=x -> x[2].parent, rev=true)
    for (key, _) in pairs
        celldata = data[key]
        # data[key] = merge(celldata, (; acc=acc[celldata.depressionpoint]))
        # celldata = data[key]
        celldata.parent == 0 && continue
        parent = data[celldata.parent]
        data[celldata.parent] = (;
            min=min(parent.min, celldata.min),
            area=parent.area + celldata.area,
            volume=(parent.spillheight - celldata.spillheight) * celldata.area + parent.volume + celldata.volume,
            spillheight=parent.spillheight,
            spillpoint=parent.spillpoint,
            depressionpoint=parent.depressionpoint,
            ex=_mergeex(parent.ex, celldata.ex),
            startex=isless(parent.startex, celldata.startex) ? parent.startex : celldata.startex,
            parent=parent.parent,
            level=parent.level,
            acc=parent.acc,
            slopein=parent.slopein,
            slopeout=parent.slopeout
        )
    end

    return data, labels, fdem, dir, outline, acc
end


function watershed(A)
    out = similar(A, Int32)
    fill!(out, 0)

    R = CartesianIndices(A)
    Δ = CartesianIndex(ntuple(x -> 3 ÷ 2, ndims(A)))
    I_first, I_last = first(R), last(R)

    label = 1

    for I in sortperm(vec(A))
        cell = R[I]
        clabel = out[cell]
        @assert clabel == 0  # Unlabeled cell
        for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
            nlabel = out[ncell]
            # Mark cell with first found label
            if nlabel > 0 && clabel == 0
                clabel = nlabel
                # If a different label is found mark as watershed
            elseif nlabel > 0 && nlabel != clabel
                clabel = -1
            end
        end
        if clabel == 0
            clabel = label
            label += 1
        end
        out[cell] = clabel
    end
    out
end



function _mergeex(a::Tuple{CartesianIndex{2},CartesianIndex{2}}, b::Tuple{CartesianIndex{2},CartesianIndex{2}})
    (min(a[1], b[1]), max(a[2], b[2]))
end


_merge(a, b, breach, spill) = (
    min(a.min, b.min),
    a.area + b.area,
    a.min < b.min ?
    a.volume + b.volume + b.area * (b.min - a.min) :
    a.volume + b.volume + a.area * (a.min - b.min),
    breach,
    spill,
    _mergeex(a.ex, b.ex)
)

_finalize(a, breach, spill) = (
    a.min,
    a.area,
    # Inverse volume
    (a.breach - a.min) * a.area - a.volume,
    breach,
    spill,
    a.ex
)

function _contours(dem::AbstractArray, loc::CartesianIndex{2}, nb::CartesianIndex{2}, h::Number, ex)
    xi, yi = min.(loc.I, nb.I)
    Δ = CartesianIndex(1, 1)

    R = CartesianIndices(dem)
    I_first, I_last = first(R), last(R)
    offset = 1
    ex = (
        (max(I_first[1], ex[1][1] - offset), max(I_first[2], ex[1][2] - offset)),
        (min(I_last[1], ex[2][1] + offset), min(I_last[2], ex[2][2] + offset))
    )

    demv = view(dem, ex[1][1]:ex[2][1], ex[1][2]:ex[2][2])
    # Adjust starting coordinates to clipped DEM coordinate system
    xi = xi - ex[1][1] + 1
    yi = yi - ex[1][2] + 1
    c = Contour.contour(axes(demv)..., demv, h, start=(xi, yi))
    ll = Contour.lines(c)
    @assert length(ll) == 1
    l = first(ll)
    if isempty(ll)
        @warn "Empty line"
        GeoInterface.Wrappers.LineString([StaticArraysCore.SVector((0, 0)), StaticArraysCore.SVector(0, 0)])
    else
        GeoInterface.Wrappers.Polygon(GeoInterface.Wrappers.LineString(dem.f.([x .+ ex[1] .- 1.5 for x in Contour.vertices(l)])))
    end
end

"""
Create watershed mask and extent from starting cell.

Uses the flow direction matrix `dir` (ldd) to trace the watershed.
"""
function _watershed!(mask, levees, queue, dir::AbstractMatrix, startcell::CartesianIndex{2}, dem::AbstractMatrix, threshold=Inf, id=0)
    R = CartesianIndices(mask)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)

    extent = (startcell, startcell)

    enqueue!(queue, startcell)
    mask[startcell] = true
    while !isempty(queue)
        cell = DataStructures.dequeue!(queue)
        for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
            # Always mark levees as the outline of the depression
            if !mask[ncell] && cell == dir[ncell] + ncell
                if dem[cell] <= threshold && dem[ncell] > threshold
                    # If we exit the depression, mark as levee
                    # levees[ncell] = id
                elseif levees[cell] == id && dem[cell] >= threshold && dem[ncell] >= threshold
                    # If we flow beyond the levee, remove it
                    # levees[cell] = 0
                end
                mask[ncell] = true
                extent = (min(extent[1], ncell), max(extent[2], ncell))
                enqueue!(queue, ncell)
            end
        end
    end
    # Expand watershed boundary in case of higher neighbors.
    # Rasters have no real crest with two drainage directions, 
    # but two cells with opposite directions.
    i = 0
    startextent = extent[1]
    v = view(mask, extent[1]:extent[2])
    for vcell in findall(v)
        cell = startextent + vcell - CartesianIndex(1, 1)
        for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
            mask[ncell] && continue
            levees[cell] = id
            if (dem[ncell] > dem[cell]) && (cell != startcell)
                mask[ncell] = true
                if levees[cell] == id
                    levees[ncell] = id
                end
                extent = (min(extent[1], ncell), max(extent[2], ncell))
            end
        end
    end

    levees[startcell] = id
    mask, extent
end

function watershed_contour(dir, loc::Vector{CartesianIndex{2}}, nb::Vector{CartesianIndex{2}}, dem, threshold, ids)

    # Watershed areas
    mask = similar(dem, Bool)
    fill!(mask, false)

    # Levee areas
    levees = similar(dem, Int32)
    fill!(levees, 0)

    queue = DataStructures.Queue{CartesianIndex{2}}()

    output = Vector{Union{Missing,GeoInterface.Wrappers.LineString}}()
    x, y = axes(mask)
    extents = Vector{Tuple{CartesianIndex{2},CartesianIndex{2}}}()
    elevations = Vector{Vector{eltype(dem)}}()

    R = CartesianIndices(mask)
    I_first, I_last = first(R), last(R)

    @showprogress "Drawing watershed boundaries" for (loci, nbi, ti, id) in zip(loc, nb, threshold, ids)
        _, extent = _watershed!(mask, levees, queue, dir, loci, dem, ti, id)
        xi, yi = min.(loci.I, nbi.I)


        c = Contour.contour(x, y, mask, 0.99, start=(xi, yi), lazy=true)
        fill!(view(mask, extent[1]:extent[2]), false)
        ll = Contour.lines(c)
        if isempty(ll)
            @warn "Empty line"
            push!(output, missing)
            push!(elevations, eltype(dem)[])
            push!(extents, extent)
        else
            @assert length(ll) == 1
            l = first(ll)
            push!(output,
                GeoInterface.Wrappers.LineString(dem.f.([x .- 0.5 for x in Contour.vertices(l)]))
            )
            push!(elevations, [dem[round(Int, i), round(Int, j)] for (i, j) in Contour.vertices(l)])
            push!(extents, extent)
        end
        isempty(queue) || error("Queue not empty after watershed")
    end
    output, levees, elevations, mask, extents
end

const order = [
    CartesianIndex(-1, -1),
    CartesianIndex(-1, 0),
    CartesianIndex(-1, 1),
    CartesianIndex(0, 1),
    CartesianIndex(1, 1),
    CartesianIndex(1, 0),
    CartesianIndex(1, -1),
    CartesianIndex(0, -1),
    CartesianIndex(-1, -1),
]

GeoInterface.isgeometry(::Type{GeoInterface.Wrappers.LineString}) = true


function limit_line(line, elevations, spillheight, threshold=3)
    firstind = findfirst(x -> x > (threshold), elevations)
    lastind = findlast(x -> x > (threshold), elevations)

    isnothing(firstind) && return line

    coords = collect(GeoInterface.getpoint(line))
    newcoords = vcat(coords[lastind+1:end], coords[begin:firstind-1])
    if length(newcoords) < 2
        return GeoInterface.Wrappers.LineString([(0, 0), (0, 0)])
    end
    return GeoInterface.Wrappers.LineString(newcoords)
end


function savegpkg(fname, data, dem)
    sort!(data, :level)

    data[!, :point] = coords.(Ref(dem), Tuple.(data.spillpoint), Ref(GeoArrays.Center()))
    fdata = select(data, Not(:levee, :spillpoint, :depressionpoint, :entrypoint, :ex, :geometry, :startex, :watershed))
    GeoDataFrames.write("$fname.gpkg", fdata; layer_name="spillpoint", geom_columns=(:point,), crs=EPSG(4326), options=Dict("FID" => "label", "OVERWRITE" => "YES"))

    fdata = select(data, Not(:levee, :spillpoint, :depressionpoint, :entrypoint, :ex, :point, :startex, :watershed))
    GeoDataFrames.write("$fname.gpkg", fdata; update=true, layer_name="depression", geom_columns=(:geometry,), crs=EPSG(4326), options=Dict("FID" => "label"))

    fdata = select(data, Not(:levee, :spillpoint, :depressionpoint, :entrypoint, :ex, :point, :startex, :geometry))
    GeoDataFrames.write("$fname.gpkg", fdata; update=true, layer_name="watershed", geom_columns=(:watershed,), crs=EPSG(4326), options=Dict("FID" => "label"))

    fdata = select(data, Not(:spillpoint, :depressionpoint, :entrypoint, :ex, :point, :startex, :geometry, :watershed))
    GeoDataFrames.write("$fname.gpkg", fdata; update=true, layer_name="levee", geom_columns=(:levee,), crs=EPSG(4326), options=Dict("FID" => "label"))

end

function makepolygon(ex::Tuple{CartesianIndex{2},CartesianIndex{2}}, dem)
    (minx, miny), (maxx, maxy) = Tuple.(ex)

    GeoInterface.Wrappers.Polygon(
        GeoInterface.Wrappers.LineString([
            dem.f((minx - 0.5, miny - 0.5)),
            dem.f((minx - 0.5, maxy - 0.5)),
            dem.f((maxx - 0.5, maxy - 0.5)),
            dem.f((maxx - 0.5, miny - 0.5)),
            dem.f((minx - 0.5, miny - 0.5)),
        ]),
    )
end

function stats(pred, actual)
    TP = pred .& actual
    FP = pred .& .!actual
    FN = .!pred .& actual
    TN = .!pred .& .!actual
    P = sum(actual)
    TPR = sum(TP) / P
    IoU = sum(TP) / (sum(TP) + sum(FP) + sum(FN))
    f1 = 2 * sum(TP) / (2sum(TP) + sum(FP) + sum(FN))
    prec = sum(TP) / (sum(TP) + sum(FP))
    rec = TPR
    return (; prec, rec, f1)
end

end  # module Breach

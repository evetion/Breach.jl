module Breach

using DataStructures
using ProgressMeter
import LocalFilters
import Contour
using GeoArrays
import GeoInterface
import Geomorphometry
import QuickHeaps
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

# Orthogonal neighbors first, so that on equal slopes (e.g. flats)
# the steepest descent in `sbreach` prefers an orthogonal direction.
const nbs =
    CartesianIndex.([(-1, 0), (0, -1), (0, 1), (1, 0), (-1, -1), (-1, 1), (1, -1), (1, 1)])


"""
    breach(dem; area_threshold=1ha, tpi_radii=(30, 90, 150), max_height=nothing)

Detect depressions in `dem` (a GeoArray, without NaN), their watersheds, and the levee part of each watershed outline.
Returns a DataFrame with a row per depression of at least `area_threshold` (m²).

The levee (`:levee`) is the part of the watershed outline connected to the spill point with a positive topographic
position index (TPI): the maximum TPI over windows with radii `tpi_radii` (m), see `multiscale_tpi`.
With `max_height` (m), the levee is also limited to the part at most `max_height` above the spill height, e.g. in
mountains, where the watershed of a depression behind a levee can follow the ridges far above it.
"""
function breach(dem; area_threshold=1ha, tpi_radii=(30, 90, 150), max_height=nothing)
    depressions, labels, fdem, dir, outline = Breach.sbreach(dem)

    dem_cellarea = abs(prod(Geomorphometry.cellsize(dem)))
    dem_celllength = sqrt(dem_cellarea)

    filter!(x -> x.second.area > dem_cellarea, depressions)
    filter!(x -> x.second.volume > 0, depressions)

    # From a typed vector, so the DataFrame has all columns, also without depressions
    T = Base.promote_op((v, k) -> merge(v, (; label=k)), _depression_type(eltype(dem)), Int)
    data = DataFrame(T[merge(v, (; label=k)) for (k, v) in depressions])

    data[!, :max_depth] = data.spillheight .- data.min  #
    data[!, :avg_depth] = data.volume ./ data.area
    data[!, :slope] = max.(data.slopein, data.slopeout)

    filter!(x -> isfinite(x.area), data)
    filter!(x -> x.area .> (area_threshold), data)

    data[!, :point] = GeoArrays.coords.(Ref(dem), Tuple.(data.spillpoint), Ref(GeoArrays.Center()))
    data[!, :geometry] = Breach.depression_contours(dem, data.depressionpoint, data.spillheight)
    data[!, :entrypoint] = data.spillpoint .+ getindex.(Ref(dir), data.spillpoint)

    filter!(x -> dir[x.spillpoint] .!= CartesianIndex(0, 0), data)
    isempty(data) && @warn "No depressions with an area above $(area_threshold) m² found, returning an empty DataFrame"

    watersheds, levees, elevations, nearest, extents = Breach.watershed_contour(dir, data.spillpoint, data.entrypoint, dem, data.spillheight, data.label)
    data[!, :watershed] = watersheds
    data[!, :roundness] = (GeoInterface.npoint.(data.watershed) * dem_celllength) ./ data.area

    # The TPI at the vertices of the watershed outlines, which lie on cell centres
    tpi = isempty(data) ? zeros(Float32, 0, 0) : multiscale_tpi(dem, tpi_radii)
    tpis = [ismissing(line) ? eltype(tpi)[] : tpi[_cells(line, dem)] for line in watersheds]
    data[!, :levee] = levee_line.(data.watershed, tpis, elevations, data.spillheight; max_height)

    return data
end

mad(x) = median(abs.(x .- median(x)))

"""The administration of a depression in `sbreach`, for a DEM with elevations of type `T`."""
_depression_type(T) = NamedTuple{(:min, :area, :volume, :spillheight, :spillpoint, :depressionpoint, :ex, :startex, :parent, :level, :slopein, :slopeout),
    Tuple{T,T,T,T,CartesianIndex{2},CartesianIndex{2},Tuple{CartesianIndex{2},CartesianIndex{2}},CartesianIndex{2},Int,Int,T,T}}

function sbreach(dem::AbstractMatrix, queued=falses(size(dem)))
    # NaN breaks the ordering of the queues (see below), replace them beforehand.
    any(isnan, dem) && throw(ArgumentError("dem contains NaN values, replace them before calling sbreach"))
    fdem = deepcopy(dem)
    # Priority is the elevation, with ties broken by insertion order (first in, first out).
    # Without it, cells of equal elevation (flats) are dequeued in arbitrary order,
    # and a cell can be processed before its orthogonal neighbors, which then can't
    # be chosen as steepest descent. With it, flats are processed breadth-first
    # from where they were found, like the edges before the interior.
    # A cell is queued only once and its priority never changes, so a plain binary heap
    # of (elevation, counter, linear index) suffices. A priority queue (as in DataStructures)
    # also maintains a lookup from cell to heap position, which made it twice as slow.
    # Note that FastMin (like <) is not consistent with NaN, so the DEM should not contain NaN.
    open = QuickHeaps.FastBinaryHeap{Tuple{eltype(dem),Int,Int}}(QuickHeaps.FastMin)
    pit = QuickHeaps.FastBinaryHeap{Tuple{eltype(dem),Int,Int}}(QuickHeaps.FastMin)
    n = 0  # insertion counter
    L = LinearIndices(dem)
    R = CartesianIndices(dem)

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
    # cellsize is (x, y) and x is the first dimension (as in GeoArrays),
    # so a step in the first dimension (a row in this matrix) has length w.
    # Only matters for non-square cells, like 1" DEMs away from the equator.
    dist = OffsetArrays.centered(SMatrix{3,3,Float64}([d w d; h 0 h; d w d]))

    dir = fill(CartesianIndex{2}(0, 0), size(dem))

    breach = similar(dem, Bool)
    fill!(breach, false)

    # Cells that have been dequeued, and thus have a path to the edge
    processed = similar(dem, Bool)
    fill!(processed, false)

    outline = similar(dem, Bool)
    fill!(outline, false)

    i = 0
    # Administration
    nt = _depression_type(eltype(dem))
    data = OrderedDict{Int,nt}()

    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)

    # Start at the edges of the DEM
    @inbounds for cell in edges(R)
        push!(open, (dem[cell], n += 1, L[cell]))
        queued[cell] = true  # queued
    end
    @inbounds while !isempty(open) || !isempty(pit)
        # Either process the next lowest cell from open
        # or continue filling a depression from pit
        ispit = false
        i += 1
        cell = if !isempty(pit)
            ispit = true
            cell = R[last(pop!(pit))]
        else
            cell = R[last(pop!(open))]
        end

        processed[cell] = true

        # Flow direction by steepest descent.
        # On discovery (below), a cell points to the neighbor that found it, which
        # is the lowest neighbor processed at that time, not the steepest one.
        # A diagonal neighbor can be lowest, but less steep than an orthogonal one.
        # Now that the cell is dequeued, more neighbors have been processed, so we pick
        # the steepest descent among them. Pointing to any processed neighbor keeps the
        # flow directions acyclic and in `order`.
        # Cells in a depression (pit) keep the direction of discovery, which follows the
        # filling from the spill point, so their flow is consistent with the (nested) labels.
        # Edge cells (no direction) remain outlets.
        if !ispit && dir[cell] != CartesianIndex(0, 0)
            steepest = -Inf
            for nb in nbs
                ncell = cell + nb
                ncell in R || continue
                processed[ncell] || continue
                slope = (dem[cell] - dem[ncell]) / dist[nb]
                # Strictly steeper, so earlier (orthogonal) neighbors win ties
                if slope > steepest
                    steepest = slope
                    dir[cell] = nb
                end
            end
        end

        # For all cells, we loop over all neighbours
        # that have not been visited yet
        # Each cell is found twice, once from the queue, once as neighbour
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
                        slopein=(dem[cell] - dem[ncell]) / dist[dir[ncell]],
                        slopeout=(dem[cell] - dem[cell+dir[cell]]) / dist[dir[cell]],
                    )
                    breach[cell] = true # mark for other neighbors
                    labels[ncell] = label  # breach
                    label += 1
                elseif breach[cell] && (labels[cell] == 0 || dem[ncell] < dem[cell])
                    # Expanding the depression from breach.
                    # Each cell is dequeued only once, so the depression that
                    # was created from this breach cell is always the last one (label - 1).
                    # For a top-level breach (labels[cell] == 0), fdem[cell] == dem[cell],
                    # so every neighbor that got here is strictly lower and belongs to it.
                    # For a nested breach, fdem[cell] is the fill level of the parent,
                    # so neighbors above the nested spill height also arrive here. Those
                    # still belong to the parent, and fall through to the branch below.
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
                            slopein=(dem[cell] - dem[ncell]) / dist[dir[ncell]],
                            slopeout=(dem[cell] - dem[cell+dir[cell]]) / dist[dir[cell]],
                        )
                        breach[cell] = true
                        labels[ncell] = label  # breach
                        label += 1
                    elseif dem[ncell] < dem[cell]
                        # Just descending
                        labels[ncell] = labels[cell]
                    elseif dem[ncell] >= data[labels[cell]].spillheight
                        # We're exiting a nested depression.
                        # Nested spill heights increase towards the top-level depression,
                        # but a single step can exit more than one level at once,
                        # so walk up the parents until ncell is below the spill height.
                        # This can't reach 0, as ncell is below the top-level fill level.
                        parent = data[labels[cell]].parent
                        while parent != 0 && dem[ncell] >= data[parent].spillheight
                            parent = data[parent].parent
                        end
                        labels[ncell] = parent
                    else
                        # Normal depression ascension
                        labels[ncell] = labels[cell]
                    end
                else
                    error("Logic error")
                end
                # Fill depression to spill height
                fdem[ncell] = fdem[cell]
                push!(pit, (dem[ncell], n += 1, L[ncell]))
            else
                # Normal climbing behaviour
                if labels[cell] != 0
                    outline[ncell] = true
                end
                push!(open, (dem[ncell], n += 1, L[ncell]))
            end
        end
        clabel = labels[cell]
        clabel <= 0 && continue  # Skip watershed cells
        celldata = data[clabel]
        # :min, :area, :volume, :spillheight, :spillpoint, :depressionpoint, :ex
        data[clabel] = (;
            min=min(celldata.min, dem[cell]),
            area=celldata.area + cellarea,
            # A cell of a parent depression can lie above the spill height of the depression
            # it was found from (see the breach branch above), so no negative volumes.
            volume=celldata.volume + max(0, (celldata.spillheight - dem[cell])) * cellarea,
            spillheight=celldata.spillheight,
            spillpoint=celldata.spillpoint,
            depressionpoint=celldata.depressionpoint,
            ex=_mergeex(celldata.ex, (cell, cell)),
            startex=isless(cell, celldata.startex) ? cell : celldata.startex,
            parent=celldata.parent,
            level=celldata.level,
            slopein=celldata.slopein,
            slopeout=celldata.slopeout
        )
    end

    # Update parents with data from children (e.g. fix volume)
    pairs = sort!(collect(data), by=x -> x[2].parent, rev=true)
    for (key, _) in pairs
        celldata = data[key]
        celldata.parent == 0 && continue
        parent = data[celldata.parent]
        data[celldata.parent] = (;
            min=min(parent.min, celldata.min),
            area=parent.area + celldata.area,
            volume=max(0, parent.spillheight - celldata.spillheight) * celldata.area + parent.volume + celldata.volume,
            spillheight=parent.spillheight,
            spillpoint=parent.spillpoint,
            depressionpoint=parent.depressionpoint,
            ex=_mergeex(parent.ex, celldata.ex),
            startex=isless(parent.startex, celldata.startex) ? parent.startex : celldata.startex,
            parent=parent.parent,
            level=parent.level,
            slopein=parent.slopein,
            slopeout=parent.slopeout
        )
    end

    return data, labels, fdem, dir, outline
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

"""
Start (lower corner of a 2x2 block of cells) for a contour on the outer outline of a set of cells,
given its first cell in column-major order.

Starting next to the spill point is not enough: the 2x2 block there can be a saddle
(two inside cells on one diagonal, two outside cells on the other) or touch a hole
(an outside cell enclosed by the set), and Contour then traces a small ring
around a single cell instead of the outline.
For the first cell (i, j), the previous column (j - 1) has no cells of the set, so it
is outside and connected to the exterior. The block with corners (i - 1, j - 1) and (i, j)
then has a single inside corner, so it's never a saddle, and its line is on the outer outline.
Edge cells are outlets and never part of a watershed or depression, so i, j > 1.
"""
_outer_start(firstcell::CartesianIndex{2}) = (firstcell - CartesianIndex(1, 1)).I

"""
    depression_contours(dem, points, spills)

Outlines of depressions, filled to their spill height, as polygons (without holes).
For each depression its depression point (`points`, e.g. `depressionpoint` of `sbreach`) and spill height (`spills`).

Replaces `_contours`, which traced the DEM contour at the spill height from a start next to the spill point,
and took the first line. That line is often only a small part of the depression: a contour on the DEM
treats a depression of which the parts only connect diagonally as separate rings, and the start can
touch other rings (the same problem as for the watershed outlines, see `_outer_start`).

Here the topology comes from the depression cells, and the geometry from the DEM:
- the depression consists of the cells below the spill height that are connected (8 neighbours) to its
  depression point. In a priority flood these are exactly the cells of the depression and its nested
  depressions, as a lower area that connects to it below the spill height would have been part of the same flood,
- its outline is traced on these cells (contour at 0.5), starting at its first cell (see `_outer_start`),
  so it is always the outer outline,
- each vertex lies on the edge between a cell inside and a cell outside the depression, and is moved along
  that edge to where the DEM (linearly interpolated) reaches the spill height, like a contour of the DEM.
"""
function depression_contours(dem, points::AbstractVector{CartesianIndex{2}}, spills::AbstractVector)
    A = parent(dem)
    n = length(points)
    out = Vector{Union{Missing,GeoInterface.Wrappers.Polygon}}(missing, n)
    # Each task has its own work arrays, the mask is as large as the DEM
    chunks = collect(Iterators.partition(1:n, max(1, cld(n, Threads.nthreads()))))
    tasks = map(chunks) do idx
        Threads.@spawn begin
            mask, stack, touched = falses(size(A)), CartesianIndex{2}[], CartesianIndex{2}[]
            for k in idx
                out[k] = _depression_contour!(mask, stack, touched, dem, A, points[k], spills[k])
            end
        end
    end
    foreach(wait, tasks)
    out
end

function _depression_contour!(mask, stack, touched, dem, A, point, spill)
    R = CartesianIndices(A)
    A[point] < spill || return missing
    # The depression cells: below the spill height, connected to the depression point
    firstcell = point
    push!(stack, point)
    mask[point] = true
    push!(touched, point)
    while !isempty(stack)
        cell = pop!(stack)
        # The first cell in column-major order (note that `min` of CartesianIndex is elementwise)
        isless(cell, firstcell) && (firstcell = cell)
        for nb in nbs
            j = cell + nb
            (j in R && !mask[j] && A[j] < spill) || continue
            mask[j] = true
            push!(touched, j)
            push!(stack, j)
        end
    end
    polygon = missing
    # A depression doesn't reach the edge of the DEM (its cells drain to the edge otherwise)
    if firstcell[1] > first(R)[1] && firstcell[2] > first(R)[2]
        x, y = axes(A)
        ll = Contour.lines(Contour.contour(x, y, mask, 0.5, start=_outer_start(firstcell), lazy=true))
        if !isempty(ll)
            vertices = map(Contour.vertices(first(ll))) do (i, j)
                # Between a cell inside and a cell outside the depression
                a = CartesianIndex(floor(Int, i), floor(Int, j))
                b = CartesianIndex(ceil(Int, i), ceil(Int, j))
                inside, outside = mask[a] ? (a, b) : (b, a)
                zi, zo = A[inside], A[outside]
                t = zo > zi ? clamp((spill - zi) / (zo - zi), 0.0, 1.0) : 0.5
                p = Tuple(inside) .+ t .* (Tuple(outside) .- Tuple(inside))
                # Cell centres are at index - 0.5 in the coordinates of the geotransform
                dem.f(p .- 0.5)
            end
            polygon = GeoInterface.Wrappers.Polygon([GeoInterface.Wrappers.LinearRing(vertices)])
        end
    end
    for t in touched
        mask[t] = false
    end
    empty!(touched)
    empty!(stack)
    polygon
end

"""
Crest of the edge between two neighboring cells `a` and `b`, the highest of the two.

Rasters have no real crest with two drainage directions, only two cells with
opposite directions, each in a different watershed. The edge between them is
shared by both watersheds, and this function only depends on that edge (not on
the order of `a` and `b`), so both watersheds pick the same cell as crest.
On equal elevations, we pick the first cell in index order.
"""
_crest(dem, a::CartesianIndex{2}, b::CartesianIndex{2}) =
    dem[a] > dem[b] ? a : dem[b] > dem[a] ? b : min(a, b)

"""
Snap the vertices of a watershed contour to the crest cells.

A contour at level 0.5 on a watershed mask puts all vertices halfway between two
4-connected neighbors, one inside and one outside the watershed. This line is
shared with the neighboring watershed, but lies between cell centres.
We move each vertex to the crest cell of its edge (see `_crest`), so the vertices
lie on cell centres, while both watersheds still end up with the same cells.
"""
function _snap(vertices, dem)
    cells = CartesianIndex{2}[]
    for (i, j) in vertices
        # One of the coordinates is halfway, the other is a whole number
        a = CartesianIndex(floor(Int, i), floor(Int, j))
        b = CartesianIndex(ceil(Int, i), ceil(Int, j))
        cell = _crest(dem, a, b)
        # Neighboring edges often share the same crest cell
        (isempty(cells) || cells[end] != cell) && push!(cells, cell)
    end
    cells
end

"""
Rotate a closed line of cells, so it starts (and ends) at the cell closest to `cell`.

The contour starts at the first cell of the watershed (see `_outer_start`), but
further processing (e.g. `levee_line`) expects the line to start at the spill point.
"""
function _rotate(cells::Vector{CartesianIndex{2}}, cell::CartesianIndex{2})
    # Only closed lines, of which the first and last cell are the same
    (length(cells) > 2 && first(cells) == last(cells)) || return cells
    ring = cells[1:end-1]
    i = argmin(sum(abs2, (c - cell).I) for c in ring)
    ring = circshift(ring, 1 - i)
    push!(ring, first(ring))
end

"""
Index of the nearest downstream spill point (in `spills`) for each cell, 0 if none.

Follows the flow directions from each cell until a spill point, an outlet, or a cell
of which the answer is already known. All cells on that path get the same answer,
so each cell is only followed once, and this takes O(n) for the whole DEM.
"""
function _nearest_spill(dir, spills)
    nearest = fill(Int32(-1), size(dir))  # -1 is unknown
    for (i, spill) in enumerate(spills)
        nearest[spill] = i
    end
    path = CartesianIndex{2}[]
    for start in CartesianIndices(dir)
        nearest[start] >= 0 && continue
        cell = start
        while nearest[cell] < 0
            if dir[cell] == CartesianIndex(0, 0)
                # Outlet at the edge without passing a spill point
                nearest[cell] = 0
                break
            end
            push!(path, cell)
            cell += dir[cell]
        end
        found = nearest[cell]
        while !isempty(path)
            nearest[pop!(path)] = found
        end
    end
    nearest
end

"""
Hierarchy of the spill points, as the entry and exit times of a depth-first traversal.

The parent of a spill point is the nearest spill point downstream of it. A spill point `j`
is then downstream of `i` (or `i` itself), iff `tin[i] <= tin[j] <= tout[i]`.
"""
function _spilltree(dir, spills, nearest)
    n = length(spills)
    children = [Int32[] for _ in 1:n]
    roots = Int32[]
    parents = zeros(Int32, n)
    for (i, spill) in enumerate(spills)
        d = dir[spill]
        parent = d == CartesianIndex(0, 0) ? 0 : nearest[spill+d]
        parents[i] = parent
        parent == 0 ? push!(roots, i) : push!(children[parent], i)
    end
    tin = zeros(Int32, n)
    tout = zeros(Int32, n)
    t = Int32(0)
    # Depth-first without recursion, as the hierarchy can be deep
    stack = [(root, false) for root in reverse(roots)]
    while !isempty(stack)
        i, visited = pop!(stack)
        if visited
            tout[i] = t
        else
            t += 1
            tin[i] = t
            push!(stack, (i, true))
            for child in children[i]
                push!(stack, (child, false))
            end
        end
    end
    tin, tout, parents
end

"""
First cell (in column-major order) of each watershed, see `_outer_start`.

A watershed consists of the cells with its spill point as nearest spill point,
and the watersheds upstream of it, so we take the first cell per nearest spill point,
and then combine these from the upstream to the downstream watersheds.
"""
function _firstcells(nearest, tin, parents)
    n = length(tin)
    unset = CartesianIndex(typemax(Int), typemax(Int))
    firstcells = fill(unset, n)
    # Iteration is in column-major order, so the first cell found per spill point is the first cell
    for cell in CartesianIndices(nearest)
        s = nearest[cell]
        s > 0 || continue
        firstcells[s] == unset && (firstcells[s] = cell)
    end
    # Upstream watersheds have a larger entry time than the ones downstream of them.
    # Note that `min` of two CartesianIndex is elementwise (a corner of their bounding box),
    # while `isless` is in column-major order, like `findfirst`.
    for i in sortperm(tin; rev=true)
        p = parents[i]
        p == 0 && continue
        isless(firstcells[i], firstcells[p]) && (firstcells[p] = firstcells[i])
    end
    firstcells
end

"""
Watershed mask of spill point `k`, evaluated lazily for each cell.

A cell is in the watershed of `k`, if its flow passes the spill point of `k`, so if
its nearest downstream spill point is `k` or upstream of `k` (see `_spilltree`).
"""
struct WatershedMask <: AbstractMatrix{Bool}
    nearest::Matrix{Int32}
    tin::Vector{Int32}
    tout::Vector{Int32}
    k::Int
end
Base.size(m::WatershedMask) = size(m.nearest)
@inline function Base.getindex(m::WatershedMask, i::Int, j::Int)
    s = @inbounds m.nearest[i, j]
    s != 0 && @inbounds(m.tin[m.k] <= m.tin[s] <= m.tout[m.k])
end

"""
    watershed_contour(dir, loc, nb, dem, threshold, ids)

Outlines of the watersheds of the spill points `loc` (unique, not at the edge of the DEM), following the flow
directions `dir`, as LineStrings with their vertices on cell centres (the crest cells, see `_snap`), starting at
the cell closest to the spill point. Watersheds are nested: the watershed of a depression contains the watersheds
of the depressions spilling into it. `nb` (the entry points) and `threshold` (the spill heights) are not used.
Returns the outlines, the crest cells of the outlines marked with `ids`, the elevations at the vertices,
the index (in `loc`) of the nearest downstream spill point of each cell (0 for none), and the extents.

Nested and cascading watersheds contain each other, so tracing each of them (e.g. by flood filling) visits
the same cells many times. Here we determine the nearest downstream spill point for all cells once, after which
the membership of a watershed is a lookup (see `WatershedMask`). Contour, given a start, only visits the cells
along the outline of the watershed, so the work is proportional to the outline instead of the area.
"""
function watershed_contour(dir, loc::Vector{CartesianIndex{2}}, nb::Vector{CartesianIndex{2}}, dem, threshold, ids)
    # Each spill point should define a single watershed
    allunique(loc) || throw(ArgumentError("spill points in `loc` should be unique"))
    # Spill points at the edge (outlets) have no entry point and no outer outline
    all(dir[p] != CartesianIndex(0, 0) for p in loc) || throw(ArgumentError("spill points in `loc` should not be outlets"))
    # Index the underlying array, as indexing a GeoArray is much slower
    A = parent(dem)

    nearest = _nearest_spill(dir, loc)
    tin, tout, parents = _spilltree(dir, loc, nearest)
    firstcells = _firstcells(nearest, tin, parents)

    # Levee areas
    levees = similar(A, Int32)
    fill!(levees, 0)

    output = Vector{Union{Missing,GeoInterface.Wrappers.LineString}}()
    x, y = axes(A)
    extents = Vector{Tuple{CartesianIndex{2},CartesianIndex{2}}}()
    elevations = Vector{Vector{eltype(A)}}()

    @showprogress desc = "Drawing watershed boundaries" for (k, (loci, id)) in enumerate(zip(loc, ids))
        mask = WatershedMask(nearest, tin, tout, k)
        # Start the contour at the first cell of the watershed, see `_outer_start`
        xi, yi = _outer_start(firstcells[k])

        # Level 0.5 puts the line on the edges between cells, which are shared
        # with the neighboring watersheds. The vertices are moved to cell centres below.
        c = Contour.contour(x, y, mask, 0.5, start=(xi, yi), lazy=true)
        ll = Contour.lines(c)
        if isempty(ll)
            @warn "Empty line"
            push!(output, missing)
            push!(elevations, eltype(A)[])
            push!(extents, (loci, loci))
        else
            @assert length(ll) == 1
            vertices = Contour.vertices(first(ll))
            # Crest cells on both sides of the edge, so the vertices are on cell centres
            cells = _rotate(_snap(vertices, A), loci)
            for cell in cells
                levees[cell] = id
            end
            push!(output,
                GeoInterface.Wrappers.LineString(dem.f.([Tuple(I) .- 0.5 for I in cells]))
            )
            # Vertices are on cell centres, so no rounding (which would round 12.5 to 12, but 13.5 to 14)
            push!(elevations, A[cells])
            # The extent of the watershed, from the cells just inside the outline.
            # Holes lie inside the outline, so don't change the extent.
            extent = (loci, loci)
            for (i, j) in vertices
                a = CartesianIndex(floor(Int, i), floor(Int, j))
                b = CartesianIndex(ceil(Int, i), ceil(Int, j))
                inside = mask[a] ? a : b
                extent = (min(extent[1], inside), max(extent[2], inside))
            end
            push!(extents, extent)
        end
    end
    output, levees, elevations, nearest, extents
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

"""Cells of the vertices of a line on cell centres (as the outlines of `watershed_contour`) of `dem`."""
function _cells(line, dem)
    finv = inv(dem.f)
    # Cell centres are at index - 0.5 in the coordinates of the geotransform
    [CartesianIndex(round.(Int, finv(SVector(GeoInterface.x(p), GeoInterface.y(p))) .+ 0.5)...) for p in GeoInterface.getpoint(line)]
end

"""
    multiscale_tpi(dem, radii=(30, 90, 150))

Topographic position index (TPI): the elevation minus the mean elevation of the other cells in a window,
as the maximum over windows with radii `radii` (m), each rounded to whole cells (at least 1).
A cell is thus higher than its surroundings (positive) if it is higher at any of these scales.
Computed with a summed area table, so large windows are as fast as small ones. At the edges of the DEM
the windows are clipped.
"""
function multiscale_tpi(dem, radii=(30, 90, 150))
    A = parent(dem)
    cs = mean(abs.(Geomorphometry.cellsize(dem)))
    cells = unique(max.(1, round.(Int, collect(radii) ./ cs)))
    # Summed area table, with a zero first row and column
    S = zeros(size(A) .+ 1)
    S[2:end, 2:end] .= cumsum(cumsum(Float64.(A); dims=1); dims=2)
    nx, ny = size(A)
    out = fill(-Inf32, size(A))
    for r in cells
        Threads.@threads for j in 1:ny
            y1, y2 = max(1, j - r), min(ny, j + r)
            for i in 1:nx
                x1, x2 = max(1, i - r), min(nx, i + r)
                s = S[x2+1, y2+1] - S[x1, y2+1] - S[x2+1, y1] + S[x1, y1]
                n = (x2 - x1 + 1) * (y2 - y1 + 1)
                out[i, j] = max(out[i, j], A[i, j] - (s - A[i, j]) / (n - 1))
            end
        end
    end
    out
end

"""
The part of a closed line connected to its start (and end), up to the first and from the last vertex not `keep`.
The line itself if all vertices are kept, `missing` if less than two are.
"""
function _connected(line, keep)
    first_ = findfirst(!, keep)
    isnothing(first_) && return line
    last_ = findlast(!, keep)
    coords = collect(GeoInterface.getpoint(line))
    # Of a closed line, the last vertex is the first one, so it's not repeated
    closed = GeoInterface.x(coords[begin]) == GeoInterface.x(coords[end]) && GeoInterface.y(coords[begin]) == GeoInterface.y(coords[end])
    newcoords = vcat(coords[last_+1:end-closed], coords[begin:first_-1])
    length(newcoords) < 2 && return missing
    GeoInterface.Wrappers.LineString(newcoords)
end

"""
    levee_line(line, tpis, elevations, spillheight; max_height=nothing)

The levee part of a watershed outline `line` (starting at the spill point, as `watershed_contour`): the part
connected to the spill point of which the TPI (`tpis`, at the vertices, see `multiscale_tpi`) is positive,
and with `max_height` (m), also at most `max_height` above `spillheight` (`elevations` at the vertices).
`missing` if there is no such part.
"""
function levee_line(line, tpis, elevations, spillheight; max_height=nothing)
    ismissing(line) && return missing
    keep = tpis .>= 0
    isnothing(max_height) || (keep .&= elevations .<= spillheight + max_height)
    _connected(line, keep)
end

"""
    limit_line(line, elevations, spillheight, threshold=3)

The part of a watershed outline `line` (starting at the spill point) connected to the spill point,
up to where the `elevations` (at the vertices) are more than `threshold` (m) above `spillheight`.
`missing` if there is no such part.
"""
function limit_line(line, elevations, spillheight, threshold=3)
    ismissing(line) && return missing
    _connected(line, elevations .<= spillheight + threshold)
end


"""
    savegpkg(fname, data, dem)

Write the output of `breach` to the geopackage `fname`.gpkg, with a layer for the spill points, the depressions,
the watersheds and the levees. Rows without a geometry (e.g. no levee) are left out of that layer.
Without depressions (an empty `data`), nothing is written.
"""
function savegpkg(fname, data, dem)
    isempty(data) && (@warn "No depressions to write to $fname.gpkg"; return nothing)
    sort!(data, :level)

    data[!, :point] = coords.(Ref(dem), Tuple.(data.spillpoint), Ref(GeoArrays.Center()))
    fdata = select(data, Not(:levee, :spillpoint, :depressionpoint, :entrypoint, :ex, :geometry, :startex, :watershed))
    GeoDataFrames.write("$fname.gpkg", fdata; layer_name="spillpoint", geom_columns=(:point,), crs=EPSG(4326), options=Dict("FID" => "label", "OVERWRITE" => "YES"))

    fdata = _withgeometry(select(data, Not(:levee, :spillpoint, :depressionpoint, :entrypoint, :ex, :point, :startex, :watershed)), :geometry)
    GeoDataFrames.write("$fname.gpkg", fdata; update=true, layer_name="depression", geom_columns=(:geometry,), crs=EPSG(4326), options=Dict("FID" => "label"))

    fdata = _withgeometry(select(data, Not(:levee, :spillpoint, :depressionpoint, :entrypoint, :ex, :point, :startex, :geometry)), :watershed)
    GeoDataFrames.write("$fname.gpkg", fdata; update=true, layer_name="watershed", geom_columns=(:watershed,), crs=EPSG(4326), options=Dict("FID" => "label"))

    fdata = _withgeometry(select(data, Not(:spillpoint, :depressionpoint, :entrypoint, :ex, :point, :startex, :geometry, :watershed)), :levee)
    GeoDataFrames.write("$fname.gpkg", fdata; update=true, layer_name="levee", geom_columns=(:levee,), crs=EPSG(4326), options=Dict("FID" => "label"))

end

"""The rows of `df` with a geometry in column `col` (GeoDataFrames can't write missing geometries)."""
function _withgeometry(df, col)
    df = df[.!ismissing.(df[!, col]), :]
    disallowmissing!(df, col)
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

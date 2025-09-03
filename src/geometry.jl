document[Symbol("Geometry")] = Symbol[]

import StaticArrays
import MillerExtendedHarmonic: MXH
import LinearAlgebra

"""
    centroid(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}

Calculate centroid of polygon
"""
function centroid(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}
    add_endpoint = !((x[1] ≈ x[end]) && (y[1] ≈ y[end]))
    A = 0.5 * sum(x[i] * y[i+1] - x[i+1] * y[i] for i in 1:length(x)-1)
    add_endpoint && (A += 0.5 * (x[end] * y[1] - x[1] * y[end]))

    x_c = -sum(0.25 * (x[i+1] - x[i]) * (y[i+1] + y[i]) * (x[i+1] + x[i]) for i in 1:length(x)-1) / A
    add_endpoint && (x_c -= 0.25 * (x[1] - x[end]) * (y[1] + y[end]) * (x[1] + x[end]) / A)

    y_c = sum(0.25 * (y[i+1] - y[i]) * (x[i+1] + x[i]) * (y[i+1] + y[i]) for i in 1:length(x)-1) / A
    add_endpoint && (y_c += 0.25 * (y[1] - y[end]) * (x[1] + x[end]) * (y[1] + y[end]) / A)

    return x_c, y_c
end

@compat public centroid
push!(document[Symbol("Geometry")], :centroid)

"""
    perimeter(r::AbstractVector{T}, z::AbstractVector{T})::T where {T<:Real}

Calculate the perimeter of a polygon
"""
function perimeter(r::AbstractVector{T}, z::AbstractVector{T})::T where {T<:Real}
    @assert length(r) == length(z) error("Vectors must be of the same length")
    n = length(r)
    perimeter = 0.0
    for i in 1:n-1
        dx = r[i+1] - r[i]
        dy = z[i+1] - z[i]
        perimeter += sqrt(dx^2 + dy^2)
    end

    # If open, add distance from last point to first point
    if is_open_polygon(r, z)
        dx = r[1] - r[end]
        dy = z[1] - z[end]
        perimeter += sqrt(dx^2 + dy^2)
    end

    return perimeter
end

@compat public perimeter
push!(document[Symbol("Geometry")], :perimeter)

"""
    area(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}

Calculate area of polygon
"""
function area(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}
    @views x1 = x[1:end-1]
    @views x2 = x[2:end]
    @views y1 = y[1:end-1]
    @views y2 = y[2:end]
    return abs.(sum(x1 .* y2) - sum(y1 .* x2)) ./ 2
end

@compat public area
push!(document[Symbol("Geometry")], :area)

"""
    revolution_volume(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}

Calculate volume of polygon revolved around x=0
"""
function revolution_volume(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}
    return area(x, y) * 2pi * centroid(x, y)[1]
end

@compat public revolution_volume
push!(document[Symbol("Geometry")], :revolution_volume)

"""
    intersection_angles(
        path1_r::AbstractVector{T},
        path1_z::AbstractVector{T},
        path2_r::AbstractVector{T},
        path2_z::AbstractVector{T},
        intersection_indexes::Vector{StaticArrays.SVector{2,Int}};
        mod_pi::Bool=true
    ) where {T<:Real}

returns angles of intersections between two paths and intersection_indexes given by intersection() function
"""
function intersection_angles(
    path1_r::AbstractVector{T},
    path1_z::AbstractVector{T},
    path2_r::AbstractVector{T},
    path2_z::AbstractVector{T},
    intersection_indexes::Vector{StaticArrays.SVector{2,Int}};
    mod_pi::Bool=true
) where {T<:Real}
    n = length(intersection_indexes)
    angles = Vector{T}(undef, n)

    for (i, index) in enumerate(intersection_indexes)
        r1, z1 = path1_r[index[1]], path1_z[index[1]]
        r1_next, z1_next = path1_r[index[1]+1], path1_z[index[1]+1]
        r2, z2 = path2_r[index[2]], path2_z[index[2]]
        r2_next, z2_next = path2_r[index[2]+1], path2_z[index[2]+1]

        angle = mod(angle_between_two_vectors((r1, z1), (r1_next, z1_next), (r2, z2), (r2_next, z2_next)), π)
        if angle > (π / 2.0) && mod_pi
            angle = π - angle
        end
        angles[i] = angle
    end

    return angles
end

@compat public intersection_angles
push!(document[Symbol("Geometry")], :intersection_angles)

"""
    intersection(
        l1_x::AbstractVector{T},
        l1_y::AbstractVector{T},
        l2_x::AbstractVector{T},
        l2_y::AbstractVector{T}
    ) where {T<:Real}

Intersections between two 2D paths, returns list of (x,y) intersection indexes and crossing points
"""
function intersection(
    l1_x::AbstractVector{T},
    l1_y::AbstractVector{T},
    l2_x::AbstractVector{T},
    l2_y::AbstractVector{T}
) where {T<:Real}

    indexes = StaticArrays.SVector{2,Int}[]
    crossings = StaticArrays.SVector{2,T}[]

    for k1 in 1:(length(l1_x)-1)
        s1_s = StaticArrays.@SVector [l1_x[k1], l1_y[k1]]
        s1_e = StaticArrays.@SVector [l1_x[k1+1], l1_y[k1+1]]
        for k2 in 1:(length(l2_x)-1)
            s2_s = StaticArrays.@SVector [l2_x[k2], l2_y[k2]]
            s2_e = StaticArrays.@SVector [l2_x[k2+1], l2_y[k2+1]]
            #crossing = _seg_intersect(s1_s, s1_e, s2_s, s2_e)
            if _intersect(s1_s, s1_e, s2_s, s2_e)
                crossing = _seg_intersect(s1_s, s1_e, s2_s, s2_e; does_intersect=true)
                push!(indexes, (k1, k2))
                push!(crossings, (crossing[1], crossing[2]))
            end
        end
    end

    return (indexes=indexes, crossings=crossings)
end

"""
    intersection(
        l1_x::AbstractVector{T},
        l1_y::AbstractVector{T},
        l2_x::AbstractVector{T},
        l2_y::AbstractVector{T},
        tolerance::Float64) where {T<:Real}

Intersections between two 2D paths, returns list of (x,y) intersection indexes and crossing points

Endpoints crossings are checked with some tolerance
"""
function intersection(
    l1_x::AbstractVector{T},
    l1_y::AbstractVector{T},
    l2_x::AbstractVector{T},
    l2_y::AbstractVector{T},
    tolerance::Float64) where {T<:Real}

    indexes, crossings = intersection(l1_x, l1_y, l2_x, l2_y)

    if all(k1 != 1 for (k1, k2) in indexes)
        for k2 in 1:(length(l2_x)-1)
            if point_to_segment_distance(l1_x[1], l1_y[1], l2_x[k2], l2_y[k2], l2_x[k2+1], l2_y[k2+1]) < tolerance
                pushfirst!(indexes, StaticArrays.SVector(1, k2))
                pushfirst!(crossings, StaticArrays.SVector(l1_x[1], l1_y[1]))
                break
            end
        end
    end

    if all(k1 != length(l1_x) - 1 for (k1, k2) in indexes)
        for k2 in 1:(length(l2_x)-1)
            if point_to_segment_distance(l1_x[end], l1_y[end], l2_x[k2], l2_y[k2], l2_x[k2+1], l2_y[k2+1]) < tolerance
                push!(indexes, StaticArrays.SVector(length(l1_x) - 1, k2))
                push!(crossings, StaticArrays.SVector(l1_x[end], l1_y[end]))
                break
            end
        end
    end

    # plot(l1_x,l1_y)
    # plot!(l2_x,l2_y)
    # scatter!([cr[1] for cr in crossings],[cr[2] for cr in crossings])
    # display(plot!())

    return (indexes=indexes, crossings=crossings)
end

@compat public intersection
push!(document[Symbol("Geometry")], :intersection)

function intersects(
    l1_x::AbstractVector{<:Real},
    l1_y::AbstractVector{<:Real},
    l2_x::AbstractVector{<:Real},
    l2_y::AbstractVector{<:Real})

    @assert length(l1_x) == length(l1_y)
    @assert length(l2_x) == length(l2_y)
    for k1 in eachindex(l1_x)[1:end-1]
        @inbounds s1_s = StaticArrays.@SVector [l1_x[k1], l1_y[k1]]
        @inbounds s1_e = StaticArrays.@SVector [l1_x[k1+1], l1_y[k1+1]]
        for k2 in eachindex(l2_x)[1:end-1]
            @inbounds s2_s = StaticArrays.@SVector [l2_x[k2], l2_y[k2]]
            @inbounds s2_e = StaticArrays.@SVector [l2_x[k2+1], l2_y[k2+1]]
            _intersect(s1_s, s1_e, s2_s, s2_e) && return true
        end
    end
    return false
end

@inline function _ccw(A, B, C)
    return (C[2] - A[2]) * (B[1] - A[1]) >= (B[2] - A[2]) * (C[1] - A[1])
end

@inline function _ccw(C_A, B_A)
    return (C_A[2] * B_A[1]) >= (B_A[2] * C_A[1])
end

@inline function _out_of_bounds(A, B, C, D)
    abxl, abxu = A[1] < B[1] ? (A[1], B[1]) : (B[1], A[1])
    cdxl, cdxu = C[1] < D[1] ? (C[1], D[1]) : (D[1], C[1])
    (abxu < cdxl || abxl > cdxu) && return true
    abyl, abyu = A[2] < B[2] ? (A[2], B[2]) : (B[2], A[2])
    cdyl, cdyu = C[2] < D[2] ? (C[2], D[2]) : (D[2], C[2])
    return abyu < cdyl || abyl > cdyu
end

@inline function _intersect(A, B, C, D)
    _out_of_bounds(A, B, C, D) && return false
    return (_ccw(A, C, D) != _ccw(B, C, D)) && (_ccw(A, B, C) != _ccw(A, B, D))
end

@inline function _intersect(A::T, B::T, C::T, D::T) where {T<:StaticArrays.StaticVector{2,<:Real}}
    _out_of_bounds(A, B, C, D) && return false
    B_A = B - A
    C_A = C - A
    D_A = D - A
    C_B = C - B
    D_B = D - B
    return (_ccw(D_A, C_A) != _ccw(D_B, C_B)) && (_ccw(C_A, B_A) != _ccw(D_A, B_A))
end

@inline function _perp(a)
    return StaticArrays.@SVector[-a[2], a[1]]
end

function _seg_intersect(a1::T, a2::T, b1::T, b2::T; does_intersect::Bool=_intersect(a1, a2, b1, b2)) where {T<:AbstractVector{<:Real}}
    if !does_intersect
        return nothing
    end
    da = a2 - a1
    db = b2 - b1
    dp = a1 - b1
    dap = _perp(da)
    denom = LinearAlgebra.dot(dap, db)
    num = LinearAlgebra.dot(dap, dp)
    return (num / denom) * db + b1
end

"""
    intersection_split(
        l1_x::AbstractVector{T},
        l1_y::AbstractVector{T},
        l2_x::AbstractVector{T},
        l2_y::AbstractVector{T}) where {T<:Real}

Returns vector of segments of l1_x,l1_y split at the intersections with l2_x,l2_y
"""
function intersection_split(
    l1_x::AbstractVector{T},
    l1_y::AbstractVector{T},
    l2_x::AbstractVector{T},
    l2_y::AbstractVector{T}) where {T<:Real}

    indexes, crossings = intersection(l1_x, l1_y, l2_x, l2_y)
    segments = Vector{@NamedTuple{r::Vector{T}, z::Vector{T}}}(undef, max(length(indexes), 1))
    if isempty(indexes)
        segments[1] = (r=l1_x, z=l1_y)
    else
        Nind = length(indexes)
        indexes1 = [(k <= Nind ? indexes[k][1] : indexes[1][1] + length(l1_x)) for k in 1:Nind+1]

        for k in 1:length(indexes)
            krange = indexes1[k]+1:indexes1[k+1]
            Nk = length(krange)
            kk = k + 1
            if kk > length(crossings)
                kk = kk - length(crossings)
            end

            r = Vector{T}(undef, Nk + 2)
            z = similar(r)

            r[1], z[1] = crossings[k]
            for (j, ind) in enumerate(krange)
                r[j+1] = getindex_circular(l1_x, ind)
                z[j+1] = getindex_circular(l1_y, ind)
            end
            r[end], z[end] = crossings[kk]
            # if segment is made only of intersections, put a point in the middle of the segment
            if length(r) == 2
                segments[k] = (r=[r[1], (r[1] + r[end]) / 2, r[end]], z=[z[1], (z[1] + z[end]) / 2, z[end]])
            else
                segments[k] = (r=r, z=z)
            end
        end
    end

    return segments
end

@compat public intersection_split
push!(document[Symbol("Geometry")], :intersection_split)

"""
    point_to_line_distance(x0::Real, y0::Real, x1::Real, y1::Real, x2::Real, y2::Real)

Distance of point (x0,y0) from line defined by points (x1,y1) and (x2,y2)
"""
function point_to_line_distance(x0::Real, y0::Real, x1::Real, y1::Real, x2::Real, y2::Real)
    return abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt((y2 - y1)^2 + (x2 - x1)^2)
end

@compat public point_to_line_distance
push!(document[Symbol("Geometry")], :point_to_line_distance)

"""
    closest_point_to_segment(x0::Real, y0::Real, x1::Real, y1::Real, x2::Real, y2::Real)

Closest point on segment defined by points (x1,y1) and (x2,y2) to point (x0,y0)
"""
function closest_point_to_segment(x0::Real, y0::Real, x1::Real, y1::Real, x2::Real, y2::Real)
    # Calculate the squared length of the segment
    segment_length_squared = (x2 - x1)^2 + (y2 - y1)^2

    if segment_length_squared == 0.0
        # The segment is just a point, return (x1,y1) [= (x2,y2)]
        return (closest_x=x1, closest_y=y1)
    end

    # Compute the projection of the point onto the line defined by the segment
    t = ((x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1)) / segment_length_squared

    # Clamp t to the range [0, 1] to stay within the segment
    t = clamp(t, 0.0, 1.0)

    # Find the closest point on the segment to the original point
    closest_x = x1 + t * (x2 - x1)
    closest_y = y1 + t * (y2 - y1)

    return (closest_x=closest_x, closest_y=closest_y)
end

@compat public closest_point_to_segment
push!(document[Symbol("Geometry")], :closest_point_to_segment)

"""
    point_to_segment_distance(x0::Real, y0::Real, x1::Real, y1::Real, x2::Real, y2::Real)

Distance of point (x0,y0) from segment defined by points (x1,y1) and (x2,y2)
"""
function point_to_segment_distance(x0::Real, y0::Real, x1::Real, y1::Real, x2::Real, y2::Real)
    closest_x, closest_y = closest_point_to_segment(x0, y0, x1, y1, x2, y2)

    # Compute the distance from the point to the closest point on the segment
    distance = hypot(x0 - closest_x, y0 - closest_y)

    return distance
end

@compat public point_to_segment_distance
push!(document[Symbol("Geometry")], :point_to_segment_distance)

"""
    point_to_path_distance(x0::Real, y0::Real, x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Distance of point (x0,y0) from path defined by vectors x and y
"""
function point_to_path_distance(x0::Real, y0::Real, x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) == length(y)
    d = Inf
    @inbounds for i in 1:length(x)-1
        x1 = x[i]
        y1 = y[i]
        x2 = x[i+1]
        y2 = y[i+1]
        dd = point_to_segment_distance(x0, y0, x1, y1, x2, y2)
        if dd < d
            d = dd
        end
    end
    return d
end

@compat public point_to_path_distance
push!(document[Symbol("Geometry")], :point_to_path_distance)

"""
    rdp_simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, epsilon::T) where {T<:Real}

Simplifies a 2D line represented by arrays of x and y coordinates using the
Ramer-Douglas-Peucker algorithm. The `epsilon` parameter controls the maximum distance
allowed between a point on the original line and its simplified representation.
"""
function rdp_simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, epsilon::T) where {T<:Real}
    #@assert x[1] != x[end] || y[1] != y[end] "p[1] = ($(x[1]),$(y[1]))  p[end] = ($(x[end]),$(y[end])) "
    closed = false
    if x[1] == x[end] && y[1] == y[end]
        closed = true
        x = x[1:end-1]
        y = y[1:end-1]
    end
    @assert length(x) == length(y) "Input arrays must have at least 3 elements"

    n = length(x)
    if n <= 3
        X, Y = x, y
    else
        # Find the point with the maximum distance from the line between the first and last points
        dmax = 0
        index = 0
        for i in 2:n-1
            d = point_to_segment_distance(x[i], y[i], x[1], y[1], x[end], y[end])
            if d > dmax
                index = i
                dmax = d
            end
        end

        # If the maximum distance is greater than epsilon, recursively simplify
        if dmax > epsilon
            # Recursive call to simplify the line segments
            left_points = rdp_simplify_2d_path(x[1:index], y[1:index], epsilon)
            right_points = rdp_simplify_2d_path(x[index:end], y[index:end], epsilon)

            # Combine the simplified line segments
            x_simplified = [left_points[1]; right_points[1][2:end]]
            y_simplified = [left_points[2]; right_points[2][2:end]]
            X, Y = x_simplified, y_simplified
        else
            # If the maximum distance is less than epsilon, return the original line segment
            X, Y = x[[1, end]], y[[1, end]]
        end
    end

    if closed
        return T[X; X[1]], T[Y; Y[1]]
    else
        return X, Y
    end
end

@compat public rdp_simplify_2d_path
push!(document[Symbol("Geometry")], :rdp_simplify_2d_path)

"""
    rwa_simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, epsilon::T) where {T<:Real}

Simplifies a 2D line represented by arrays of x and y coordinates using the Reumann-Witkam Algorithm algorithm.
This algorithm uses a threshold value to determine which points to keep in the path.
Points are kept if the angle between the previous and next line segments is greater than the threshold,
and removed if it is less than or equal to the threshold.
"""
function rwa_simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, threshold::T) where {T<:Real}
    points = [(x[i], y[i]) for i in eachindex(x)]
    simplified_points = [points[1]]
    prev_angle = 0
    for i in 2:length(points)-1
        angle = calculate_angle(points[i-1], points[i], points[i+1])
        if abs(angle - prev_angle) > threshold
            push!(simplified_points, points[i])
            prev_angle = angle
        end
    end
    push!(simplified_points, points[end])
    simplified_x = [p[1] for p in simplified_points]
    simplified_y = [p[2] for p in simplified_points]
    return simplified_x, simplified_y
end


@compat public rwa_simplify_2d_path
push!(document[Symbol("Geometry")], :rwa_simplify_2d_path)

"""
    calculate_angle(p1::T, p2::T, p3::T) where {T}

Calculate the angle between three points
"""
function calculate_angle(p1::Tuple{T,T}, p2::Tuple{T,T}, p3::Tuple{T,T}) where {T<:Real}
    v1 = [p2[1] - p1[1], p2[2] - p1[2]]
    v2 = [p3[1] - p2[1], p3[2] - p2[2]]
    dot_product = dot(v1, v2)
    magnitude_product = norm(v1) * norm(v2)
    return acosd(min(dot_product / magnitude_product, one(T)))
end

@compat public calculate_angle
push!(document[Symbol("Geometry")], :calculate_angle)

"""
    simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, simplification_factor::T; model::Symbol=:distance)

Simplify 2D path by `:curvature` (Reumann-Witkam Algorithm) or `:distance` (Ramer-Douglas-Peucker) algorithms
"""
function simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, simplification_factor::T; model::Symbol=:distance) where {T<:Real}
    if model == :curvature
        return rwa_simplify_2d_path(x, y, simplification_factor)
    elseif model == :distance
        return rdp_simplify_2d_path(x, y, simplification_factor)
    else
        error("simplify_2d_line model can be either :curvature or :distance")
    end
end

@compat public simplify_2d_path
push!(document[Symbol("Geometry")], :simplify_2d_path)

"""
    resample_2d_path(
        x::AbstractVector{T},
        y::AbstractVector{T};
        step::Float64=0.0,
        n_points::Integer=0,
        curvature_weight::Float64=0.0,
        retain_extrema::Bool=false,
        retain_original_xy::Bool=false,
        method::Symbol=:cubic) where {T<:Real}

Resample 2D line with uniform stepping (or number of points)
with option to add more points where curvature is highest
and option to retain extrema in x and y (in these cases stepping is not constant anymore!)
"""
function resample_2d_path(
    x::AbstractVector{T},
    y::AbstractVector{T};
    step::Float64=0.0,
    n_points::Integer=0,
    curvature_weight::Float64=0.0,
    retain_extrema::Bool=false,
    retain_original_xy::Bool=false,
    method::Symbol=:cubic) where {T<:Real}

    t = similar(x)
    t[1] = zero(T)
    for i in 2:length(t)
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        t[i] = t[i-1] + sqrt(dx^2 + dy^2)
    end

    if curvature_weight != 0.0
        @assert 0.0 < curvature_weight < 1.0
        c = moving_average(abs.(curvature(x, y)), round(Int, length(x) / 2.0 * (1.0 - curvature_weight), RoundUp))
        c = c ./ maximum(c)
        c = cumsum((1.0 - curvature_weight) .+ c * curvature_weight)
        t = (c .- c[1]) ./ (c[end] - c[1]) .* (t[end] - t[1]) .+ t[1]
    end

    if n_points === 0
        if step !== 0.0
            n_points = ceil(Int, t[end] / step)
        else
            n_points = length(x)
        end
    end

    # points of interest
    ti = range(t[1], t[end], n_points)
    if retain_original_xy
        ti = sort!(unique!(vcat(t, ti)))
    end

    # interpolate
    xi = interp1d(t, x, method).(ti)
    yi = interp1d(t, y, method).(ti)

    # retain extrema in x and y
    if retain_extrema
        ti = collect(ti)
        for k in (argmax(x), argmax(y), argmin(x), argmin(y))
            index = argmin_abs(ti, t[k])
            ti[index] = t[k]
            xi[index] = x[k]
            yi[index] = y[k]
        end
    end

    # if original path closed, make sure resampled path closes too, independently of interpolation method used
    if is_closed_polygon(x, y)
        xi[end] = xi[1]
        yi[end] = yi[1]
    end

    return xi, yi
end

@compat public resample_2d_path
push!(document[Symbol("Geometry")], :resample_2d_path)

"""
    resample_plasma_boundary(
        x::AbstractVector{T},
        y::AbstractVector{T};
        step::Float64=0.0,
        n_points::Integer=0,
        curvature_weight::Float64=0.0,
        retain_extrema::Bool=true,
        retain_original_xy::Bool=false,
        method::Symbol=:linear) where {T<:Real}

Like resample_2d_path but with `retain_extrema=true` and `method=:linear` as defaults
"""
function resample_plasma_boundary(
    x::AbstractVector{T},
    y::AbstractVector{T};
    step::Float64=0.0,
    n_points::Integer=0,
    curvature_weight::Float64=0.0,
    retain_extrema::Bool=true,
    retain_original_xy::Bool=false,
    method::Symbol=:linear) where {T<:Real}
    x, y = resample_2d_path(closed_polygon(x, y).rz...; step, n_points, curvature_weight, retain_extrema, retain_original_xy, method)
    return x, y
end

@compat public resample_plasma_boundary
push!(document[Symbol("Geometry")], :resample_plasma_boundary)

"""
    is_z_offset(pr::Vector{T}, pz::Vector{T}; order::Int=4, precision::Float64=1E-3) where {T<:Real}

Returns true the shape is offset from z=0 (ie. does not have mxh.z0=0)
"""
function is_z_offset(pr::Vector{T}, pz::Vector{T}; order::Int=4, precision::Float64=1E-3) where {T<:Real}
    if abs(maximum(pz) + minimum(pz)) / 2 > precision
        return true
    end
    pr = deepcopy(pr)
    pz = deepcopy(pz)
    IMAS.reorder_flux_surface!(pr, pz)
    return is_z_offset(MXH(pr, pz, order); precision)
end

"""
    is_z_offset(mxh::MXH; precision::Float64=1E-3)
"""
function is_z_offset(mxh::MXH; precision::Float64=1E-3)
    return abs(mxh.Z0) > precision
end

@compat public is_z_offset
push!(document[Symbol("Geometry")], :is_z_offset)

"""
    is_updown_symmetric(pr::Vector{T}, pz::Vector{T}; order::Int=4, precision::Float64=1E-3) where {T<:Real}

Returns true if boundary is updown symmetric (independent of mxh.z0)
"""
function is_updown_symmetric(pr::Vector{T}, pz::Vector{T}; order::Int=4, precision::Float64=1E-3) where {T<:Real}
    pr = deepcopy(pr)
    pz = deepcopy(pz)
    IMAS.reorder_flux_surface!(pr, pz)
    return is_updown_symmetric(MXH(pr, pz, order); precision)
end

"""
    is_updown_symmetric(mxh::MXH; precision::Float64=1E-3)
"""
function is_updown_symmetric(mxh::MXH; precision::Float64=1E-3)
    return sum(abs.(mxh.c)) / length(mxh.c) < precision
end

@compat public is_updown_symmetric
push!(document[Symbol("Geometry")], :is_updown_symmetric)

"""
    minimum_distance_polygons_vertices(
        R_obj1::AbstractVector{<:T},
        Z_obj1::AbstractVector{<:T},
        R_obj2::AbstractVector{<:T},
        Z_obj2::AbstractVector{<:T};
        return_index::Bool=false) where {T<:Real}

Returns minimum distance between two polygons vertices and index of points on the two polygons
"""
function minimum_distance_polygons_vertices(
    R_obj1::AbstractVector{<:T},
    Z_obj1::AbstractVector{<:T},
    R_obj2::AbstractVector{<:T},
    Z_obj2::AbstractVector{<:T}) where {T<:Real}

    distance2 = Inf
    ik1 = 0
    ik2 = 0
    for k1 in eachindex(R_obj1)
        for k2 in eachindex(R_obj2)
            @inbounds d = (R_obj1[k1] - R_obj2[k2])^2 + (Z_obj1[k1] - Z_obj2[k2])^2
            if distance2 > d
                ik1 = k1
                ik2 = k2
                distance2 = d
            end
        end
    end
    return (distance=sqrt(distance2), k1=ik1, k2=ik2)
end

@compat public minimum_distance_polygons_vertices
push!(document[Symbol("Geometry")], :minimum_distance_polygons_vertices)

"""
    minimum_distance_polygons(
        R_obj1::AbstractVector{<:T},
        Z_obj1::AbstractVector{<:T},
        R_obj2::AbstractVector{<:T},
        Z_obj2::AbstractVector{<:T}) where {T<:Real}

Returns minimum distance between two polygons

NOTE: this is the actual distance, not the distance between the vertices
"""
function minimum_distance_polygons(
    R_obj1::AbstractVector{T},
    Z_obj1::AbstractVector{T},
    R_obj2::AbstractVector{T},
    Z_obj2::AbstractVector{T}) where {T<:Real}

    distance = Inf
    for k1 in eachindex(R_obj1)
        d = point_to_path_distance(R_obj1[k1], Z_obj1[k1], R_obj2, Z_obj2)
        if distance > d
            distance = d
        end
    end

    return distance
end

@compat public minimum_distance_polygons
push!(document[Symbol("Geometry")], :minimum_distance_polygons)

"""
    min_mean_distance_polygons(
        R_obj1::AbstractVector{<:T},
        Z_obj1::AbstractVector{<:T},
        R_obj2::AbstractVector{<:T},
        Z_obj2::AbstractVector{<:T}) where {T<:Real}

Calculate the minimum and mean distances between two polygons in 2D space.

NOTE: this is the actual distance, not the distance between the vertices
"""
function min_mean_distance_polygons(
    R_obj1::AbstractVector{<:T},
    Z_obj1::AbstractVector{<:T},
    R_obj2::AbstractVector{<:T},
    Z_obj2::AbstractVector{<:T}) where {T<:Real}

    mean_distance = 0.0
    min_distance = Inf
    for k1 in eachindex(R_obj1)
        actual_distance = point_to_path_distance(R_obj1[k1], Z_obj1[k1], R_obj2, Z_obj2)
        # Update global minimum distance
        if min_distance > actual_distance
            min_distance = actual_distance
        end
        # Accumulate the difference from the target distance
        mean_distance += actual_distance
    end

    # Get the mean distance error
    mean_distance = mean_distance / length(R_obj1)

    return (min_distance=min_distance, mean_distance=mean_distance)
end

@compat public min_mean_distance_polygons
push!(document[Symbol("Geometry")], :min_mean_distance_polygons)

"""
    curvature(pr::AbstractVector{T}, pz::AbstractVector{T}) where {T<:Real}

Calculate the curvature of a 2D path defined by `pr` and `pz` using a finite difference approximation.

The path is assumed to be closed if the first and last points are the same, and open otherwise.

# Arguments

  - `pr`: Real abstract vector representing the r-coordinates of the path.
  - `pz`: Real abstract vector representing the z-coordinates of the path.

# Returns

  - A vector of the same length as `pr` and `pz` with the calculated curvature values.
"""
function curvature(pr::AbstractVector{T}, pz::AbstractVector{T}) where {T<:Real}
    n = length(pr)
    curvature_res = Vector{T}(undef, n)
    if is_closed_polygon(pr, pz)
        dr1 = pr[end-1] - pr[1]
        dz1 = pz[end-1] - pz[1]
    else
        dr1 = 0.0
        dz1 = 0.0
    end
    a1 = sqrt(dr1^2 + dz1^2) + 1E-32
    for i in 1:n
        dr2 = if i < n
            pr[i+1] - pr[i]
        else
            pr[1] - pr[end]
        end
        dz2 = if i < n
            pz[i+1] - pz[i]
        else
            pz[1] - pz[end]
        end
        a2 = sqrt(dr2^2 + dz2^2) + 1E-32
        curvature_res[i] = (dr1 / a1) * (dz2 / a2) - (dr2 / a2) * (dz1 / a1)
        dr1, dz1, a1 = dr2, dz2, a2
    end
    return curvature_res
end

@compat public curvature
push!(document[Symbol("Geometry")], :curvature)

"""
    angle_between_two_vectors(
        v1_p1::Tuple{T,T},
        v1_p2::Tuple{T,T},
        v2_p1::Tuple{T,T},
        v2_p2::Tuple{T,T}) where {T<:Real}

Returns angle in radiants between two vectors defined by their start and end points
"""
function angle_between_two_vectors(
    v1_p1::Tuple{T,T},
    v1_p2::Tuple{T,T},
    v2_p1::Tuple{T,T},
    v2_p2::Tuple{T,T}) where {T<:Real}

    v1_x = v1_p2[1] - v1_p1[1]
    v1_y = v1_p2[2] - v1_p1[2]
    v2_x = v2_p2[1] - v2_p1[1]
    v2_y = v2_p2[2] - v2_p1[2]

    arg = (v1_x * v2_x + v1_y * v2_y) / (sqrt(v1_x^2 + v1_y^2) * sqrt(v2_x^2 + v2_y^2))
    # limit arg to [-1.0, 1.0], which is guaranteed mathematically by (a · b) / (|a| |b|)
    return acos(max(min(arg, 1.0), -1.0))
end

@compat public angle_between_two_vectors
push!(document[Symbol("Geometry")], :angle_between_two_vectors)

"""
    bisector(v1, v2, v3)

Returns signed unitary bisector given three vertices
"""
function bisector(v1, v2, v3)
    # Convert points to vectors
    a = [v1[1] - v2[1], v1[2] - v2[2]]
    b = [v3[1] - v2[1], v3[2] - v2[2]]

    # Normalize vectors
    @assert norm(a) > 0.0 "$a"
    a /= norm(a)
    @assert norm(b) > 0.0 "$b"
    b /= norm(b)

    # Check if vectors are opposite to each other, indicating a straight line
    if isapprox(a, -b; atol=1e-8)
        # Choose a vector that is perpendicular to a (or b)
        bisect = [-a[2], a[1]]
    else
        # Otherwise, calculate the bisector
        bisect = a + b
        bisect /= norm(bisect)
    end

    # cross product to make bisector always point in/out
    # depending if v1,v2,v3 are clockwise or counter-clockwise
    cross_prod = a[1] * bisect[2] - a[2] * bisect[1]

    return bisect * sign(cross_prod)
end

@compat public bisector
push!(document[Symbol("Geometry")], :bisector)

"""
    polygon_rays(vertices::AbstractVector, extent_a::Float64, extent_b::Float64)

Returns bisecting "rays" (lines) that radiate from vertices of a polygon
"""
function polygon_rays(vertices::AbstractVector, extent_a::Float64, extent_b::Float64)
    @assert is_open_polygon(vertices)

    # Create rays
    rays = []
    for i in eachindex(vertices)
        # Get the current vertex and the adjacent ones
        v1 = vertices[mod(i - 1, length(vertices))+1]
        v2 = vertices[mod(i + 0, length(vertices))+1]
        v3 = vertices[mod(i + 1, length(vertices))+1]

        # Calculate bisector
        bisect = bisector(v1, v2, v3)

        # Define the ray
        ray = [(v2[1] + extent_a * bisect[1], v2[2] + extent_a * bisect[2]), (v2[1] + extent_b * bisect[1], v2[2] + extent_b * bisect[2])]
        push!(rays, ray)
    end

    return rays
end

@compat public polygon_rays
push!(document[Symbol("Geometry")], :polygon_rays)

"""
    split_long_segments(R::AbstractVector{T}, Z::AbstractVector{T}, max_length::Float64) where {T<:Real}

Split long segments of a polygon so that each resulting segment is always <= max_length
"""
function split_long_segments(R::AbstractVector{T}, Z::AbstractVector{T}, max_length::Float64) where {T<:Real}
    opoly = open_polygon(R, Z)

    Rout = [opoly.R[1]]
    Zout = [opoly.Z[1]]
    for k in eachindex(opoly.R)
        r1 = IMAS.getindex_circular(opoly.R, k)
        z1 = IMAS.getindex_circular(opoly.Z, k)
        r2 = IMAS.getindex_circular(opoly.R, k + 1)
        z2 = IMAS.getindex_circular(opoly.Z, k + 1)
        d = sqrt((r2 - r1)^2 + (z2 - z1)^2)
        if d > max_length
            # linear interpolation
            n = round(Int, d / max_length, RoundUp) + 1
            append!(Rout, collect(range(r1 + (r2 - r1) / (n - 1), r2, n - 1)))
            append!(Zout, collect(range(z1 + (z2 - z1) / (n - 1), z2, n - 1)))
        else
            push!(Rout, r2)
            push!(Zout, z2)
        end
    end
    if opoly.was_closed
        return Rout, Zout
    else
        return Rout[1:end-1], Zout[1:end-1]
    end
end

"""
    split_long_segments(R::AbstractVector{T}, Z::AbstractVector{T}, n_points::Int) where {T<:Real}

Split long segments of a polygon so that there are at least n_points in it
"""
function split_long_segments(R::AbstractVector{T}, Z::AbstractVector{T}, n_points::Int) where {T<:Real}
    L = perimeter(R, Z)
    max_length = L / n_points
    return split_long_segments(R, Z, max_length)
end

@compat public split_long_segments
push!(document[Symbol("Geometry")], :split_long_segments)

"""
    thick_line_polygon(r1, z1, r2, z2, thickness1, thickness2)

Generates a closed polygon from a thick line. Returns points of the quadrilateral polygon
"""
function thick_line_polygon(r1::Float64, z1::Float64, r2::Float64, z2::Float64, thickness1::Float64, thickness2::Float64)
    direction = normalize([z2 - z1, -(r2 - r1)]) # Perpendicular direction
    offset1 = direction .* thickness1 / 2
    offset2 = direction .* thickness2 / 2
    p1 = [r1, z1] + offset1
    p2 = [r2, z2] + offset2
    p3 = [r2, z2] - offset2
    p4 = [r1, z1] - offset1
    return [p1, p2, p3, p4, p1]
end

"""
    thick_line_polygon(pr::AbstractVector{Float64}, pz::AbstractVector{Float64}, thickness::AbstractVector{Float64})
"""
function thick_line_polygon(pr::AbstractVector{Float64}, pz::AbstractVector{Float64}, thickness::AbstractVector{Float64})
    PTS = [thick_line_polygon(pr[i], pz[i], pr[i+1], pz[i+1], thickness[i], thickness[i+1]) for i in 1:length(pr)-1]
    x1 = map(x -> x[1][1], PTS)
    y1 = map(x -> x[1][2], PTS)
    x2 = map(x -> x[2][1], PTS)
    y2 = map(x -> x[2][2], PTS)
    x3 = map(x -> x[3][1], PTS)
    y3 = map(x -> x[3][2], PTS)
    x4 = map(x -> x[4][1], PTS)
    y4 = map(x -> x[4][2], PTS)

    r1 = [x1; x2[end]; x3[end]; reverse!(x4); x1[1]]
    z1 = [y1; y2[end]; y3[end]; reverse!(y4); y1[1]]
    r2 = [x1[1]; x2; reverse!(x3); x4[end]; x1[1]]
    z2 = [y1[1]; y2; reverse!(y3); y4[end]; y1[1]]

    return (r=(r1 .+ r2) / 2.0, z=(z1 .+ z2) / 2.0)
end

@compat public thick_line_polygon
push!(document[Symbol("Geometry")], :thick_line_polygon)

@inline function points_isless(p::AbstractVector{T}, q::AbstractVector{T}) where {T}
    return p[1] < q[1] || (p[1] == q[1] && p[2] < q[2])
end

@inline function points_isless(p::Tuple{T,T}, q::Tuple{T,T}) where {T}
    return p[1] < q[1] || (p[1] == q[1] && p[2] < q[2])
end

@inline function isrightturn(p::AbstractVector{T}, q::AbstractVector{T}, r::AbstractVector{T}) where {T}
    return (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1]) < 0.0
end

@inline function isrightturn(p::Tuple{T,T}, q::Tuple{T,T}, r::Tuple{T,T}) where {T}
    return (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1]) < 0.0
end

function halfhull_indices(points::AbstractVector, sorted_indices::AbstractVector{Int})
    hull_indices = Vector{Int}(undef, length(sorted_indices))
    n = 0
    for idx in sorted_indices
        while n > 1 && !isrightturn(points[hull_indices[n-1]], points[hull_indices[n]], points[idx])
            n -= 1
        end
        n += 1
        hull_indices[n] = idx
    end
    return view(hull_indices, 1:n)
end

function grahamscan_indices(points::AbstractVector)
    # Create array of indices and sort them based on point comparison
    indices = collect(1:length(points))
    sort!(indices; lt=(i, j) -> points_isless(points[i], points[j]))

    # Compute upper hull
    upperhull_indices = halfhull_indices(points, indices)

    # Reverse indices for lower hull
    reverse!(indices)
    lowerhull_indices = halfhull_indices(points, indices)

    # Combine hulls (excluding duplicate endpoints)
    return [upperhull_indices; lowerhull_indices[2:end-1]]
end

"""
    convex_hull_indices(xy_points::AbstractVector)

Compute the indices of points that form the convex hull of a set of 2D points.
Returns the indices in counter-clockwise order.

This is the internal function that works with indices only.
"""
function convex_hull_indices(xy_points::AbstractVector)
    if length(xy_points) <= 1
        return collect(1:length(xy_points))
    end
    return grahamscan_indices(xy_points)
end

# Original halfhull function (kept for compatibility)
function halfhull(points::AbstractVector)
    halfhull = similar(points)
    n = 0
    for p in points
        while n > 1 && !isrightturn(halfhull[n-1], halfhull[n], p)
            n -= 1
        end
        n += 1
        halfhull[n] = p
    end
    return view(halfhull, 1:n)
end

# Modified grahamscan! to use index-based approach internally
function grahamscan!(points::AbstractVector)
    hull_indices = grahamscan_indices(points)
    return points[hull_indices]
end

"""
    convex_hull!(xy_points::AbstractVector; closed_polygon::Bool)

Compute the convex hull of a set of 2D points, sorted in counter-clockwise order.
The resulting convex hull forms a closed polygon by appending the first point at the end.

NOTE: The input vector is not modified anymore (uses index-based approach internally)
"""
function convex_hull!(xy_points::AbstractVector; closed_polygon::Bool=false)
    hull_indices = convex_hull_indices(xy_points)
    hull_points = view(xy_points, hull_indices)

    if closed_polygon && !isempty(hull_points)
        return [hull_points; hull_points[1:1]]  # Use view for first element too
    else
        return collect(hull_points)  # Convert view to array for return
    end
end

"""
    convex_hull(xy_points::AbstractVector; closed_polygon::Bool)

Compute the convex hull of a set of 2D points, sorted in counter-clockwise order.
The resulting convex hull forms a closed polygon by appending the first point at the end.
"""
function convex_hull(xy_points::AbstractVector; closed_polygon::Bool=false)
    hull_indices = convex_hull_indices(xy_points)
    hull_points = xy_points[hull_indices]

    if closed_polygon && !isempty(hull_points)
        return push!(hull_points, hull_points[1])
    else
        return hull_points
    end
end

"""
    convex_hull(x::AbstractVector{T}, y::AbstractVector{T}; closed_polygon::Bool) where {T}
"""
function convex_hull(x::AbstractVector{T}, y::AbstractVector{T}; closed_polygon::Bool=false) where {T}
    xy_points = [(xx, yy) for (xx, yy) in zip(x, y)]
    return convex_hull(xy_points; closed_polygon)
end

# Export both the regular function and the index-based function
@compat public convex_hull, convex_hull_indices
push!(document[Symbol("Geometry")], :convex_hull)
push!(document[Symbol("Geometry")], :convex_hull_indices)

"""
    inscribed_polygon(pr::AbstractVector{Float64}, pz::AbstractVector{Float64}; target_area_fraction::Float64=0.9, max_iterations::Int=length(pr))

Build an inscribed polygon by starting with all boundary points and iteratively removing
the point that causes the least area loss, while ensuring no original boundary points
end up inside the inscribed polygon.
"""
function inscribed_polygon(pr::AbstractVector{Float64}, pz::AbstractVector{Float64}; target_area_fraction::Float64=0.9, max_iterations::Int=length(pr))
    # Ensure we have a proper polygon
    @assert length(pr) == length(pz) "pr and pz must have the same length"
    @assert length(pr) >= 3 "Polygon must have at least 3 vertices"
    @assert 0.0 < target_area_fraction <= 1.0 "target_area_fraction must be between 0 and 1"

    # Work with open polygon to avoid duplication
    open_poly = open_polygon(pr, pz)
    original_r = open_poly.R
    original_z = open_poly.Z
    n_original = length(original_r)

    # Pre-allocate buffers
    max_size = n_original + 1  # +1 for closing point
    candidate_r_buffer = Vector{Float64}(undef, max_size)
    candidate_z_buffer = Vector{Float64}(undef, max_size)
    candidate_indices_buffer = Vector{Int}(undef, n_original)

    # Pre-allocate polygon tuple buffer for inpolygon checks
    polygon_tuple_buffer = Vector{Tuple{Float64,Float64}}(undef, max_size)

    # Pre-create original points as tuples for inpolygon checks
    original_points = [(r, z) for (r, z) in zip(original_r, original_z)]

    # Calculate original area and target
    original_r_closed = [original_r; original_r[1]]
    original_z_closed = [original_z; original_z[1]]
    original_area = area(original_r_closed, original_z_closed)
    target_area = target_area_fraction * original_area

    # Start with all points in the inscribed polygon
    inscribed_r = copy(original_r)
    inscribed_z = copy(original_z)
    current_indices = collect(1:n_original)

    # Pre-allocate for area calculation
    area_calc_r = Vector{Float64}(undef, max_size)
    area_calc_z = Vector{Float64}(undef, max_size)

    iteration = 0
    total_removed = 0

    while length(inscribed_r) > 3 && iteration < max_iterations
        iteration += 1

        # Calculate current area using pre-allocated buffer
        n_current = length(inscribed_r)
        @inbounds for i in 1:n_current
            area_calc_r[i] = inscribed_r[i]
            area_calc_z[i] = inscribed_z[i]
        end
        area_calc_r[n_current+1] = inscribed_r[1]  # Close polygon
        area_calc_z[n_current+1] = inscribed_z[1]
        current_area = area(@view(area_calc_r[1:n_current+1]), @view(area_calc_z[1:n_current+1]))

        # Check if we've reached target area
        if current_area <= target_area
            break
        end

        # Try removing each point and find the best candidate
        best_removal_idx = -1
        best_area_loss = Inf

        @inbounds for remove_idx in 1:length(inscribed_r)
            candidate_length = length(inscribed_r) - 1

            # Skip if polygon becomes too small
            if candidate_length < 3
                continue
            end

            # Create candidate polygon in pre-allocated buffer
            # Copy elements before removal point
            for i in 1:(remove_idx-1)
                candidate_r_buffer[i] = inscribed_r[i]
                candidate_z_buffer[i] = inscribed_z[i]
                candidate_indices_buffer[i] = current_indices[i]
            end

            # Copy elements after removal point
            for i in remove_idx:(candidate_length)
                src_idx = i + 1  # Skip the removed element
                candidate_r_buffer[i] = inscribed_r[src_idx]
                candidate_z_buffer[i] = inscribed_z[src_idx]
                candidate_indices_buffer[i] = current_indices[src_idx]
            end

            # Create polygon tuple buffer for inpolygon check (closed polygon)
            for i in 1:candidate_length
                polygon_tuple_buffer[i] = (candidate_r_buffer[i], candidate_z_buffer[i])
            end
            polygon_tuple_buffer[candidate_length+1] = (candidate_r_buffer[1], candidate_z_buffer[1])  # Close

            # Check if candidate polygon contains any original boundary points
            contains_original = false
            for original_point in original_points
                if PolygonOps.inpolygon(original_point, @view(polygon_tuple_buffer[1:candidate_length+1])) == 1
                    contains_original = true
                    break
                end
            end

            if contains_original
                continue  # Invalid removal - skip
            end

            # Calculate area loss using pre-allocated area calculation buffer
            for i in 1:candidate_length
                area_calc_r[i] = candidate_r_buffer[i]
                area_calc_z[i] = candidate_z_buffer[i]
            end
            area_calc_r[candidate_length+1] = candidate_r_buffer[1]  # Close polygon
            area_calc_z[candidate_length+1] = candidate_z_buffer[1]

            candidate_area = area(@view(area_calc_r[1:candidate_length+1]), @view(area_calc_z[1:candidate_length+1]))
            area_loss = current_area - candidate_area

            # Update best candidate if this is better (smaller area loss)
            if area_loss < best_area_loss
                best_area_loss = area_loss
                best_removal_idx = remove_idx
            end
        end

        # Check if we found a valid removal
        if best_removal_idx == -1
            break
        end

        # Apply the best removal by updating inscribed arrays in-place
        remove_idx = best_removal_idx
        new_length = length(inscribed_r) - 1

        # Shift elements to remove the selected point
        @inbounds for i in remove_idx:new_length
            inscribed_r[i] = inscribed_r[i+1]
            inscribed_z[i] = inscribed_z[i+1]
            current_indices[i] = current_indices[i+1]
        end

        # Resize arrays to new length
        resize!(inscribed_r, new_length)
        resize!(inscribed_z, new_length)
        resize!(current_indices, new_length)

        total_removed += 1

        # Calculate new area for early exit check
        @inbounds for i in 1:new_length
            area_calc_r[i] = inscribed_r[i]
            area_calc_z[i] = inscribed_z[i]
        end
        area_calc_r[new_length+1] = inscribed_r[1]
        area_calc_z[new_length+1] = inscribed_z[1]
        new_area = area(@view(area_calc_r[1:new_length+1]), @view(area_calc_z[1:new_length+1]))

        # Early exit check
        if new_area <= target_area
            break
        end
    end

    # Create final closed polygon
    final_length = length(inscribed_r)
    inscribed_r_closed = Vector{Float64}(undef, final_length + 1)
    inscribed_z_closed = Vector{Float64}(undef, final_length + 1)

    @inbounds for i in 1:final_length
        inscribed_r_closed[i] = inscribed_r[i]
        inscribed_z_closed[i] = inscribed_z[i]
    end
    inscribed_r_closed[final_length+1] = inscribed_r[1]
    inscribed_z_closed[final_length+1] = inscribed_z[1]

    final_area = area(inscribed_r_closed, inscribed_z_closed)
    final_fraction = final_area / original_area

    # println("\nOptimized shrinking algorithm completed:")
    # println("  Iterations: $(iteration)")
    # println("  Points removed: $(total_removed)")
    # println("  Final vertices: $(length(inscribed_r))")
    # println("  Final area: $(final_area)")
    # println("  Area fraction: $(100*final_fraction)%")
    # println("  Remaining indices: $(current_indices)")

    return inscribed_r_closed, inscribed_z_closed
end

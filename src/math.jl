import LinearAlgebra
import StaticArrays
import DataInterpolations

"""
    norm01(x::T)::T where {T<:AbstractVector{<:Real}}

Normalize a vector so that the first item in the array is 0 and the last one is 1
This is handy where psi_norm should be used (and IMAS does not define a psi_norm array)
"""
function norm01(x::T)::T where {T<:AbstractVector{<:Real}}
    return (x .- x[1]) ./ (x[end] .- x[1])
end

"""
    to_range(vector::AbstractVector)

Turn a vector into a range (if possible)
"""
function to_range(vector::AbstractVector{<:Real})
    N = length(vector)
    dv = vector[2] - vector[1]
    if !(1 - sum(abs(vector[k] - vector[k-1] - dv) for k in 2:N) / (N - 1) ≈ 1.0)
        error("to_range requires vector data to be equally spaced: $(vector)")
    end
    return range(vector[1], vector[end], N)
end

"""
    meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}

Return coordinate matrices from coordinate vectors
"""
function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    return last.(Iterators.product(y, x)), first.(Iterators.product(y, x))
end

"""
    centroid(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}

Calculate centroid of polygon
"""
function centroid(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}
    A = sum((y[i+1] - y[i]) * 0.5 * (x[i+1] + x[i]) for i in 1:length(x)-1)
    x_c = -sum((x[i+1] - x[i]) * 0.5 * (y[i+1] + y[i]) * 0.5 * (x[i+1] + x[i]) for i in 1:length(x)-1) / A
    y_c = sum((y[i+1] - y[i]) * 0.5 * (x[i+1] + x[i]) * 0.5 * (y[i+1] + y[i]) for i in 1:length(x)-1) / A
    return x_c, y_c
end

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

"""
    revolution_volume(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}

Calculate volume of polygon revolved around x=0
"""
function revolution_volume(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}
    return area(x, y) * 2pi * centroid(x, y)[1]
end

"""
    intersection_angles(path1_r::AbstractVector{T}, path1_z::AbstractVector{T}, path2_r::AbstractVector{T}, path2_z::AbstractVector{T}, intersection_indexes::Vector{Tuple{Int,Int}}) where {T<:Real}

returns angles of intersections between two paths and intersection_indexes given by intersection() function
"""
function intersection_angles(
    path1_r::AbstractVector{T},
    path1_z::AbstractVector{T},
    path2_r::AbstractVector{T},
    path2_z::AbstractVector{T},
    intersection_indexes::Vector{Tuple{Int,Int}};
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

"""
    intersection(
        l1_x::AbstractVector{T},
        l1_y::AbstractVector{T},
        l2_x::AbstractVector{T},
        l2_y::AbstractVector{T}) where {T<:Real}

Intersections between two 2D paths, returns list of (x,y) intersection indexes and crossing points
"""
function intersection(
    l1_x::AbstractVector{T},
    l1_y::AbstractVector{T},
    l2_x::AbstractVector{T},
    l2_y::AbstractVector{T}) where {T<:Real}

    indexes = NTuple{2,Int}[]
    crossings = NTuple{2,T}[]

    for k1 in 1:(length(l1_x)-1)
        s1_s = StaticArrays.@SVector [l1_x[k1], l1_y[k1]]
        s1_e = StaticArrays.@SVector [l1_x[k1+1], l1_y[k1+1]]
        for k2 in 1:(length(l2_x)-1)
            s2_s = StaticArrays.@SVector [l2_x[k2], l2_y[k2]]
            s2_e = StaticArrays.@SVector [l2_x[k2+1], l2_y[k2+1]]
            crossing = _seg_intersect(s1_s, s1_e, s2_s, s2_e)
            if crossing !== nothing
                push!(indexes, (k1, k2))
                push!(crossings, (crossing[1], crossing[2]))
            end
        end
    end

    return (indexes=indexes, crossings=crossings)
end

function _ccw(A, B, C)
    return (C[2] - A[2]) * (B[1] - A[1]) >= (B[2] - A[2]) * (C[1] - A[1])
end

function _intersect(A, B, C, D)
    return (_ccw(A, C, D) != _ccw(B, C, D)) && (_ccw(A, B, C) != _ccw(A, B, D))
end

function _perp(a)
    return [-a[2], a[1]]
end

function _seg_intersect(a1::T, a2::T, b1::T, b2::T) where {T<:AbstractVector{<:Real}}
    if !_intersect(a1, a2, b1, b2)
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
    point_to_line_distance(x0::T, y0::T, x1::T, y1::T, x2::T, y2::T) where {T<:Real}

Distance of point (x0,y0) from line defined by points (x1,y1) and (x2,y2)
"""
function point_to_line_distance(x0::T, y0::T, x1::T, y1::T, x2::T, y2::T) where {T<:Real}
    return abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt((y2 - y1)^2 + (x2 - x1)^2)
end

"""
    point_to_path_distance(x0::T, y0::T, x::Vector{T}, y::Vector{T}) where {T<:Real}

Distance of point (x0,y0) from path defined by vectors x and y
"""
function point_to_path_distance(x0::T, y0::T, x::Vector{T}, y::Vector{T}) where {T<:Real}
    d = Inf
    for i in 1:length(x)-1
        x1 = x[i]
        y1 = y[i]
        x2 = x[i+1]
        y2 = y[i+1]
        dd = point_to_line_distance(x0, y0, x1, y1, x2, y2)
        if dd < d
            d = dd
        end
    end
    return d
end

"""
    rdp_simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, epsilon::T) where {T<:Real}

Simplifies a 2D line represented by arrays of x and y coordinates using the
Ramer-Douglas-Peucker algorithm. The `epsilon` parameter controls the maximum distance
allowed between a point on the original line and its simplified representation.
"""
function rdp_simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, epsilon::T) where {T<:Real}
    @assert x[1] != x[end] || y[1] != y[end]

    n = length(x)
    if n != length(y)
        error("Input arrays must have at least 3 elements and the same length")
    end

    if n <= 3
        return x, y
    else
        # Find the point with the maximum distance from the line between the first and last points
        dmax = 0
        index = 0
        for i in 2:n-1
            d = point_to_line_distance(x[i], y[i], x[1], y[1], x[end], y[end])
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
            return (x_simplified, y_simplified)
        else
            # If the maximum distance is less than epsilon, return the original line segment
            return x[[1, end]], y[[1, end]]
        end
    end
end

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

"""
    moving_average(data::Vector{<:Real}, window_size::Int)

Calculate the moving average of a data vector using a specified window size.
The window size is always rounded up to the closest odd number to maintain symmetry around each data point.
"""
function moving_average(data::Vector{<:Real}, window_size::Int)
    smoothed_data = copy(data)
    window_size = isodd(window_size) ? window_size : window_size + 1
    pad_size = div(window_size - 1, 2)
    if pad_size < 1
        return smoothed_data
    end
    for i in eachindex(data)
        window_start = max(1, i - pad_size)
        window_end = min(length(data), i + pad_size)
        smoothed_data[i] = sum(data[window_start:window_end]) / (window_end - window_start + 1)
    end
    return smoothed_data
end

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
        c = moving_average(abs.(curvature(x, y)), Int(ceil(length(x) / 2.0 * (1.0 - curvature_weight))))
        c = c ./ maximum(c)
        c = cumsum((1.0 - curvature_weight) .+ c * curvature_weight)
        #c = moving_average(c, Int(ceil(length(x) / 20)))
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
            index = argmin(abs.(ti .- t[k]))
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
    x, y = resample_2d_path(x, y; step, n_points, curvature_weight, retain_extrema, retain_original_xy, method)
    return x, y
end

"""
    minimum_distance_two_shapes(
        R_obj1::AbstractVector{<:T},
        Z_obj1::AbstractVector{<:T},
        R_obj2::AbstractVector{<:T},
        Z_obj2::AbstractVector{<:T};
        return_index::Bool=false) where {T<:Real}

Returns minimum distance between two shapes
"""
function minimum_distance_two_shapes(
    R_obj1::AbstractVector{<:T},
    Z_obj1::AbstractVector{<:T},
    R_obj2::AbstractVector{<:T},
    Z_obj2::AbstractVector{<:T};
    return_index::Bool=false) where {T<:Real}

    distance = Inf
    ik1 = 0
    ik2 = 0
    for k1 in eachindex(R_obj1)
        for k2 in eachindex(R_obj2)
            @inbounds d = (R_obj1[k1] - R_obj2[k2])^2 + (Z_obj1[k1] - Z_obj2[k2])^2
            if distance > d
                ik1 = k1
                ik2 = k2
                distance = d
            end
        end
    end
    if return_index
        return ik1, ik2
    else
        return sqrt(distance)
    end
end

"""
    mean_distance_error_two_shapes(
        R_obj1::AbstractVector{<:T},
        Z_obj1::AbstractVector{<:T},
        R_obj2::AbstractVector{<:T},
        Z_obj2::AbstractVector{<:T},
        target_distance::T,
        above_target::Bool=true)

Returns mean error distance between two shapes and a target distance
"""
function mean_distance_error_two_shapes(
    R_obj1::AbstractVector{<:T},
    Z_obj1::AbstractVector{<:T},
    R_obj2::AbstractVector{<:T},
    Z_obj2::AbstractVector{<:T},
    target_distance::T;
    above_target::Bool=false,
    below_target::Bool=false) where {T<:Real}

    mean_distance_error = 0.0
    n = 0
    for k1 in eachindex(R_obj1)
        distance = Inf
        for k2 in eachindex(R_obj2)
            @inbounds d = (R_obj1[k1] - R_obj2[k2])^2 + (Z_obj1[k1] - Z_obj2[k2])^2
            if d < distance
                distance = d
            end
        end
        if above_target && distance > target_distance
            mean_distance_error += (distance - target_distance)^2
            n += 1
        elseif below_target && distance < target_distance
            mean_distance_error += (distance - target_distance)^2
            n += 1
        else
            mean_distance_error += (distance - target_distance)^2
            n += 1
        end
    end
    return sqrt(mean_distance_error) / n
end

function min_mean_distance_error_two_shapes(
    R_obj1::AbstractVector{<:T},
    Z_obj1::AbstractVector{<:T},
    R_obj2::AbstractVector{<:T},
    Z_obj2::AbstractVector{<:T},
    target_distance::T;
    above_target::Bool=false,
    below_target::Bool=false) where {T<:Real}

    min_distance = Inf
    mean_distance_error = 0.0
    n = 0

    for k1 in eachindex(R_obj1)
        for k2 in eachindex(R_obj2)
            @inbounds d = (R_obj1[k1] - R_obj2[k2])^2 + (Z_obj1[k1] - Z_obj2[k2])^2

            # Calculate minimum distance
            if min_distance > d
                min_distance = d
            end

            # Calculate mean error distance
            if above_target && d > target_distance
                mean_distance_error += (d - target_distance)^2
                n += 1
            elseif below_target && d < target_distance
                mean_distance_error += (d - target_distance)^2
                n += 1
            else
                mean_distance_error += (d - target_distance)^2
                n += 1
            end
        end
    end

    # Return results
    min_distance = sqrt(min_distance)
    mean_distance_error = sqrt(mean_distance_error) / n

    return min_distance, mean_distance_error
end

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

"""
    calc_z(x::AbstractVector{<:Real}, f::AbstractVector{<:Real}, method::Symbol)

Returns the gradient scale lengths of vector f on x

The finite difference `method` of the gradient can be one of [:third_order, :second_order, :central, :backward, :forward]

NOTE: the inverse scale length is NEGATIVE for typical density/temperature profiles
"""
function calc_z(x::AbstractVector{<:Real}, f::AbstractVector{<:Real}, method::Symbol)
    g = gradient(x, f; method)
    return g ./ f
end

"""
    integ_z(rho::AbstractVector{<:Real}, z_profile::AbstractVector{<:Real}, bc::Real)

Backward integration of inverse scale length vector with given edge boundary condition
"""
function integ_z(rho::AbstractVector{<:Real}, z_profile::AbstractVector{<:Real}, bc::Real)
    profile_new = similar(rho)
    profile_new[end] = bc
    for i in length(rho)-1:-1:1
        a = rho[i]
        fa = z_profile[i]
        b = rho[i+1]
        fb = z_profile[i+1]
        trapz_integral = (b - a) * (fa + fb) / 2.0
        profile_new[i] = profile_new[i+1] * exp(trapz_integral)
    end
    return profile_new
end

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

"""
    unique_indices(vec::AbstractVector)::Vector{Int}

Return the indices of the first occurrence of each unique element in the input vector `vec`
"""
function unique_indices(vec::AbstractVector)::Vector{Int}
    uniq_elements = unique(vec)
    return [findfirst(==(elem), vec) for elem in uniq_elements]
end

"""
    getindex_circular(vec::AbstractVector{T}, idx::Int)::T where {T}

Return the element of the vector `vec` at the position `idx`.

If `idx` is beyond the length of `vec` or less than 1, it wraps around in a circular manner.
"""
function getindex_circular(vec::AbstractVector{T}, idx::Int)::T where {T}
    cidx = mod(idx - 1, length(vec)) + 1
    return vec[cidx]
end

"""
    pack_grid_gradients(x::AbstractVector{T}, y::AbstractVector{T}; n_points::Int=length(x), l::Float64=1E-2) where {T<:Float64}

Returns grid between `minimum(x)` and `maximum(x)` with `n_points` points positioned to
sample `y(x)` in such a way to pack more points where gradients are greates.

`l` controls how much the adaptive gradiant sampling should approach linear sampling.
"""
function pack_grid_gradients(x::AbstractVector{T}, y::AbstractVector{T}; n_points::Int=length(x), l::Float64=1E-2) where {T<:Float64}
    tmp = abs.(gradient(x, y))
    m, M = extrema(tmp)
    tmp .+= (M - m) * l
    cumsum!(tmp, tmp)
    tmp .-= tmp[1]
    tmp ./= tmp[end]
    return interp1d(tmp, x).(range(0.0, 1.0, n_points))
end

"""
    chunk_indices(dims::Tuple{Vararg{Int}}, N::Int)

Split the indices of an array with dimensions `dims` into `N` chunks of similar size.
Each chunk is a generator of `CartesianIndex` objects.

# Arguments

  - `dims::Tuple{Vararg{Int}}`: A tuple specifying the dimensions of the array. For a 2D array, this would be `(rows, cols)`.
  - `N::Int`: The number of chunks to split the indices into.

# Returns

  - Vector with chunks of cartesian indices. Each chunk can be iterated over to get the individual `CartesianIndex` objects.
"""
function chunk_indices(dims::Tuple{Vararg{Int}}, N::Int)
    # Total number of elements
    total_elements = prod(dims)

    # Calculate chunk size
    chunk_size, remainder = divrem(total_elements, N)

    # Split indices into N chunks
    chunks = []
    start_idx = 1
    for i in 1:N
        end_idx = start_idx + chunk_size - 1
        # Distribute the remainder among the first few chunks
        if i <= remainder
            end_idx += 1
        end
        chunk_1d = start_idx:end_idx
        chunk_multi = (CartesianIndices(dims)[i] for i in chunk_1d)
        push!(chunks, chunk_multi)
        start_idx = end_idx + 1
    end

    return chunks
end

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

"""
    split_long_segments(R::AbstractVector{T},Z::AbstractVector{T},max_length::Float64) where {T<:Real}

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
            n = Int(ceil(d / max_length)) + 1
            rr = range(r1, r2, n)
            zz = range(z1, z2, n)
            append!(Rout, collect(rr)[2:end])
            append!(Zout, collect(zz)[2:end])
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
    L = sum(sqrt.(diff(R) .^ 2.0 + diff(Z) .^ 2.0))
    max_length = L / n_points
    return split_long_segments(R, Z, max_length)
end

"""
    is_open_polygon(R::AbstractVector{T}, Z::AbstractVector{T})::Bool where {T<:Real}

Determine if a polygon, defined by separate vectors for R and Z coordinates, is open.
"""
function is_open_polygon(R::AbstractVector{T}, Z::AbstractVector{T})::Bool where {T<:Real}
    return R[1] != R[end] || Z[1] != Z[end]
end

function is_open_polygon(vertices::AbstractVector)::Bool
    r1 = vertices[1][1]
    z1 = vertices[1][2]
    rend = vertices[end][1]
    zend = vertices[end][2]
    return r1 != rend || z1 != zend
end

"""
    is_closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T})::Bool where {T<:Real}

Determine if a polygon, defined by separate vectors for R and Z coordinates, is closed.
"""
function is_closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T})::Bool where {T<:Real}
    return !is_open_polygon(R, Z)
end

function is_closed_polygon(vertices::AbstractVector)::Bool
    return !is_open_polygon(vertices)
end

"""
    open_polygon(R::AbstractVector{T}, Z::AbstractVector{T}) where {T<:Real}

Convert a closed polygon into an open polygon by removing the last vertex if it is the same as the first.

Returns a named tuple containing the status of the polygon (was_closed, was_open) and the modified R and Z vectors.
"""
function open_polygon(R::AbstractVector{T}, Z::AbstractVector{T}) where {T<:Real}
    if is_open_polygon(R, Z)
        was_closed = false
    else
        was_closed = true
        R = R[1:end-1]
        Z = Z[1:end-1]
    end
    return (was_closed=was_closed, was_open=!was_closed, R=R, Z=Z, r=R, z=Z)
end

"""
    closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T}) where {T<:Real}

Convert an open polygon into a closed polygon by adding the first vertex to the end if it is not already closed.

Returns a named tuple containing the status of the polygon (was_closed, was_open) and the modified R and Z vectors.
"""
function closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T}) where {T<:Real}
    if is_open_polygon(R, Z)
        was_closed = false
        R = [R; R[1]]
        Z = [Z; Z[1]]
    else
        was_closed = true
    end
    return (was_closed=was_closed, was_open=!was_closed, R=R, Z=Z, r=R, z=Z)
end

"""
    closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T}, closed::Bool) where {T<:Real}

Returns a closed polygon depending on `closed`
"""
function closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T}, closed::Bool) where {T<:Real}
    if is_open_polygon(R, Z) && closed
        was_closed = false
        R = [R; R[1]]
        Z = [Z; Z[1]]
    else
        was_closed = true
    end
    return (was_closed=was_closed, was_open=!was_closed, R=R, Z=Z, r=R, z=Z)
end

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

"""
    is_clockwise(r::AbstractVector{T}, z::AbstractVector{T})::Bool where {T<:Real}

Returns true/false if polygon is defined clockwise
"""
function is_clockwise(r::AbstractVector{T}, z::AbstractVector{T})::Bool where {T<:Real}
    # Check if the vectors are of the same length
    if length(r) != length(z)
        error("Vectors must be of the same length")
    end

    area = zero(T)
    n = length(r)

    for i in 1:n-1
        area += (r[i] * z[i+1] - r[i+1] * z[i])
    end
    area += (r[n] * z[1] - r[1] * z[n])

    return area < 0
end

"""
    is_counterclockwise(r::AbstractVector{T}, z::AbstractVector{T})::Bool where {T<:Real}

Returns true/false if polygon is defined counterclockwise
"""
function is_counterclockwise(r::AbstractVector{T}, z::AbstractVector{T})::Bool where {T<:Real}
    return !is_clockwise(r, z)
end

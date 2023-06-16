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
    tmp = diff(vector)
    if !(1 - sum(abs.(tmp .- tmp[1])) / length(vector) ≈ 1.0)
        error("to_range requires vector data to be equally spaced")
    end
    return range(vector[1], vector[end], length=length(vector))
end

function gradient(arr::AbstractVector; method::Symbol=:central)
    return gradient(1:length(arr), arr; method)
end

"""
    gradient(coord::AbstractVector, arr::AbstractVector; method::Symbol=:central)

Finite difference method of the gradient: [:central, :backward, :forward]

The returned gradient hence has the same shape as the input array. https://numpy.org/doc/stable/reference/generated/numpy.gradient.html

For central difference, the gradient is computed using second order accurate central differences in the interior points and first order accurate one-sides (forward or backward) differences at the boundaries
"""
function gradient(coord::AbstractVector, arr::AbstractVector; method::Symbol=:central)
    if length(coord) != length(arr)
        error("The length of your coord (length = $(length(coord))) is not equal to the length of your arr (length = $(length(arr)))")
    end

    np = size(arr)[1]
    out = similar(arr)
    dcoord = diff(coord)

    # Forward difference at the beginning
    out[1] = (arr[2] - arr[1]) / dcoord[1]

    # Central difference in interior using numpy method
    if method == :central
        for p = 2:np-1
            dp1 = dcoord[p-1]
            dp2 = dcoord[p]
            a = -dp2 / (dp1 * (dp1 + dp2))
            b = (dp2 - dp1) / (dp1 * dp2)
            c = dp1 / (dp2 * (dp1 + dp2))
            out[p] = a * arr[p-1] + b * arr[p] + c * arr[p+1]
        end
    elseif method == :backward
        for p = 2:np-1
            out[p] = (arr[p] - arr[p-1]) / dcoord[p]
        end
    elseif method == :forward
        for p = 2:np-1
            out[p] = (arr[p+1] - arr[p]) / dcoord[p]
        end
    else
        error("difference method $(method) doesn't exist in gradient function")
    end

    # backward difference at the end
    out[end] = (arr[end] - arr[end-1]) / dcoord[end]

    return out
end

function gradient(arr::Matrix; method::Symbol=:central)
    return gradient(1:size(arr)[1], 1:size(arr)[2], arr; method)
end

"""
    gradient(coord1::AbstractVector, coord2::AbstractVector, arr::Matrix; method::Symbol=:central, dim::Int=0)

Finite difference method of the gradient: [:central, :backward, :forward]
applied to a matrix
on both dimensions (dim=0) or only (dim=1) or (dim=2)
"""
function gradient(coord1::AbstractVector, coord2::AbstractVector, arr::Matrix; method::Symbol=:central, dim::Int=0)
    if dim ∈ (0, 1)
        d1 = hcat(map(x -> gradient(coord1, x; method), eachcol(arr))...)
        if dim == 1
            return d1
        end
    end
    if dim ∈ (0, 2)
        d2 = transpose(hcat(map(x -> gradient(coord2, x; method), eachrow(arr))...))
        if dim == 2
            return d2
        end
    end
    return d1, d2
end

"""
    centraldiff(y::AbstractVector{<:Real})

Calculates central difference of a vector assuming that the data is equally-spaced
"""
function centraldiff(y::AbstractVector{<:T}) where {T<:Real}
    n = length(y)
    a = similar(y, T, n)
    a[1] = (y[2] - y[1]) / 2
    a[n] = (y[n] - y[n-1]) / 2
    for i in 2:n-1
        a[i] = (y[i+1] - y[i-1]) / 2
    end
    return a
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
    dy = diff(y)
    dx = diff(x)
    x_shift = view(x, 2:length(x))
    y_shift = view(y, 2:length(y))
    x0 = (x_shift .+ x[1:end-1]) .* 0.5
    y0 = (y_shift .+ y[1:end-1]) .* 0.5
    A = sum(dy .* x0)
    x_c = -sum(dx .* y0 .* x0) ./ A
    y_c = sum(dy .* x0 .* y0) ./ A
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
    intersection_angles(path1_r::T, path1_z::T, path2_r::T, path2_z::T, intersection_indexes::Vector{Tuple{Int, Int}}) where {T<:AbstractVector{<:Real}}

returns angles of intersections between two paths and intersection_indexes given by intersection() function
"""
function intersection_angles(path1_r::T, path1_z::T, path2_r::T, path2_z::T, intersection_indexes::Vector{Tuple{Int,Int}}) where {T<:AbstractVector{<:Real}}
    angles = Float64[]
    for index in intersection_indexes
        path1_p1 = [path1_r[index[1]], path1_z[index[1]]]
        path1_p2 = [path1_r[index[1]+1], path1_z[index[1]+1]]
        path2_p1 = [path2_r[index[2]], path2_z[index[2]]]
        path2_p2 = [path2_r[index[2]+1], path2_z[index[2]+1]]
        angle = mod(IMAS.angle_between_two_vectors(path1_p1, path1_p2, path2_p1, path2_p2), π)
        if angle > (π / 2.0)
            angle = π - angle
        end
        push!(angles, angle)
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

    for k1 = 1:(length(l1_x)-1)
        s1_s = StaticArrays.@SVector [l1_x[k1], l1_y[k1]]
        s1_e = StaticArrays.@SVector [l1_x[k1+1], l1_y[k1+1]]
        for k2 = 1:(length(l2_x)-1)
            s2_s = StaticArrays.@SVector [l2_x[k2], l2_y[k2]]
            s2_e = StaticArrays.@SVector [l2_x[k2+1], l2_y[k2+1]]
            crossing = _seg_intersect(s1_s, s1_e, s2_s, s2_e)
            if crossing !== nothing
                push!(indexes, (k1, k2))
                push!(crossings, (crossing[1], crossing[2]))
            end
        end
    end

    return indexes, crossings
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
    n = length(x)
    if n < 3 || n != length(y)
        error("Input arrays must have at least 3 elements and the same length")
    end

    if n == 3
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
    calculate_angle(p1::T, p2::T, p3::T) where {T<:AbstractVector{<:Real}}

Calculate the angle between three points
"""
function calculate_angle(p1::T, p2::T, p3::T) where {T<:AbstractVector{<:Real}}
    v1 = [p2[1] - p1[1], p2[2] - p1[2]]
    v2 = [p3[1] - p2[1], p3[2] - p2[2]]
    dot_product = dot(v1, v2)
    magnitude_product = norm(v1) * norm(v2)
    return acosd(dot_product / magnitude_product)
end

"""
    simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, simplification_factor::T; model::Symbol=:curvature)

Simplify 2D path by `:curvature` (Reumann-Witkam Algorithm) or `:distance` (Ramer-Douglas-Peucker) algorithms
"""
function simplify_2d_path(x::AbstractArray{T}, y::AbstractArray{T}, simplification_factor::T; model::Symbol=:curvature) where {T<:Real}
    if model == :curvature
        return rwa_simplify_2d_path(x, y, simplification_factor)
    elseif model == :distance
        return rdp_simplify_2d_path(x, y, simplification_factor)
    else
        error("simplify_2d_line model can be either :curvature or :distance")
    end
end

"""
    resample_2d_path(
        x::AbstractVector{T},
        y::AbstractVector{T};
        step::Float64=0.0,
        n_points::Integer=0,
        curvature_weight::Float64=0.0,
        method::Symbol=:cubic) where {T<:Real}

Resample 2D line with uniform stepping (or number of points)
and with option to add more points where curvature is highest
"""
function resample_2d_path(
    x::AbstractVector{T},
    y::AbstractVector{T};
    step::Float64=0.0,
    n_points::Integer=0,
    curvature_weight::Float64=0.0,
    method::Symbol=:cubic) where {T<:Real}

    s = similar(x)
    s[1] = zero(T)
    for i in 2:length(s)
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        s[i] = s[i-1] + sqrt(dx^2 + dy^2)
    end

    if curvature_weight > 0.0
        s0 = s[end]
        s ./= s0
        c = similar(s)
        c = abs.(curvature(x, y))
        c ./= maximum(c)
        s .+= (c .* curvature_weight)
        s ./= s[end]
        s .*= s0
    end

    if n_points === 0
        if step !== 0.0
            n_points = ceil(Int, s[end] / step)
        else
            n_points = length(x)
        end
    end

    t = range(s[1], s[end]; length=n_points)
    return interp1d(s, x, method).(t), interp1d(s, y, method).(t)
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
    dr1, dz1 = if pr[1] == pr[end] && pz[1] == pz[end]
        pr[end-1] - pr[end], pz[end-1] - pz[end]
    else
        0.0, 0.0
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
    calc_z(x::AbstractVector{<:Real}, f::AbstractVector{<:Real})

Returns the gradient scale lengths of vector f on x

NOTE: negative inverse scale length for typical density/temperature profiles
"""
function calc_z(x::AbstractVector{<:Real}, f::AbstractVector{<:Real})
    f[findall(ff -> (ff < 1e-32), f)] .= 1e-32
    itp = DataInterpolations.CubicSpline(f, x)
    g = [DataInterpolations.derivative(itp, x0) for x0 in x]
    return g ./ f
end

"""
    integ_z(rho::AbstractVector{<:Real}, z_profile::AbstractVector{<:Real}, bc::Real)

Backward integration of inverse scale length vector with given edge boundary condition
"""
function integ_z(rho::AbstractVector{<:Real}, z_profile::AbstractVector{<:Real}, bc::Real)
    f = interp1d(rho, z_profile, :quadratic)
    profile_new = similar(rho)
    profile_new[end] = bc
    for i in length(rho)-1:-1:1
        a = rho[i]
        fa = z_profile[i]
        b = rho[i+1]
        fb = z_profile[i+1]
        simpson_integral = (b - a) / 6 * (fa + 4 * f((a + b) * 0.5) + fb)
        profile_new[i] = profile_new[i+1] * exp(simpson_integral)
    end
    return profile_new
end

"""
    angle_between_two_vectors(
        v1_p1::Vector{T},
        v1_p2::Vector{T},
        v2_p1::Vector{T},
        v2_p2::Vector{T}) where {T<:Real}

Returns angle in radiants between two vectors defined by their start and end points
"""
function angle_between_two_vectors(
    v1_p1::Vector{T},
    v1_p2::Vector{T},
    v2_p1::Vector{T},
    v2_p2::Vector{T}) where {T<:Real}

    v1_x = v1_p2[1] - v1_p1[1]
    v1_y = v1_p2[2] - v1_p1[2]
    v2_x = v2_p2[1] - v2_p1[1]
    v2_y = v2_p2[2] - v2_p1[2]

    return acos((v1_x * v2_x + v1_y * v2_y) / (sqrt(v1_x^2 + v1_y^2) * sqrt(v2_x^2 + v2_y^2)))
end

function unique_indices(arr)
    uniq_elements = unique(arr)
    return [findfirst(==(elem), arr) for elem in uniq_elements]
end
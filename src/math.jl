import LinearAlgebra
import StaticArrays

"""
    norm01(x::Vector{<:Real})

Normalize a vector so that the first item in the array is 0 and the last one is 1
This is handy where psi_norm should be used (and IMAS does not define a psi_norm array)
"""
function norm01(x::AbstractVector{<:Real})
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

function gradient(arr::AbstractVector)
    return gradient(1:length(arr), arr)
end

"""
    gradient(coord::AbstractVector, arr::AbstractVector)

Gradient of a vector computed using second order accurate central differences in the interior points and first order accurate one-sides (forward or backwards) differences at the boundaries
The returned gradient hence has the same shape as the input array.
https://numpy.org/doc/stable/reference/generated/numpy.gradient.html
"""
function gradient(coord::AbstractVector, arr::AbstractVector)
    np = size(arr)[1]
    out = similar(arr)
    dcoord = diff(coord)

    # Forward difference at the beginning
    out[1] = (arr[2] - arr[1]) / dcoord[1]

    # Central difference in interior using numpy method
    for p = 2:np-1
        dp1 = dcoord[p-1]
        dp2 = dcoord[p]
        a = -dp2 / (dp1 * (dp1 + dp2))
        b = (dp2 - dp1) / (dp1 * dp2)
        c = dp1 / (dp2 * (dp1 + dp2))
        out[p] = a * arr[p-1] + b * arr[p] + c * arr[p+1]
    end

    # Backwards difference at the end
    out[end] = (arr[end] - arr[end-1]) / dcoord[end]

    return out
end

function gradient(arr::Matrix)
    return gradient(1:size(arr)[1], 1:size(arr)[2], arr)
end

function gradient(coord1::AbstractVector, coord2::AbstractVector, arr::Matrix)
    d1 = hcat(map(x -> gradient(coord1, x), eachcol(arr))...)
    d2 = transpose(hcat(map(x -> gradient(coord2, x), eachrow(arr))...))
    return d1, d2
end

"""
    centraldiff(y::AbstractVector)

Calculates central difference of a vector assuming that the data is equally-spaced
"""
function centraldiff(y::AbstractVector)
    dy = diff(y) / 2
    a = [dy[1]; dy]
    a .+= [dy; dy[end]]
    return a
end

"""
    meshgrid(x1::Union{Number,AbstractVector}, x2::Union{Number,AbstractVector})

Return coordinate matrices from coordinate vectors
"""
function meshgrid(x1::Union{Number,AbstractVector}, x2::Union{Number,AbstractVector})
    return x1' .* ones(length(x2)), ones(length(x1))' .* x2
end

"""
    centroid(x::Vector{<:Real}, y::Vector{<:Real})

Calculate centroid of polygon
"""
function centroid(x::Vector{<:Real}, y::Vector{<:Real})
    dy = diff(y)
    dx = diff(x)
    x0 = (x[2:end] .+ x[1:end-1]) .* 0.5
    y0 = (y[2:end] .+ y[1:end-1]) .* 0.5
    A = sum(dy .* x0)
    x_c = -sum(dx .* y0 .* x0) ./ A
    y_c = sum(dy .* x0 .* y0) ./ A
    return x_c, y_c
end

"""
    area(x::Vector{<:Real}, y::Vector{<:Real})

Calculate area of polygon
"""
function area(x::Vector{<:Real}, y::Vector{<:Real})
    x1 = x[1:end-1]
    x2 = x[2:end]
    y1 = y[1:end-1]
    y2 = y[2:end]
    return abs.(sum(x1 .* y2) - sum(y1 .* x2)) ./ 2
end

"""
    revolution_volume(x::Vector{<:Real}, y::Vector{<:Real})

Calculate volume of polygon revolved around x=0
"""
function revolution_volume(x::Vector{<:Real}, y::Vector{<:Real})
    return area(x, y) * 2pi * centroid(x, y)[1]
end

"""
    intersection(
        l1_x::AbstractVector{<:Real},
        l1_y::AbstractVector{<:Real},
        l2_x::AbstractVector{<:Real},
        l2_y::AbstractVector{<:Real};
        as_list_of_points::Bool=true)

Intersections between two 2D paths, returns list of (x,y) intersection points
"""
function intersection(
    l1_x::AbstractVector{T},
    l1_y::AbstractVector{T},
    l2_x::AbstractVector{T},
    l2_y::AbstractVector{T};
    as_list_of_points::Bool=true,
    return_indexes::Bool=false) where {T<:Real}

    indexes = NTuple{2,Int}[]
    if as_list_of_points
        crossings = NTuple{2,T}[]
    else
        crossings_x = T[]
        crossings_y = T[]
    end

    for k1 = 1:(length(l1_x)-1)
        s1_s = StaticArrays.@SVector [l1_x[k1], l1_y[k1]]
        s1_e = StaticArrays.@SVector [l1_x[k1+1], l1_y[k1+1]]
        for k2 = 1:(length(l2_x)-1)
            s2_s = StaticArrays.@SVector [l2_x[k2], l2_y[k2]]
            s2_e = StaticArrays.@SVector [l2_x[k2+1], l2_y[k2+1]]
            crossing = _seg_intersect(s1_s, s1_e, s2_s, s2_e)
            if crossing !== nothing
                push!(indexes, (k1,k2))
                if as_list_of_points
                    push!(crossings, (crossing[1], crossing[2]))
                else
                    push!(crossings_x, crossing[1])
                    push!(crossings_y, crossing[2])
                end
            end
        end
    end
    if as_list_of_points
        if return_indexes
            return indexes, crossings
        else
            return crossings
        end
    else
        if return_indexes
            return indexes, crossings_x, crossings_y
        else
            return crossings_x, crossings_y
        end
    end
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

function _seg_intersect(a1, a2, b1, b2)
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
    resample_2d_line(x, y, step)

Resample 2D line with uniform stepping
"""
function resample_2d_line(x::Vector{T}, y::Vector{T}, step::Union{Nothing,T}=nothing) where {T<:Real}
    s = cumsum(sqrt.(diff(x) .^ 2 + diff(y) .^ 2))
    s = vcat(0.0, s)
    if step !== nothing
        n = Integer(ceil(s[end] / step))
    else
        n = length(x)
    end
    t = range(s[1], s[end]; length=n)
    return interp1d(s, x).(t), interp1d(s, y).(t)
end

"""
    minimum_distance_two_shapes(R_obj1, Z_obj1, R_obj2, Z_obj2; return_index=false)

Returns minimum distance between two shapes
"""
function minimum_distance_two_shapes(R_obj1, Z_obj1, R_obj2, Z_obj2; return_index=false)
    R_obj1, Z_obj1, R_obj2, Z_obj2 = promote(R_obj1, Z_obj1, R_obj2, Z_obj2)
    distance = Inf
    ik1 = 0
    ik2 = 0
    for k1 in 1:length(R_obj1)
        for k2 in 1:length(R_obj2)
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
    mean_distance_error_two_shapes(R_obj1, Z_obj1, R_obj2, Z_obj2, target_distance)

Returns mean error distance between two shapes and a target distance
"""
function mean_distance_error_two_shapes(R_obj1, Z_obj1, R_obj2, Z_obj2, target_distance)
    R_obj1, Z_obj1, R_obj2, Z_obj2 = promote(R_obj1, Z_obj1, R_obj2, Z_obj2)
    mean_distance_error = 0.0
    for k1 in 1:length(R_obj1)
        distance = Inf
        for k2 in 1:length(R_obj2)
            @inbounds d = (R_obj1[k1] - R_obj2[k2])^2 + (Z_obj1[k1] - Z_obj2[k2])^2
            if distance > d
                distance = d
            end
        end
        mean_distance_error += (distance - target_distance)^2
    end
    return sqrt(mean_distance_error) / length(R_obj1)
end
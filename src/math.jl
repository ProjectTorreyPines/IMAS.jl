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
    centraldiff(y::AbstractVector{<:Real})

Calculates central difference of a vector assuming that the data is equally-spaced
"""
function centraldiff(y::AbstractVector{<:Real})
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
    centroid(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}

Calculate centroid of polygon
"""
function centroid(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}
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
    area(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}

Calculate area of polygon
"""
function area(x::AbstractVector{<:T}, y::AbstractVector{<:T}) where {T<:Real}
    x1 = x[1:end-1]
    x2 = x[2:end]
    y1 = y[1:end-1]
    y2 = y[2:end]
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
                push!(indexes, (k1, k2))
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
    resample_2d_line(x::Vector{T}, y::Vector{T}; step::Union{Nothing,T}=nothing, n_points=::Union{Nothing,Integer}=nothing)

Resample 2D line with uniform stepping
"""
function resample_2d_line(
    x::AbstractVector{T},
    y::AbstractVector{T};
    step::Union{Nothing,T}=nothing,
    n_points::Union{Nothing,Integer}=nothing) where {T<:Real}

    s = cumsum(sqrt.(diff(x) .^ 2 + diff(y) .^ 2))
    s = vcat(0.0, s)
    if n_points === nothing
        if step !== nothing
            n_points = Integer(ceil(s[end] / step))
        else
            n_points = length(x)
        end
    end
    t = range(s[1], s[end]; length=n_points)
    return interp1d(s, x).(t), interp1d(s, y).(t)
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
    mean_distance_error_two_shapes(
        R_obj1::AbstractVector{<:T},
        Z_obj1::AbstractVector{<:T},
        R_obj2::AbstractVector{<:T},
        Z_obj2::AbstractVector{<:T},
        target_distance::T) where {T<:Real}

Returns mean error distance between two shapes and a target distance
"""
function mean_distance_error_two_shapes(
    R_obj1::AbstractVector{<:T},
    Z_obj1::AbstractVector{<:T},
    R_obj2::AbstractVector{<:T},
    Z_obj2::AbstractVector{<:T},
    target_distance::T) where {T<:Real}

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

"""
    curvature(pr::AbstractVector{<:T}, pz::AbstractVector{<:T}) where {T<:Real}

Returns curvature of a circular 2D line
"""
function curvature(pr::AbstractVector{<:T}, pz::AbstractVector{<:T}) where {T<:Real}
    @assert (pr[1] == pr[end]) & (pz[1] == pz[end]) "curvature: 1st and last point must be the same"
    dr = diff(vcat(pr[end-1], pr, pr[2]))
    dz = diff(vcat(pz[end-1], pz, pz[2]))
    a = sqrt.(dr .^ 2.0 .+ dz .^ 2.0)
    dr1 = dr[1:end-1] ./ a[1:end-1]
    dz1 = dz[1:end-1] ./ a[1:end-1]
    dr2 = dr[2:end] ./ a[2:end]
    dz2 = dz[2:end] ./ a[2:end]
    return dr1 .* dz2 .- dr2 .* dz1
end

"""
Returns determinant of 2x2 matrix
"""
function det(a, b)
    return a[1] * b[2] - a[2] * b[1]
end

"""
    line_intersection(line1::Array{Float64, 2}, line2::Array{Float64, 2})

Returns the intersection of two lines. 
Each line is defined by two x,y points. 
"""
function line_intersection(line1::Array{Float64, 2}, line2::Array{Float64, 2})
    xdiff = (line1[1,1] - line1[2,1], line2[1,1] - line2[2,1])
    ydiff = (line1[1,2] - line1[2,2], line2[1,2] - line2[2,2])

    xmax = round(max(line2[1,1],line2[2,1]), digits=9)
    xmin = round(min(line2[1,1],line2[2,1]), digits=9)
    ymax = round(max(line2[1,2],line2[2,2]), digits=9)
    ymin = round(min(line2[1,2],line2[2,2]), digits=9)

    seg_xmax = round(max(line1[1,1],line1[2,1]), digits=9)
    seg_xmin = round(min(line1[1,1],line1[2,1]), digits=9)
    seg_ymax = round(max(line1[1,2],line1[2,2]), digits=9)
    seg_ymin = round(min(line1[1,2],line1[2,2]), digits=9)

    div = det(xdiff, ydiff)
    if div == 0
       return nothing
    end
    d = (det(line1[1,:],line1[2,:]), det(line2[1,:],line2[2,:]))
    x = round(det(d, xdiff) / div, digits=9)
    y = round(det(d, ydiff) / div, digits=9)
    if x==xmax && x==xmin && y<ymax && y>ymin && x<=seg_xmax && x>=seg_xmin && y<=seg_ymax && y>=seg_ymin
        return [x y]
    elseif x<xmax && x>xmin && y==ymax && y==ymin && x<=seg_xmax && x>=seg_xmin &&y<=seg_ymax && y>=seg_ymin
        return [x y]
    elseif x<=xmax && x>xmin && y<=ymax && y>ymin && x<=seg_xmax && x>=seg_xmin && y<=seg_ymax && y>=seg_ymin
        return [x y]
    else
        return nothing
    end
end


"""
Iterates through two lists of coordinates, finds intersection points between
those two coordinates.
"""
function layer_intersections(coord_list1, coord_list2)
    intersections=[]
    line1_already_done=false
    for (coord1_index, coord1) in enumerate(coord_list1[:,1])
        if coord1_index == length(coord_list1[:,1]) || length(coord_list1[:,1]) <= 2
            if line1_already_done==true
                continue
            end
            line1=transpose(hcat(coord_list1[1,:], coord_list1[length(coord_list1[:,1]),:]))
            line1_already_done=true
        else
            line1 = transpose(hcat(coord_list1[coord1_index,:], coord_list1[coord1_index+1,:]))
        end
        for (coord2_index, coord2) in enumerate(coord_list2[:,1])
            if coord2_index == length(coord_list2[:,1]) || length(coord_list2[:,1]) <= 2
                line2=transpose(hcat(coord_list2[1,:], coord_list2[length(coord_list2[:,1]),:]))
            else
                line2 = transpose(hcat(coord_list2[coord2_index,:], coord_list2[coord2_index+1,:]))
            end
            result_point=line_intersection(line1, line2)
            if !(result_point === nothing)
                intersections=push!(intersections,result_point)
            end
        end
    end
    return intersections
end

"""
Finds distance between two points.
"""
function find_distance(point1, point2)
    if length(point1)==0 || length(point2)==0
        return nothing
    else
        return sqrt((point1[1]-point2[1])^2 + (point1[2]-point2[2])^2)
    end
end
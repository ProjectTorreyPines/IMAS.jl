import Interpolations
import LinearAlgebra
import StaticArrays

"""
    common_base_string(s1::String, s2::String)::Vector{String}

given two strings it returns a tuple of 3 strings that is the common initial part, and then the remaining parts
"""
function common_base_string(s1::String, s2::String)::Tuple{String,String,String}
    index = nothing
    for k = 1:min(length(s1), length(s2))
        sub = SubString(s2, 1, k)
        if startswith(s1, sub)
            index = k
        end
    end
    if index === nothing
        return "", s1, s2
    else
        return string(SubString(s1, 1, index)), string(SubString(s1, index + 1, length(s1))), string(SubString(s2, index + 1, length(s2)))
    end
end

"""
    ismissing(ids::IDS, leaf)::Bool

returns true/false if field is missing in IDS
"""
function Base.ismissing(ids::IDS, field::Symbol)::Bool
    try
        getproperty(ids, field)
        return false
    catch
        return true
    end
end

function Base.ismissing(ids::IDSvector, field::Int)::Bool
    return length(ids) < field
end

function Base.ismissing(ids::Union{IDS,IDSvector}, field::Vector)::Bool
    if length(field) == 1
        return ismissing(ids, field[1])
    end
    if typeof(field[1]) <: Integer
        if field[1] <= length(ids)
            return ismissing(ids[field[1]], field[2:end])
        else
            return true
        end
    else
        return ismissing(getfield(ids, field[1]), field[2:end])
    end
end

function is_missing(ids::IDS, field::Symbol)::Bool
    @warn "Use of `IMAS.is_missing(ids, field)` is deprecated. Use `IMAS.ismissing(ids, field)`` instead."
    return ismissing(ids, field)
end

"""
    iscallable(f)

returns true if argument is callable
"""
function iscallable(f)
    return !isempty(methods(f))
end

"""
    norm01(x::Vector{T} where T<:Real)

Normalize a vector so that the first item in the array is 0 and the last one is 1
This is handy where psi_norm should be used (and IMAS does not define a psi_norm array)
"""
function norm01(x::AbstractVector{T} where {T<:Real})
    return (x .- x[1]) ./ (x[end] .- x[1])
end

"""
    to_range(vector::AbstractVector)

Turn a vector into a range (if possible)
"""
function to_range(vector::AbstractVector{T} where {T<:Real})
    tmp = diff(vector)
    if !(1 - sum(abs.(tmp .- tmp[1])) / length(vector) ≈ 1.0)
        error("to_range requires vector data to be equally spaced")
    end
    return range(vector[1], vector[end], length = length(vector))
end

"""
    Interpolations.CubicSplineInterpolation(x::AbstractVector, y::AbstractVector, args...; kw...)

Attempt to convert x::Vector to Range to feed to Interpolations.CubicSplineInterpolation
"""
function Interpolations.CubicSplineInterpolation(x::AbstractVector, y::AbstractVector, args...; kw...)
    if x[end] < x[1]
        return Interpolations.CubicSplineInterpolation(reverse(to_range(x)), reverse(y), args...; kw...)
    else
        return Interpolations.CubicSplineInterpolation(to_range(x), y, args...; kw...)
    end
end

"""
    Interpolations.LinearInterpolation(x::AbstractVector, y::AbstractVector, args...; kw...)

Attempt to convert x::Vector to Range to feed to Interpolations.LinearInterpolation
"""
function Interpolations.LinearInterpolation(x::AbstractVector, y::AbstractVector, args...; kw...)
    if x[end] < x[1]
        return Interpolations.LinearInterpolation(reverse(to_range(x)), reverse(y), args...; kw...)
    else
        return Interpolations.LinearInterpolation(to_range(x), y, args...; kw...)
    end
end


function interp(ids::IDS, field::Symbol; kw...)
    coord = coordinates(ids, field)
    if length(coord[:values]) > 1
        error("Cannot interpolate multi-dimensional $(f2i(ids)).$field that has coordinates $([k for k in coord[:names]])")
    end
    return interp(coord[:values][1], getproperty(ids, field); kw...)
end


"""
    interp(xs, y, scheme::Symbol=:linear, extrapolate::Symbol=:linear; kw...)

Interface to Interpolations.jl that makes it similar to Scipy.interpolate
"""
function interp(xs, y, scheme::Symbol = :linear, extrapolate::Symbol = :linear; kw...)

    if isa(xs, Union{AbstractVector,AbstractRange})
        xs = (xs,)
    end

    if extrapolate == :throw
        extrapolation_bc = Interpolations.Throw()
    elseif extrapolate == :linear
        extrapolation_bc = Interpolations.Line()
    elseif extrapolate == :flat
        extrapolation_bc = Interpolations.Flat()
    elseif extrapolate == :periodic
        extrapolation_bc = Interpolations.Periodic()
    elseif isa(extrapolate, Number)
        extrapolation_bc = extrapolate
    else
        error("interp extrapolation_bc can only be :throw, :flat, :linear, :periodic, or a number")
    end

    # Interpolate.jl does not handle arrays of length one
    if length(size(y)) == 1 && size(y)[1] == 1
        if extrapolate != :throw
            xs = ([xs[1][1] - 1, xs[1][1] + 1],)
            y = [y[1], y[1]]
            scheme = :constant
        end
    end

    if scheme == :constant
        itp = Interpolations.ConstantInterpolation(xs, y; extrapolation_bc = extrapolation_bc)
    elseif scheme == :linear
        itp = Interpolations.LinearInterpolation(xs, y; extrapolation_bc = extrapolation_bc)
    elseif scheme == :cubic
        itp = Interpolations.CubicSplineInterpolation(xs, y; extrapolation_bc = extrapolation_bc)
    else
        error("interp scheme can only be :constant, :linear, or :cubic ")
    end

    return itp
end

"""
    gradient(arr::AbstractVector, coord=1:length(arr))

Gradient of a vector computed using second order accurate central differences in the interior points and first order accurate one-sides (forward or backwards) differences at the boundaries
The returned gradient hence has the same shape as the input array.
https://numpy.org/doc/stable/reference/generated/numpy.gradient.html
"""
function gradient(arr::AbstractVector, coord = 1:length(arr))
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

function gradient(arr::Matrix, coord1 = 1:size(arr)[1], coord2 = 1:size(arr)[2])
    d1 = hcat(map(x -> gradient(x, coord1), eachcol(arr))...)
    d2 = transpose(hcat(map(x -> gradient(x, coord2), eachrow(arr))...))
    return d1, d2
end

"""
    meshgrid(x1::Union{Number,AbstractVector}, x2::Union{Number,AbstractVector})

Return coordinate matrices from coordinate vectors
"""
function meshgrid(x1::Union{Number,AbstractVector}, x2::Union{Number,AbstractVector})
    return x1' .* ones(length(x2)), ones(length(x1))' .* x2
end

"""
intersection(l1_x::AbstractVector{T},
             l1_y::AbstractVector{T},
             l2_x::AbstractVector{T},
             l2_y::AbstractVector{T};
             as_list_of_points::Bool=true) where T

Intersections between two 2D paths, returns list of (x,y) intersection points
"""
function intersection(
    l1_x::AbstractVector{T},
    l1_y::AbstractVector{T},
    l2_x::AbstractVector{T},
    l2_y::AbstractVector{T};
    as_list_of_points::Bool = true,
) where {T}
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
        return crossings
    else
        return crossings_x, crossings_y
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
    same_length_vectors(args...)

Returns scalars and vectors as vectors of the same lengths
For example:

    same_length_vectors(1, [2], [3,3,6], [4,4,4,4,4,4])
    
    4-element Vector{Vector{Int64}}:
    [1, 1, 1, 1, 1, 1]
    [2, 2, 2, 2, 2, 2]
    [3, 3, 6, 3, 3, 6]
    [4, 4, 4, 4, 4, 4]
"""
function same_length_vectors(args...)
    n = maximum(map(length, args))
    args = collect(map(x -> isa(x, Vector) ? x : [x], args))
    args = map(x -> vcat([x for k = 1:n]...)[1:n], args)
end

diff_function(v1, v2, tol) = maximum(abs.(v1 .- v2) ./ (tol + sum(abs.(v1) .+ abs.(v2)) / 2.0 / length(v1)))


"""
    diff(ids1::IDS, ids2::IDS, path = Any[]; tol = 1E-2, plot_function = nothing, recursive = true, differences = Dict(), verbose=true)

Compares two IDSs and returns dictionary with differences
"""
function Base.diff(ids1::IDS, ids2::IDS, path = Any[]; tol = 1E-2, plot_function = nothing, recursive = true, differences = Dict(), verbose = true)
    @assert typeof(ids1) === typeof(ids2) "IDS types are different:  $(typeof(ids1))  and  $(typeof(ids2))"
    for field in sort(collect(fieldnames(typeof(ids2))))
        if field === :_parent
            continue
        end
        v1 = missing
        v2 = missing
        try
            v1 = getproperty(ids1, field)
        catch e
            if typeof(e) <: IMASmissingDataException
                v1 = missing
            else
                rethrow
            end
        end
        try
            v2 = getproperty(ids2, field)
        catch e
            if typeof(e) <: IMASmissingDataException
                v2 = missing
            else
                rethrow
            end
        end

        pathname = p2i(vcat(path, field))
        if v1 === v2 === missing
            continue
        elseif typeof(v1) != typeof(v2)
            differences[pathname] = "types  $(typeof(v1)) --  $(typeof(v2))"
        elseif typeof(v1) <: IDS
            if recursive
                diff(v1, v2, vcat(path, field); tol, plot_function, recursive, differences)
            end
        elseif typeof(v1) <: IDSvector
            if recursive
                if length(v1) != length(v2)
                    differences[pathname] = "length  $(length(v1)) --  $(length(v2))"
                else
                    for k = 1:length(v1)
                        diff(v1[k], v2[k], vcat(path, field, k); tol, plot_function, recursive, differences)
                    end
                end
            end
        elseif typeof(v1) <: AbstractArray
            if length(v1) != length(v2)
                differences[pathname] = "length  $(length(v1)) --  $(length(v2))"
            elseif diff_function(v1, v2, tol) > tol
                if (typeof(v1) <: AbstractVector) && (plot_function !== nothing)
                    display(plot_function(1:length(v1), [v1, v2], title = pathname, label = ""))
                end
                differences[pathname] = @sprintf("value %g", diff_function(v1, v2, tol))
            end
        elseif typeof(v1) <: Number
            if diff_function(v1, v2, tol) > tol
                differences[pathname] = "value  $v1 --  $v2"
            end
        elseif typeof(v1) <: Union{String,Symbol}
            if v1 != v2
                differences[pathname] = "value  $v1 --  $v2"
            end
        else
            println((pathname, typeof(v1), typeof(v2)))
            asdadas # raise error to force use to handle this explicitly
            differences[pathname] = "value  $v1 --  $v2"
        end
        if verbose && (pathname in keys(differences))
            printstyled(pathname; bold = true)
            printstyled(" ➡ "; color = :red)
            println(differences[pathname])
        end
    end
    return differences
end

"""
    resample_2d_line(x, y, step)

Resample 2D line with uniform stepping
"""
function resample_2d_line(x::Vector{T}, y::Vector{T}, step::Union{Nothing,T} = nothing) where {T<:Real}
    s = cumsum(sqrt.(gradient(x) .^ 2 + gradient(y) .^ 2))
    if step !== nothing
        n = Integer(ceil(s[end] / step))
    else
        n = length(x)
    end
    t = range(s[1], s[end]; length = n)
    return interp(s, x)(t), interp(s, y)(t)
end

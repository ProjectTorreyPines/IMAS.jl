import Interpolations

"""
    set_field_time_array(ids::IDS, field::Symbol, time_index::Integer, value::Real)

Set the value at a given time_index of a time dependent array in a given IDS
"""
function set_field_time_array(ids::IDS, field::Symbol, time_index::Integer, value::Real)
    array = missing
    try
        array = getproperty(ids, field)
    catch
        # pass
    end
    if (array === missing)
        if time_index == 1
            return setproperty!(ids, field, [value])
        else
            error("$(f2i(ids)).$(field) is empty: cannot set data at index $(time_index)")
        end
    end
    if time_index <= length(array)
        array[time_index] = value
    elseif time_index == (length(array) + 1)
        push!(array, value)
    else
        error("$(f2i(ids)).$(field) has length $(length(array)): cannot set data at index $(time_index)")
    end
    return setproperty!(ids, field, array)
end

"""
    get_time_index(ids_time_slice::IDS, time::Real)::Integer

Return index of a given time and resize array of structures accordingly if necessary
"""
function get_time_index(ids_time_slice::IDSvector, time::Real)::Integer
    ids = top(ids_time_slice)
    if is_missing(ids, :time)
        time_index = 1
        set_field_time_array(ids, :time, time_index, time)
        resize!(ids_time_slice, time_index)
    elseif time in ids.time
        time_index = findall(x -> x == time, ids.time)[1]
    elseif maximum(ids.time) < time
        time_index = length(ids.time) + 1
        set_field_time_array(ids, :time, time_index, time)
        resize!(ids_time_slice, time_index)
    else
        error("Cannot append time slice of $(f2i(ids_time_slice)) at $(time) seconds")
    end
    return time_index
end


"""
    common_base_string(s1::String, s2::String)::Vector{String}

given two strings it returns a tuple of 3 strings that is the common initial part, and then the remaining parts
"""
function common_base_string(s1::String, s2::String)::Tuple{String,String,String}
    index = nothing
    for k in 1:min(length(s1), length(s2))
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
    is_missing(ids::IDS, leaf)::Bool

returns true/false if field is missing in IDS
"""
function is_missing(ids::IDS, field)::Bool
    try
        getproperty(ids, field)
        return false
    catch
        return true
    end
end

"""
    iscallable(f)

returns true if argument is callable
"""
function iscallable(f)
    return !isempty(methods(f))
end

"""
    norm01(x::Vector{T} where T<:Real)::Vector{Real}

Normalize a vector so that the first item in the array is 0 and the last one is 1
This is handy where psi_norm should be used (and IMAS does not define a psi_norm array)
"""
function norm01(x::AbstractVector{T} where T <: Real)::Vector{Real}
    return (x .- x[1]) ./ (x[end] .- x[1])
end

"""
    to_range(vector::AbstractVector)

Turn a vector into a range (if possible)
"""
function to_range(vector::AbstractVector{T} where T <: Real)
    tmp = diff(vector)
    if ! (1 - sum(abs.(tmp .- tmp[1])) / length(vector) ≈ 1.0)
        error("to_range requires vector data to be equally spaced")
    end
    return range(vector[1], vector[end], length=length(vector))
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

"""
    interp(xs, y, scheme::Symbol=:linear, extrapolate::Symbol=:linear; kw...)

Interface to Interpolations.jl that makes it similar to Scipy.interpolate
"""
function interp(xs, y, scheme::Symbol=:linear, extrapolate::Symbol=:linear; kw...)

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
            y = [y[1],y[1]]
            scheme = :constant
        end
    end

    if scheme == :constant
        itp = Interpolations.ConstantInterpolation(xs, y; extrapolation_bc=extrapolation_bc)
    elseif scheme == :linear
        itp = Interpolations.LinearInterpolation(xs, y; extrapolation_bc=extrapolation_bc)
    elseif scheme == :cubic
        itp = Interpolations.CubicSplineInterpolation(xs, y; extrapolation_bc=extrapolation_bc)
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
function gradient(arr::AbstractVector, coord=1:length(arr))
    np = size(arr)[1]
    out = similar(arr)
    dcoord = diff(coord)

    # Forward difference at the beginning
    out[1] = (arr[2] - arr[1]) / dcoord[1]

    # Central difference in interior using numpy method
    for p in 2:np - 1
        dp1 = dcoord[p - 1]
        dp2 = dcoord[p]
        a = -dp2 / (dp1 * (dp1 + dp2))
        b = (dp2 - dp1) / (dp1 * dp2)
        c = dp1 / (dp2 * (dp1 + dp2))
        out[p] = a * arr[p - 1] + b * arr[p] + c * arr[p + 1]
    end

    # Backwards difference at the end
    out[end] = (arr[end] - arr[end - 1]) / dcoord[end]

    return out
end

function gradient(arr::Matrix, coord1=1:size(arr)[1], coord2=1:size(arr)[2])
    d1 = hcat(map(x -> gradient(x, coord1), eachcol(arr))...)
    d2 = transpose(hcat(map(x -> gradient(x, coord2), eachrow(arr))...))
    return d1, d2
end

"""
    upsample(dim1::Union{AbstractVector,AbstractRange},
        dim2::Union{AbstractVector,AbstractRange},
        matrix::AbstractMatrix,
        upsample_factor::Integer)

Interpolate a matrix onto a grid upsample_factor larger in each dimension
"""
function upsample(dim1::Union{AbstractVector,AbstractRange},
    dim2::Union{AbstractVector,AbstractRange},
    matrix::AbstractMatrix,
    upsample_factor::Integer)

    if upsample_factor == 1
        return dim1, dim2, matrix
    elseif upsample_factor < 1
        error("Cannot upsample with factor $(upsample_factor)")
    end

    r = range(dim1[1], dim1[end], length=length(dim1))
    z = range(dim2[1], dim2[end], length=length(dim2))
    matrix_interpolant = Interpolations.CubicSplineInterpolation((r, z), matrix)
    r_upsampled = range(dim1[1], dim1[end], length=length(dim1)*upsample_factor)
    z_upsampled = range(dim2[1], dim2[end], length=length(dim2)*upsample_factor)
    return r_upsampled, z_upsampled, matrix_interpolant(r_upsampled,z_upsampled)

end
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
    norm(x::Vector{T} where T<:Real)::Vector{Real}

Normalize a vector so that the first item in the array is 0 and the last one is 1
This is handy where psi_norm should be used (and IMAS does not define a psi_norm array)
"""
function norm(x::AbstractVector{T} where T <: Real)::Vector{Real}
    return (x - x[1]) / (x[end] - x[1])
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
    Interpolations.CubicSplineInterpolation(x::AbstractVector{T} where T<:Real, y::AbstractVector where T<:Real, args...; kw...)

Attempt to convert x::Vector to Range to feed to Interpolations.CubicSplineInterpolation
"""
function Interpolations.CubicSplineInterpolation(x::AbstractVector{T} where T <: Real, y::AbstractVector where T<:Real, args...; kw...)
    if x[end]<x[1]
        return Interpolations.CubicSplineInterpolation(to_range(x[end:-1:1]), y[end:-1:1], args...; kw...)
    else
        return Interpolations.CubicSplineInterpolation(to_range(x), y, args...; kw...)
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

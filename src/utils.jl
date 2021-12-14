import Interpolations


"""
    set_timedep_value!(time_ids::IDS, value_ids::IDS, value_symbol::Symbol, time0::Real, value0::Real; time_symbol::Symbol=:time)

Set the value of a time-depentent quantity
"""
function set_timedep_value!(time_ids::IDS, value_ids::IDS, value_symbol::Symbol, time0::Real, value0::Real; time_symbol::Symbol=:time)
    if is_missing(value_ids, value_symbol) || is_missing(time_ids, time_symbol)
        setproperty!(time_ids, time_symbol, [time0])
        setproperty!(value_ids, value_symbol, [value0])
    else
        time = getproperty(time_ids, time_symbol)
        value = getproperty(value_ids, value_symbol)
        if any(time .== time0)
            time_index=argmin(abs.(time.-time0))
            value[time_index]=value0
        elseif time0 > maximum(time)
            push!(time,time0)
            push!(value,value0)
        elseif time0 < minimum(time)
            pushfirst!(time,time0)
            pushfirst!(value,value0)
        else
            for (k,t) in enumerate(time)
                if t<time0
                    continue
                else
                    setproperty!(time_ids, time_symbol, vcat(time[1:k-1],time0,time[k:end]))
                    setproperty!(value_ids, value_symbol, vcat(value[1:k-1],value0,value[k:end]))
                    break
                end
            end
        end
    end
    return getproperty(time_ids, time_symbol), getproperty(value_ids, value_symbol)
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
    meshgrid(x1::Union{Number,AbstractVector}, x2::Union{Number,AbstractVector})

Return coordinate matrices from coordinate vectors
"""
function meshgrid(x1::Union{Number,AbstractVector}, x2::Union{Number,AbstractVector})
    return x1' .* ones(length(x2)), ones(length(x1))' .* x2
end
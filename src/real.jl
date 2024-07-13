# =============== #
# Measurements.jl #
# =============== #
import Ratios
using Measurements

function Base.Int(x::Measurement)
    return Int(x.val)
end

function Base.convert(::Type{Measurement{T}}, x::Ratios.SimpleRatio{S}) where {T<:AbstractFloat,S}
    return x.num / x.den ± 0.0
end

function Base.unsafe_trunc(::Type{Int64}, x::Measurement{Float64})
    return Int(x.val)
end

function Base.convert(t::Type{T}, x::Measurement{T}) where {T<:AbstractFloat}
    if x.err == 0.0
        return convert(t, x.val)
    else
        error("Conversion of $(Base.summary(x)) $x to $t cannot be automatically handled")
    end
end

function force_float(x::Real)
    return x
end

"""
    force_float(x::Measurement)

Returns Float64
"""
function force_float(x::Measurement)
    ## we purposly do not do it recursively since generally
    ## Measurement of Measurement is an indication of someghing going wrong
    # return force_float(x.val)
    return x.val
end

"""
    value ± uncertainty

Equivalent to Measurement.measurement(value, uncertainty)
"""
±
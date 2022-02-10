# =============== #
# Measurements.jl #
# =============== #

import Ratios
using Measurements

function Base.Int(x::Measurement)
    return Int(x.val)
end

function Base.convert(::Type{Measurement{T}}, x::Ratios.SimpleRatio{S}) where {T<:AbstractFloat, S}
    return x.num/x.den ± 0.0
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

function no_error(x::Real)
    return x
end

function no_error(x::Measurement)
    return no_error(x.val)
end

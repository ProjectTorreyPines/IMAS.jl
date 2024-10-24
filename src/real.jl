# =============== #
# Measurements.jl #
# =============== #
using Measurements

"""
    value ± uncertainty

Equivalent to Measurement.measurement(value, uncertainty)
"""
±

function IMASdd.setraw!(ids::IMASdd.IDS{Float64}, field::Symbol, v::Measurements.Measurement)
    IMASdd.setraw!(ids, field, v.val)
end

function Base.convert(::Type{Array{<:Measurements.Measurement{Float64}, N}}, a::AbstractArray{Float64, N}) where N
    return a .± 0.0
end
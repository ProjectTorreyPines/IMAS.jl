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

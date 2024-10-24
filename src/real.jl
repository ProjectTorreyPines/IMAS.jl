# =============== #
# Measurements.jl #
# =============== #
import Measurements
import Measurements: ±

"""
    convert(::Type{Array{<:Measurements.Measurement{Float64},N}}, a::AbstractArray{Float64,N}) where {N}

convert a Vector of Float64 to a Measurements
"""
function Base.convert(::Type{Array{<:Measurements.Measurement{Float64},N}}, a::AbstractArray{Float64,N}) where {N}
    return a .± 0.0
end

function Base.convert(::Type{Float64}, v::Measurements.Measurement{Float64})
    return v.val
end

"""
    val ± err

Equivalent to Measurement.measurement(val, err)
"""
±

"""
    ± IDS

Unary operator that converts an IDS to Measurements
"""
function ±(@nospecialize(ids::IDS))
    return Measurements.Measurement(ids)
end

function Measurements.Measurement(@nospecialize(ids::IDS))
    ids_new = typeof(ids).name.wrapper{Measurements.Measurement{Float64}}()
    return fill!(ids_new, ids)
end

"""
    fill!(@nospecialize(ids_new::IDS{<:Measurements.Measurement{Float64}}), @nospecialize(ids::IDS{<:Float64}), field::Symbol)

Go from Float64 to Measurements
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:Measurements.Measurement{Float64}}), @nospecialize(ids::IDS{<:Float64}), field::Symbol)
    if endswith(string(field), "__error")
        return nothing
    else
        if !(fieldtype(typeof(ids), field) <: eltype(ids))
            value = getraw(ids, field)
            setraw!(ids_new, field, value)
        else
            efield = Symbol("$(field)__error")
            val = getraw(ids, field)
            if !ismissing(ids, efield)
                err = getraw(ids, efield)
                value = val ± err
            else
                value = val .± 0.0
            end
            setraw!(ids_new, field, value)
        end
    end
    return nothing
end

"""
    fill!(@nospecialize(ids_new::IDS{<:Float64}), @nospecialize(ids::IDS{<:Measurements.Measurement{Float64}}), field::Symbol)

Go from Measurements to Float64
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:Float64}), @nospecialize(ids::IDS{<:Measurements.Measurement{Float64}}), field::Symbol)
    if endswith(string(field), "__error")
        return nothing
    else
        value = getraw(ids, field)
        if !(fieldtype(typeof(ids), field) <: eltype(ids))
            setraw!(ids_new, field, value)
        else
            setraw!(ids_new, field, value.val)
            setraw!(ids_new, Symbol("$(field)__error"), value.err)
        end
    end
    return nothing
end

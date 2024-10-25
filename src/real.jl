# =============== #
# Measurements.jl #
# =============== #
import Measurements
import Measurements: ±

"""
    convert(::Type{Array{<:Measurements.Measurement{T},N}}, a::AbstractArray{T,N}) where {N}

convert AbstractArray{T,N} to a Array{<:Measurements.Measurement{T},N}
"""
function Base.convert(::Type{Array{<:Measurements.Measurement{T},N}}, a::AbstractArray{T,N}) where {T<:Real,N}
    return a .± 0.0
end

function Base.convert(::Type{T}, v::Measurements.Measurement{T}) where {T<:Real}
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

function Measurements.Measurement(@nospecialize(ids::IDS{T})) where {T<:Real}
    ids_new = typeof(ids).name.wrapper{Measurements.Measurement{T}}()
    return fill!(ids_new, ids)
end

"""
    fill!(@nospecialize(ids_new::IDS{<:Measurements.Measurement{T}}), @nospecialize(ids::IDS{<:T}), field::Symbol) where {T<:Real}

Go from IDS{T} to IDS{Measurements.Measurement{T}}
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:Measurements.Measurement{T}}), @nospecialize(ids::IDS{<:T}), field::Symbol) where {T<:Real}
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
    fill!(@nospecialize(ids_new::IDS{<:T}), @nospecialize(ids::IDS{<:Measurements.Measurement{T}}), field::Symbol)

Go from IDS{Measurements.Measurement{T}} to IDS{T}
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:T}), @nospecialize(ids::IDS{<:Measurements.Measurement{T}}), field::Symbol) where {T<:Real}
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

# =============== #
# Measurements.jl #
# =============== #
import Measurements
import Measurements: Measurement, ±
import Ratios

"""
    convert(::Type{Array{<:Measurement{T},N}}, a::AbstractArray{T,N}) where {N}

convert AbstractArray{T,N} to a Array{<:Measurement{T},N}
"""
function Base.convert(::Type{Array{<:Measurement{T},N}}, a::AbstractArray{T,N}) where {T<:Real,N}
    return a .± 0.0
end

function Base.convert(::Type{T}, v::Measurement{T}) where {T<:Real}
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
    return Measurement(ids)
end

function Base.convert(::Type{Measurement{T}}, x::Ratios.SimpleRatio{S}) where {T<:AbstractFloat,S}
    return x.num / x.den ± 0.0
end

function Base.unsafe_trunc(::Type{Int64}, x::Measurement{T}) where {T<:Real}
    return Int(x.val)
end

function Measurement(@nospecialize(ids::IDS{T})) where {T<:Real}
    ids_new = typeof(ids).name.wrapper{Measurement{T}}()
    return fill!(ids_new, ids)
end

"""
    fill!(@nospecialize(ids_new::IDS{<:Measurement{T}}), @nospecialize(ids::IDS{<:T}), field::Symbol) where {T1<:Real,T2<:Real}

Go from IDS{T} to IDS{Measurement{T}}
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:Measurement{T1}}), @nospecialize(ids::IDS{<:T2}), field::Symbol) where {T1<:Real,T2<:Real}
    if endswith(string(field), "_σ")
        return nothing
    else
        value = getraw(ids, field)
        if field == :time || !(eltype(value) <: T2)
            setraw!(ids_new, field, value)
        else
            efield = Symbol("$(field)_σ")
            if !ismissing(ids, efield)
                error = getraw(ids, efield)
                uncer = value ± error
            else
                uncer = value .± 0.0
            end
            setraw!(ids_new, field, uncer)
        end
    end
    return nothing
end

"""
    fill!(@nospecialize(ids_new::IDS{<:T1}), @nospecialize(ids::IDS{<:Measurement{T2}}), field::Symbol) where {T1<:Real,T2<:Real}

Go from IDS{Measurement{T}} to IDS{T}
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:T1}), @nospecialize(ids::IDS{<:Measurement{T2}}), field::Symbol) where {T1<:Real,T2<:Real}
    if endswith(string(field), "_σ")
        return nothing
    else
        value = getraw(ids, field)
        if !(eltype(value) <: T2)
            setraw!(ids_new, field, value)
        else
            setraw!(ids_new, field, value.val)
            setraw!(ids_new, Symbol("$(field)_σ"), value.err)
        end
    end
    return nothing
end

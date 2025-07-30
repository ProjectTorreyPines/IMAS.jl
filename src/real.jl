document[:Real] = Symbol[]

# =============== #
# Measurements.jl #
# =============== #
import Measurements
import Measurements: Measurement, ±
import Ratios
import ForwardDiff
import ForwardDiff: Dual

function Base.convert(::Type{Array{<:Measurement{T},N}}, a::AbstractArray{T,N}) where {T<:Real,N}
    return a .± 0.0
end

function Base.convert(::Type{T}, v::Measurement{T}) where {T<:Real}
    return v.val
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
    fill!(ids_new::IDS{<:T1}, ids::IDS{<:T2}, field::Symbol) where {T1<:Measurement{<:Real},T2<:Real}

Function used to map fields in `IDS{T}` to `IDS{Measurement{T}}`
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:T1}), @nospecialize(ids::IDS{<:T2}), field::Symbol) where {T1<:Measurement{<:Real},T2<:Real}
    if endswith(string(field), "_σ")
        return nothing
    else
        value = getfield(ids, field)
        if field == :time || !(eltype(value) <: T2)
            _setproperty!(ids_new, field, value; from_cocos=internal_cocos)
        else
            efield = Symbol("$(field)_σ")
            if !ismissing(ids, efield)
                error = getfield(ids, efield)
                uncer = value .± error
            else
                uncer = value .± 0.0
            end
            _setproperty!(ids_new, field, uncer; from_cocos=internal_cocos)
        end
    end
    return nothing
end

"""
    fill!(ids_new::IDS{<:T1}, ids::IDS{<:T2}, field::Symbol) where {T1<:Real,T2<:Measurement{<:Real}}

Function used to map fields in `IDS{Measurement{T}}` to `IDS{T}`
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:T1}), @nospecialize(ids::IDS{<:T2}), field::Symbol) where {T1<:Real,T2<:Measurement{<:Real}}
    if endswith(string(field), "_σ")
        return nothing
    else
        value = getraw(ids, field)
        if eltype(value) <: T2
            if typeof(value) <: AbstractArray
                _setproperty!(ids_new, field, [v.val for v in value]; from_cocos=internal_cocos)
                _setproperty!(ids_new, Symbol("$(field)_σ"), [v.err for v in value]; from_cocos=internal_cocos)
            else
                _setproperty!(ids_new, field, value.val; from_cocos=internal_cocos)
                _setproperty!(ids_new, Symbol("$(field)_σ"), value.err; from_cocos=internal_cocos)
            end
        else
            _setproperty!(ids_new, field, value; from_cocos=internal_cocos)
        end
    end
    return nothing
end

@compat public fill!
push!(document[:Real], :fill!)

function Base.setproperty!(ids::IDS{T}, field::Symbol, value::AbstractVector{<:Measurements.Measurement}) where {T<:Float64}
    setproperty!(ids, field, [v.val for v in value])
    setproperty!(ids, Symbol("$(field)_σ"), [v.err for v in value])
    return value
end

function Base.setproperty!(ids::IDS{T}, field::Symbol, value::Measurements.Measurement) where {T<:Float64}
    setproperty!(ids, field, value.val)
    setproperty!(ids, Symbol("$(field)_σ"), value.err)
    return value
end

"""
    fill!(@nospecialize(ids_new::IDS{<:T1}), @nospecialize(ids::IDS{<:T2}), field::Symbol) where {T1<:Float64,T2<:ForwardDiff.Dual{Float64,Float64,0}}

Function used to map fields in `IDS{Dual{T,T,0}}` to `IDS{T}`
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:T1}), @nospecialize(ids::IDS{<:T2}), field::Symbol) where {T1<:Float64,T2<:ForwardDiff.Dual{Float64,Float64,0}}
    value = getraw(ids, field)
    if eltype(value) <: T2
        _setproperty!(ids_new, field, ForwardDiff.value.(value); from_cocos=internal_cocos)
    else
        _setproperty!(ids_new, field, value; from_cocos=internal_cocos)
    end
    return nothing
end

function Base.setproperty!(ids::IDS{T}, field::Symbol, value::AbstractVector{<:ForwardDiff.Dual}) where {T<:Float64}
    setproperty!(ids, field, [ForwardDiff.value(v) for v in value])
    return value
end

function Base.setproperty!(ids::IDS{T}, field::Symbol, value::ForwardDiff.Dual) where {T<:Float64}
    setproperty!(ids, field, ForwardDiff.value(value))
    return value
end

function Base.round(::Type{Int64}, x::ForwardDiff.Dual, r::RoundingMode)
    return round(Int, ForwardDiff.value(x), r)
end

# Default method for Float64 - no conversion needed
force_float64(x::Float64) = x

# Method for Dual - extract the value part
force_float64(x::ForwardDiff.Dual) = Float64(ForwardDiff.value(x))

# Method for Measurements.Measurement - extract the value part
force_float64(x::Measurements.Measurement) = Float64(Measurements.value(x))

# Fallback method for other Real types
force_float64(x::Real) = Float64(x)

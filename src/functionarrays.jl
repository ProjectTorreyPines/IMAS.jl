#= ======== =#
abstract type FDS end

mutable struct FDSvector{T} <: AbstractVector{T}
    value::Vector{T}
    _parent::Union{Nothing, WeakRef}
    function FDSvector(x::Vector{T}) where {T<:FDS}
        return new{T}(x, nothing)
    end
end

struct bla <: FDS end

function Base.getindex(x::FDSvector, i::Int64)
    x.value[i]
end

Base.size(x::FDSvector) = size(x.value)

function Base.setindex!(x::FDSvector, v, i::Int64)
    x.value[i] = v
    v._parent = x._parent
end

import Base: push!, pop!

function push!(x::FDSvector, v)
    push!(x.value, v)
end

function pop!(x::Vector{T}) where {T<:FDS}
    pop!(x.value)
end

#= ======== =#

using Interpolations

abstract type AbstractFunctionArray{T,N} <: AbstractArray{T,N} end
const AFA = AbstractFunctionArray{T,N} where {T,N}

#= ======== =#

struct NumericalFunctionVector{T} <: AbstractFunctionArray{T,1}
    domain::AbstractVector{T}
    value::AbstractVector{T}
end

function NumericalFunctionVector(value)
    NumericalFunctionVector(1:length(value), value)
end

Base.broadcastable(x::NumericalFunctionVector) = Base.broadcastable(x.value)

function Base.getindex(x::NumericalFunctionVector, i::Int64)
    x.value[i]
end

Base.size(p::NumericalFunctionVector) = size(p.domain)

function Base.setindex!(x::NumericalFunctionVector, v, i::Int64)
    x.value[i] = v
end

function (x::NumericalFunctionVector)(y)
    LinearInterpolation(x.domain, x.value, )(y)
end

#= ======== =#

struct AnalyticalFunctionVector{T} <: AbstractFunctionArray{T,1}
    domain::AbstractVector{T}
    funct::Function
end

Base.broadcastable(x::AnalyticalFunctionVector) = Base.broadcastable(x.funct(x.domain))

function Base.getindex(x::AnalyticalFunctionVector, i::Int64)
    x.funct(x.domain[i])
end

Base.size(p::AnalyticalFunctionVector) = size(p.domain)

function Base.setindex!(x::AnalyticalFunctionVector, v, i::Int64)
    error("Cannot setindex! of a AnalyticalFunctionVector")
end

function (x::AnalyticalFunctionVector)(y)
    x.funct(y)
end

#= ======== =#
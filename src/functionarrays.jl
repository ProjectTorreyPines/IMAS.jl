"""
    struct_field_type(structure::DataType, field::Symbol)

Return the typeof of a given `field` witin a `structure`
"""
function struct_field_type(structure::DataType, field::Symbol)
    names = fieldnames(structure)
    index = findfirst(isequal(field), names)
    return structure.types[index]
end

#= === =#
#  FDS  # 
#= === =#

abstract type FDS end

function Base.getproperty(fds::FDS, field::Symbol)
    return getfield(fds, field)
end

function Base.setproperty!(fds::FDS, field::Symbol, v)
    if typeof(v) <: WeakRef
        return setfield!(fds, field, v)
    end

    if ! (typeof(v) <: AbstractFunctionArray)
        target_type = typeintersect(convertsion_types, struct_field_type(typeof(fds), field))
        if target_type <: AbstractFunctionArray
            v = NumericalFunctionVector(1:length(v), v)
        end
    end

    return setfield!(fds, field, v)
end

#= ========= =#
#  FDSvector  # 
#= ========= =#

mutable struct FDSvector{T} <: AbstractVector{T}
    value::Vector{T}
    _parent::WeakRef
    function FDSvector(x::Vector{T}) where {T <: FDS}
        return new{T}(x, WeakRef(nothing))
    end
end

function Base.getindex(x::FDSvector{T}, i::Int64) where {T <: FDS}
    x.value[i]
end

function Base.size(x::FDSvector{T}) where {T <: FDS}
    size(x.value)
end

function Base.setindex!(x::FDSvector{T}, v, i::Int64) where {T <: FDS}
    x.value[i] = v
    v._parent = x._parent
end

import Base: push!, pop!

function push!(x::FDSvector{T}, v) where {T <: FDS}
    push!(x.value, v)
end

function pop!(x::Vector{T}) where {T <: FDS}
    pop!(x.value)
end

#= ===================== =#
#  AbstractFunctionArray  #
#= ===================== =#

using Interpolations

abstract type AbstractFunctionArray{T,N} <: AbstractArray{T,N} end
const AFA = AbstractFunctionArray{T,N} where {T,N}

#= ====================== =#
#  NumericalFunctionArray  #
#= ====================== =#

struct NumericalFunctionVector{T} <: AbstractFunctionArray{T,1}
    domain::AbstractVector
    value::AbstractVector{T}
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

#= ======================= =#
#  AnalyticalFunctionArray  #
#= ======================= =#

struct AnalyticalFunctionVector <: AbstractFunctionArray{Float64,1}
    domain::AbstractVector
    func::Function
end

Base.broadcastable(x::AnalyticalFunctionVector) = Base.broadcastable(x.func(x.domain))

function Base.getindex(x::AnalyticalFunctionVector, i::Int64)
    x.func(x.domain[i])
end

Base.size(p::AnalyticalFunctionVector) = size(p.domain)

function Base.setindex!(x::AnalyticalFunctionVector, v, i::Int64)
    error("Cannot setindex! of a AnalyticalFunctionVector")
end

function (x::AnalyticalFunctionVector)(y)
    x.func(y)
end

#= ======== =#
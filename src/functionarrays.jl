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
    if typeof(v) <: FDS
        setfield!(v, :_parent, WeakRef(fds))
    elseif typeof(v) <: AbstractFDArray
        setfield!(v, :_field, field)
    end

    # if the value is not an AbstractFDArray...
    if ! (typeof(v) <: AbstractFDArray)
        target_type = typeintersect(conversion_types, struct_field_type(typeof(fds), field))
        # ...but the target should be one
        if target_type <: AbstractFDArray
            # figure out the coordinates
            coords = coordinates(fds, field)
            for (k, (c_name, c_value)) in enumerate(zip(coords[:names], coords[:values]))
                # do not allow assigning data before coordinates
                if c_value === missing
                    error("Assign data to `$c_name` before assigning `$(f2(fds)).$(field)`")
                end
                # generate indexes for data that does not have coordinates
                if c_value === nothing
                    coords[:values][k] = Vector{Float64}(collect(1.0:float(size(v)[k])))
                end
            end
            coords[:values] = Vector{Vector{Float64}}(coords[:values])
            # convert value to AbstractFDArray type
            if typeof(v) <: Function
                v = AnalyticalFDVector(WeakRef(fds), field, v)
            else
                v = NumericalFDVector(WeakRef(fds), field, coords[:values], v)
            end
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
        return new{T}(x, WeakRef(missing))
    end
end

function Base.getindex(x::FDSvector{T}, i::Int64) where {T <: FDS}
    x.value[i]
end

function Base.size(x::FDSvector{T}) where {T <: FDS}
    size(x.value)
end

function Base.length(x::FDSvector{T}) where {T <: FDS}
    length(x.value)
end

function Base.setindex!(x::FDSvector{T}, v, i::Int64) where {T <: FDS}
    x.value[i] = v
    setfield!(v, :_parent, WeakRef(x))
end

import Base: push!, pop!

function push!(x::FDSvector{T}, v) where {T <: FDS}
    setfield!(v, :_parent, WeakRef(x))
    push!(x.value, v)
end

function pop!(x::Vector{T}) where {T <: FDS}
    pop!(x.value)
end

function iterate(fds::FDSvector{T}) where {T <: FDS}
    return fds[1], 2
end

function iterate(fds::FDSvector{T}, state) where {T <: FDS}
    if isempty(state)
        nothing
    else
        fds[state], state + 1
    end
end

#= ======= =#
#  FDArray  #
#= ======= =#

using Interpolations

abstract type AbstractFDArray{T,N} <: AbstractArray{T,N} end
const AbstractFDVector = AbstractFDArray{T,1} where T

function coordinates(fdv::AbstractFDVector)
    return coordinates(fdv._parent.value, fdv._field)
end

function Base.setindex!(fdv::AbstractFDVector, v, i::Int64)
    error("Cannot setindex! of a $(typeof(fdv))")
end

#= ================= =#
#  NumericalFDVector  #
#= ================= =#

struct NumericalFDVector <: AbstractFDVector{Float64}
    _parent::WeakRef
    _field::Symbol
    coord_values::Vector{Vector{Float64}}
    value::Vector{Float64}
end

function Base.broadcastable(fdv::NumericalFDVector)
    coords = coordinates(fdv)
    if coords[:values][1] === nothing
        value = fdv.value
    else
        value = LinearInterpolation(fdv.coord_values[1], fdv.value)(coords[:values][1])
    end
    return Base.broadcastable(value)
end

function Base.getindex(fdv::NumericalFDVector, i::Int64)
    coords = coordinates(fdv)
    if coords[:values][1] === nothing
        value = fdv.value
    else
        value = LinearInterpolation(fdv.coord_values[1], fdv.value)(coords[:values][1])
    end
    return value[i]
end

function Base.size(fdv::NumericalFDVector)
    coords = coordinates(fdv)
    if coords[:values][1] === nothing
        value = fdv.value
    else
        value = LinearInterpolation(fdv.coord_values[1], fdv.value)(coords[:values][1])
    end
    return size(value)
end

function (fdv::NumericalFDVector)(y)
    LinearInterpolation(fdv.coord_values[1], fdv.value)(y)
end

#= ================== =#
#  AnalyticalFDVector  #
#= ================== =#

struct AnalyticalFDVector <: AbstractFDVector{Float64}
    _parent::WeakRef
    _field::Symbol
    func::Function
end

function Base.broadcastable(fdv::AnalyticalFDVector)
    y = fdv.func(coordinates(fdv)[:values][1])
    Base.broadcastable(y)
end

function Base.getindex(fdv::AnalyticalFDVector, i::Int64)
    fdv.func(coordinates(fdv)[:values][1][i])
end

Base.size(fdv::AnalyticalFDVector) = size(coordinates(fdv)[:values][1])

function (fdv::AnalyticalFDVector)(y)
    fdv.func(y)
end

include("derived_quantities.jl")

#= ============ =#
#  IMAS DD read  #
#= ============ =#

import JSON
using Memoize

"""
    imas_load_dd(ids; imas_version=imas_version)

Read the IMAS data structures in the OMAS JSON format
"""
@memoize function imas_load_dd(ids)
    JSON.parsefile(joinpath(dirname(dirname(@__FILE__)), "data_structures", "$ids.json"))  # parse and transform data
end

"""
    imas_info(location::String)

Return information of a node in the IMAS data structure
"""
function imas_info(location::String)
    location = replace(location, r"\[[0-9]+\]$" => "[:]")
    location = replace(location, r"\[:\]$" => "")
    return imas_load_dd(split(location, ".")[1])[location]
end

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
abstract type FDSvectorElement <: FDS end

function Base.getproperty(fds::FDS, field::Symbol)
    value = getfield(fds, field)
    if value === missing
        error("$(f2fs(fds)).$(field) is missing")
    else
        return value
    end
end

function Base.setproperty!(fds::FDS, field::Symbol, v)
    if typeof(v) <: FDS
        setfield!(v, :_parent, WeakRef(fds))
    elseif typeof(v) <: AbstractFDArray
        setfield!(v, :_field, field)
    end

    # if the value is not an AbstractFDArray...
    if typeof(v) <: Function
        v = FDVector(v, WeakRef(fds), field, false)
    elseif ! (typeof(v) <: AbstractFDArray)
        target_type = typeintersect(conversion_types, struct_field_type(typeof(fds), field))
        # ...but the target should be one
        if target_type <: AbstractFDArray
            # figure out the coordinates
            coords = coordinates(fds, field)
            # do not allow assigning data before coordinates
            for (c_name, c_value) in zip(coords[:names], coords[:values])
                if c_value === missing
                    error("Assign data to `$c_name` before assigning `$(f2u(fds)).$(field)`")
                end
            end
            # convert value to FDVector
            v = FDVector(coords[:values], v, WeakRef(fds), field)
        end
    end

    return setfield!(fds, field, v)
end

"""
    fds_ancestors(fds::FDS)::Dict{Symbol,Union{Missing,FDS}}

Return dictionary with pointers to ancestors to an FDS
"""
function fds_ancestors(fds::FDS)::Dict{Symbol,Union{Missing,FDS}}
    ancestors = Dict()
    # initialize ancestors to missing
    path = f2p(fds)
    for p in path
        if typeof(p) <: String
            ancestors[Symbol(p)] = missing
        end
    end
    # traverse ancestors and assign pointers
    h = fds
    while h !== missing
        path = f2p(h)
        if typeof(path[end]) <: String
            ancestors[Symbol(path[end])] = h
        elseif path[end] != 0
            ancestors[Symbol(path[end - 1])] = h
        end
        h = h._parent.value
    end
    return ancestors
end

"""
    exec_function_in_fds_namespace(fds::FDS, func::Function, func_args)

Execute a function passing the FDS stack as arguments to the function

# Arguments
- `fds::FDS`: FDS structure
- `func::Function`: function to be executed should accept in inputs the list of symbols that are traversed to get to the `fds`
- `func_args`: arguments passed to func, in addition to the list of symbols that are traversed to get to the `fds`

# Example `func`

    function pressure(x; core_profiles, electrons, profiles_1d)
        return electrons.temperature.*electrons.density * 1.60218e-19
    end

# Other example of valid `func`, now using argument splatting

    function pressure(x; electrons, _...)
        return electrons.temperature.*electrons.density * 1.60218e-19
    end
"""
function exec_function_with_fds_ancestors(fds::FDS, func::Function, func_args)
    ancestors = fds_ancestors(fds)
    func(func_args...;ancestors...)
end

#= ========= =#
#  FDSvector  #
#= ========= =#

mutable struct FDSvector{T} <: AbstractVector{T}
    value::Vector{T}
    _parent::WeakRef
    function FDSvector(x::Vector{T}) where {T <: FDSvectorElement}
        return new{T}(x, WeakRef(missing))
    end
end

function Base.getindex(x::FDSvector{T}, i::Int64) where {T <: FDSvectorElement}
    x.value[i]
end

function Base.size(x::FDSvector{T}) where {T <: FDSvectorElement}
    size(x.value)
end

function Base.length(x::FDSvector{T}) where {T <: FDSvectorElement}
    length(x.value)
end

function Base.setindex!(x::FDSvector{T}, v, i::Int64) where {T <: FDSvectorElement}
    x.value[i] = v
    setfield!(v, :_parent, WeakRef(x))
end

function Base.push!(x::FDSvector{T}, v) where {T <: FDSvectorElement}
    setfield!(v, :_parent, WeakRef(x))
    push!(x.value, v)
end

function Base.pop!(x::FDSvector{T}) where {T <: FDSvectorElement}
    pop!(x.value)
end

function iterate(fds::FDSvector{T}) where {T <: FDSvectorElement}
    return fds[1], 2
end

function iterate(fds::FDSvector{T}, state) where {T <: FDSvectorElement}
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

#= ======== =#
#  FDVector  #
#= ======== =#

struct FDVector <: AbstractFDVector{Float64}
    func_or_value::Any
    _parent::WeakRef
    _field::Symbol
    _raw::Bool
end

function FDVector(coord_values, value, _parent::WeakRef, _field::Symbol)
    # coordinate arrays
    if coord_values[1] === nothing
        func_or_value = value
        raw = true
    # arrays that have non-sorted coordinates (eg. Z limiter Vs R limiter)
    elseif ! issorted(coord_values[1]) | ! allunique(coord_values[1])
        func_or_value = value
        raw = true
    # user defined functions
    elseif iscallable(value)
        func_or_value = value
        raw = false
    # vector that have a coordinate
    else
        if length(coord_values[1]) == 1
            func_or_value = x -> value[1]
        else
            func_or_value = LinearInterpolation(coord_values[1], value)
        end
        raw = false
    end
    return FDVector(func_or_value, _parent, _field, raw)
end

function Base.broadcastable(fdv::FDVector)
    if fdv._raw
        value = fdv.func_or_value
    elseif typeof(fdv.func_or_value) <: Function # user defined function
        value = exec_function_with_fds_ancestors(fdv._parent.value, fdv.func_or_value, coordinates(fdv)[:values])
    else # interpolation
        value = fdv.func_or_value(coordinates(fdv)[:values][1])
    end
    Base.broadcastable(value)
end

function Base.getindex(fdv::FDVector, i::Int64)
    if fdv._raw
        return fdv.func_or_value[i]
    elseif typeof(fdv.func_or_value) <: Function # user defined function
        return exec_function_with_fds_ancestors(fdv._parent.value, fdv.func_or_value, coordinates(fdv)[:values])[i]
    else # interpolation
        return fdv.func_or_value(coordinates(fdv)[:values][1])[i]
    end
end

function Base.size(fdv::FDVector)
    if fdv._raw
        return size(fdv.func_or_value)
    else
        return Tuple{Int}([size(c)[1] for c in coordinates(fdv)[:values]])
    end
end

function (fdv::FDVector)(x)
    if fdv._raw
        return LinearInterpolation(range(0.0, 1.0, length=length(fdv.func_or_value)), fdv.func_or_value)(x)
    elseif typeof(fdv.func_or_value) <: Function # user defined function
        return exec_function_with_fds_ancestors(fdv._parent.value, fdv.func_or_value, [x])
    else # interpolation
        return fdv.func_or_value(x)
    end
end

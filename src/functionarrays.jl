include("expressions.jl")

#= ============ =#
#  IMAS DD read  #
#= ============ =#

import JSON
import Memoize

"""
    imas_load_dd(ids; imas_version=imas_version)

Read the IMAS data structures in the OMAS JSON format
"""
@Memoize.memoize function imas_load_dd(ids)
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
    # raise a nice error for missing values
    if value === missing
        error("$(f2fs(fds)).$(field) is missing")
    # interpolate functions on given coordinates
    elseif typeof(value) <: Function
        x = coordinates(fds,field)[:values]
        return exec_expression_with_ancestor_args(fds, field, value, x)
    else
        return value
    end
end

"""
    Base.setproperty!(fds::FDS, field::Symbol, v; skip_non_coordinates=false)

# Arguments
- `skip_non_coordinates::Bool=false`: don't do anything for fields that not a IMAS coordinate
"""
function Base.setproperty!(fds::FDS, field::Symbol, v; skip_non_coordinates=false)
    if typeof(v) <: FDS
        setfield!(v, :_parent, WeakRef(fds))
    end
    if typeof(v) <: StepRangeLen
        v = collect(v)
    end
    if typeof(v) <: AbstractArray
        # figure out the coordinates
        coords = coordinates(fds, field)
        # do not allow assigning data before coordinates
        for (c_name, c_value) in zip(coords[:names], coords[:values])
            if c_value === missing
                if skip_non_coordinates
                    return
                else
                    error("Assign data to `$c_name` before assigning `$(f2u(fds)).$(field)`")
                end
            end
        end
    end
    setfield!(fds, field, v)
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

const expression_call_stack = String[]

"""
    exec_expression_with_ancestor_args(fds::FDS, field::Symbol, func::Function, func_args)

Execute a function passing the FDS stack as arguments to the function

# Arguments
- `fds::FDS`: FDS structure
- `field::Symbol`: Field in the `fds` that is being called 
- `func::Function`: Function to be executed should accept in inputs the list of symbols that are traversed to get to the `fds`
- `func_args`: Arguments passed to func, in addition to the list of symbols that are traversed to get to the `fds`

# Example `func`

    function pressure(x; core_profiles, electrons, profiles_1d)
        return electrons.temperature.*electrons.density * 1.60218e-19
    end

# Other example of valid `func`, now using argument splatting

    function pressure(x; electrons, _...)
        return electrons.temperature.*electrons.density * 1.60218e-19
    end
"""
function exec_expression_with_ancestor_args(fds::FDS, field::Symbol, func::Function, func_args)
    structure_name = f2u(fds)
    # keep track of recursion
    if ! (structure_name in expression_call_stack)
        push!(expression_call_stack, "$(structure_name).$(field)")
    else
        culprits = join(expression_call_stack, "\n    * ")
        error("These expressions are calling themselves recursively:\n    * $(culprits)\nAssign a numerical value to one of them to break the cycle.")
    end
    # find ancestors to this fds
    ancestors = fds_ancestors(fds)
    # execute and in all cases pop the call_stack
    try
        func(func_args...;ancestors...)
    finally
        pop!(expression_call_stack)
    end
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
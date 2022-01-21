include("expressions.jl")

#= ============ =#
#  IMAS DD read  #
#= ============ =#

import JSON
import Memoize

"""
    imas_dd_ids_names(dirs::Vector{String}=["data_structures", "data_structures_extra"])

Return list of IDS names
"""
function imas_dd_ids_names(dirs::Vector{String}=["data_structures", "data_structures_extra"])
    filenames = String[]
    for dir in dirs
        append!(filenames, [filename for filename in readdir(joinpath(dirname(dirname(@__FILE__)), dir)) if ! startswith(filename, "_") && endswith(filename, ".json")])
    end
    filenames = [replace(filename, ".json" => "") for filename in sort(filenames)]
end

"""
    imas_dd_ids(ids; imas_version=imas_version)

Read the IMAS data structures in the OMAS JSON format
"""
Memoize.@memoize function imas_dd_ids(ids)
    tmp = Dict()
    for dir in ["data_structures", "data_structures_extra"]
        filename = joinpath(dirname(dirname(@__FILE__)), dir, "$ids.json")
        if ispath(filename)
            tmp = merge(tmp, JSON.parsefile(filename))
        end
    end
    return tmp
end

"""
    imas_info(location::String)

Return information of a node in the IMAS data structure
"""
function imas_info(location::String)
    location = replace(location, r"\[[0-9]+\]$" => "[:]")
    location = replace(location, r"\[:\]$" => "")
    return imas_dd_ids(split(location, ".")[1])[location]
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
#  IDS  #
#= === =#

abstract type IDS end
abstract type IDSvectorElement <: IDS end

function Base.getproperty(ids::IDS, field::Symbol)
    value = getfield(ids, field)
    # raise a nice error for missing values
    if typeof(value) <: IDS
        return value
    elseif value === missing
        error("$(f2i(ids)).$(field) is missing")
    # interpolate functions on given coordinates
    elseif typeof(value) <: Function
        x = coordinates(ids, field)[:values]
        try
            return exec_expression_with_ancestor_args(ids, field, value, x)
        catch e
            println("Error with expression in $(f2u(ids)).$field")
            rethrow(e)
        end
    else
        return value
    end
end

"""
    Base.setproperty!(ids::IDS, field::Symbol, v; skip_non_coordinates=false)

# Arguments
- `skip_non_coordinates::Bool=false`: don't do anything for fields that not a IMAS coordinate
"""
function Base.setproperty!(ids::IDS, field::Symbol, v; skip_non_coordinates=false)
    if typeof(v) <: IDS
        setfield!(v, :_parent, WeakRef(ids))
    end
    if typeof(v) <: StepRangeLen
        v = collect(v)
    end
    if typeof(v) <: AbstractArray
        # figure out the coordinates
        coords = coordinates(ids, field)
        # do not allow assigning data before coordinates
        for (c_name, c_value) in zip(coords[:names], coords[:values])
            if c_value === missing
                if skip_non_coordinates
                    return
                else
                    error("Assign data to `$c_name` before assigning `$(f2u(ids)).$(field)`")
                end
            end
        end
    end
    try
        setfield!(ids, field, v)
    catch e
        typeof(e) <: TypeError && error("$(typeof(v)) is the wrong type for $(f2u(ids)).$(field), it should be $(fieldtype(typeof(ids), field)))")
        rethrow()
    end
end

"""
    ids_ancestors(ids::IDS)::Dict{Symbol,Union{Missing,IDS,Int}}

Return dictionary with pointers to ancestors to an IDS
"""
function ids_ancestors(ids::IDS)::Dict{Symbol,Union{Missing,IDS,Int}}
    ancestors = Dict()
    # initialize ancestors to missing
    path = f2p(ids)
    for (k, p) in enumerate(path)
        if typeof(p) <: String
            ancestors[Symbol(p)] = missing
        elseif p != 0
            ancestors[Symbol(path[k - 1] * "_index")] = missing
        end
    end
    # traverse ancestors and assign pointers
    h = ids
    while h !== missing
        path = f2p(h)
        if typeof(path[end]) <: String
            ancestors[Symbol(path[end])] = h
        elseif path[end] != 0
            ancestors[Symbol(path[end - 1])] = h
            ancestors[Symbol(path[end - 1] * "_index")] = path[end]
        end
        h = h._parent.value
    end
    return ancestors
end

const expression_call_stack = String[]

"""
    exec_expression_with_ancestor_args(ids::IDS, field::Symbol, func::Function, func_args)

Execute a function passing the IDS stack as arguments to the function

# Arguments
- `ids::IDS`: IDS structure
- `field::Symbol`: Field in the `ids` that is being called 
- `func::Function`: Function to be executed should accept in inputs the list of symbols that are traversed to get to the `ids`
- `func_args`: Arguments passed to func, in addition to the list of symbols that are traversed to get to the `ids`

# Example `func`

    function pressure(x; core_profiles, electrons, profiles_1d)
        return electrons.temperature.*electrons.density * 1.60218e-19
    end

# Other example of valid `func`, now using argument splatting

    function pressure(x; electrons, _...)
        return electrons.temperature.*electrons.density * 1.60218e-19
    end
"""
function exec_expression_with_ancestor_args(ids::IDS, field::Symbol, func::Function, func_args)
    structure_name = "$(f2u(ids)).$(field)"
    # keep track of recursion
    if ! (structure_name in expression_call_stack)
        push!(expression_call_stack, structure_name)
    else
        culprits = join(expression_call_stack, "\n    * ")
        error("These expressions are calling themselves recursively:\n    * $(culprits)\nAssign a numerical value to one of them to break the cycle.")
    end
    # find ancestors to this ids
    ancestors = ids_ancestors(ids)
    # execute and in all cases pop the call_stack
    try
        func(func_args...;ancestors...)
    finally
        pop!(expression_call_stack)
    end
end

#= ========= =#
#  IDSvector  #
#= ========= =#

mutable struct IDSvector{T} <: AbstractVector{T}
    value::Vector{T}
    _parent::WeakRef
    function IDSvector(x::Vector{T}) where {T <: IDSvectorElement}
        return new{T}(x, WeakRef(missing))
    end
end

function Base.getindex(x::IDSvector{T}, i::Int) where {T <: IDSvectorElement}
    x.value[i]
end

function Base.size(x::IDSvector{T}) where {T <: IDSvectorElement}
    size(x.value)
end

function Base.length(x::IDSvector{T}) where {T <: IDSvectorElement}
    length(x.value)
end

function Base.setindex!(x::IDSvector{T}, v, i::Int) where {T <: IDSvectorElement}
    x.value[i] = v
    setfield!(v, :_parent, WeakRef(x))
end

function Base.push!(x::IDSvector{T}, v) where {T <: IDSvectorElement}
    setfield!(v, :_parent, WeakRef(x))
    push!(x.value, v)
end

function Base.pop!(x::IDSvector{T}) where {T <: IDSvectorElement}
    pop!(x.value)
end

function iterate(ids::IDSvector{T}) where {T <: IDSvectorElement}
    return ids[1], 2
end

function iterate(ids::IDSvector{T}, state) where {T <: IDSvectorElement}
    if isempty(state)
        nothing
    else
        ids[state], state + 1
    end
end

function Base.deleteat!(x::IDSvector{T}, i::Int) where {T <: IDSvectorElement}
    return Base.deleteat!(x.value, i)
end
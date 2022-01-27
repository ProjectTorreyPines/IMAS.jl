include("expressions.jl")

struct GlobalTime <: AbstractFloat end
const τ = GlobalTime()

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
abstract type IDSvectorStaticElement <: IDSvectorElement end
abstract type IDSvectorTimeElement <: IDSvectorElement end

function Base.getproperty(ids::IDS, field::Symbol)
    value = getfield(ids, field)
    # raise a nice error for missing values
    if typeof(value) <: Union{IDS,IDSvector}
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
function _setproperty!(ids::IDS, field::Symbol, v)
    try
        setfield!(ids, field, v)
    catch e
        typeof(e) <: TypeError && error("$(typeof(v)) is the wrong type for $(f2u(ids)).$(field), it should be $(fieldtype(typeof(ids), field)))")
        rethrow()
    end
end

function Base.setproperty!(ids::IDS, field::Symbol, v; skip_non_coordinates=false)
    _setproperty!(ids, field, v)
end

function Base.setproperty!(ids::IDS, field::Symbol, v::StepRangeLen; skip_non_coordinates=false)
    v = collect(v)
    setproperty!(ids, field, v; skip_non_coordinates=skip_non_coordinates)
end

function Base.setproperty!(ids::IDS, field::Symbol, v::AbstractArray; skip_non_coordinates=false)
    # figure out the coordinates
    coords = coordinates(ids, field)
    # do not allow assigning data before coordinates
    for (c_name, c_value) in zip(coords[:names], coords[:values])
        if c_value === missing
            if skip_non_coordinates
                return
            else
                error("Can't assign data to `$(f2u(ids)).$(field)` before `$c_name`")
            end
        end
    end
    _setproperty!(ids, field, v)
end

function Base.setproperty!(ids::IDSvector{T}, field::Symbol, v::T; skip_non_coordinates=false) where {T <: IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    _setproperty!(ids, field, v)
end

function Base.setproperty!(ids::IDS, field::Symbol, v::Union{IDS,IDSvector}; skip_non_coordinates=false)
    setfield!(v, :_parent, WeakRef(ids))
    _setproperty!(ids, field, v)
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
    try
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
        func(func_args...; ancestors...)
    finally
        if length(expression_call_stack) > 0
            pop!(expression_call_stack)
        end
    end
end

#= ========= =#
#  IDSvector  #
#= ========= =#

mutable struct IDSvector{T} <: AbstractVector{T}
    _value::Vector{T}
    _parent::WeakRef
    function IDSvector(ids::Vector{T}) where {T <: IDSvectorElement}
        return new{T}(ids, WeakRef(missing))
    end
end

function Base.size(ids::IDSvector{T}) where {T <: IDSvectorElement}
    size(ids._value)
end

function Base.length(ids::IDSvector{T}) where {T <: IDSvectorElement}
    length(ids._value)
end

function Base.getindex(ids::IDSvector{T}) where {T <: IDSvectorTimeElement}
    return getindex(ids, global_time(ids))
end

function Base.getindex(ids::IDSvector{T}, time0::GlobalTime) where {T <: IDSvectorTimeElement}
    return getindex(ids, global_time(ids))
end

function Base.getindex(ids::IDSvector{T}, time0::AbstractFloat) where {T <: IDSvectorTimeElement}
    if length(ids) == 0
        ids[1]
    end
    time = time_array(ids) 
    i = argmin(abs.(time .- time0))
    ids._value[i]
end

function Base.getindex(ids::IDSvector{T}, i::Int) where {T <: IDSvectorElement}
    try
        ids._value[i]
    catch
        error("Attempt to access $(length(ids))-element $(typeof(ids)) at index [$i]. Need to `resize!(ids, $i)`")
    end
end

function Base.setindex!(ids::IDSvector{T}, v::T) where {T <: IDSvectorTimeElement}
    setindex!(ids, v, global_time(ids))
end

function Base.setindex!(ids::IDSvector{T}, v::T, time0::GlobalTime) where {T <: IDSvectorTimeElement}
    setindex!(ids, v, global_time(ids))
end

function Base.setindex!(ids::IDSvector{T}, v::T, time0::AbstractFloat) where {T <: IDSvectorTimeElement}
    time = time_array(ids)
    if time0 < minimum(time)
        pushfirst!(time, time0)
        pushfirst!(ids, v)
    elseif time0 > maximum(time)
        push!(time, time0)
        push!(ids, v)
    elseif minimum(abs.(time .- time0)) == 0
        i = argmin(abs.(time .- time0))
        # perfect match --> overwrite
        ids._value[i] = v
    else
        error("Cannot insert data at time $time0 in middle of a time array structure ranging between $(time[1]) and $(time[end])")
    end
    if hasfield(typeof(v), :time)
        v.time=time0
    end
    setfield!(v, :_parent, WeakRef(ids))
end

function Base.setindex!(ids::IDSvector{T}, v::T, i::Int) where {T <: IDSvectorElement}
    ids._value[i] = v
    setfield!(v, :_parent, WeakRef(ids))
end

function Base.push!(ids::IDSvector{T}, v::T) where {T <: IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    push!(ids._value, v)
end

function Base.pushfirst!(ids::IDSvector{T}, v::T) where {T <: IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    pushfirst!(ids._value, v)
end

function Base.pop!(ids::IDSvector{T}) where {T <: IDSvectorElement}
    pop!(ids._value)
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

function Base.insert!(ids::IDSvector{T}, i, v::T) where {T <: IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    insert!(ids._value, i, v)
end

function Base.deleteat!(ids::IDSvector{T}, i::Int) where {T <: IDSvectorElement}
    return Base.deleteat!(ids._value, i)
end

function Base.resize!(ids::IDSvector{T}) where {T <: IDSvectorTimeElement}
    return Base.resize!(ids, global_time(ids))
end

function Base.resize!(ids::IDSvector{T}, time0::GlobalTime) where {T <: IDSvectorTimeElement}
    return Base.resize!(ids, global_time(ids))
end

function Base.resize!(ids::IDSvector{T}, time0::AbstractFloat) where {T <: IDSvectorTimeElement}
    time = time_array(ids)
    if (length(ids) == 0) || (time0 > maximum(time))
        k = length(ids) + 1
        resize!(ids, k)
        push!(time, time0)
        if hasfield(typeof(ids[k]), :time)
            ids[k].time = time0
        end
    elseif time0 < maximum(time)
        error("Cannot resize structure at time $time0 for a time array structure already ranging between $(time[1]) and $(time[end])")
    end
    return ids
end

function Base.resize!(ids::IDSvector{T}, n::Int) where {T <: IDSvectorElement}
    if n > length(ids)
        for k in length(ids):n - 1
            obj = eltype(ids)()
            setfield!(obj, :_parent, WeakRef(ids))
            push!(ids._value, obj)
        end
    elseif n < length(ids)
        for k in n:length(ids) - 1
            pop!(ids._value)
        end
    end
    return ids
end

function Base.empty!(ids::IDS)
    tmp = typeof(ids)()
    for item in fieldnames(typeof(ids))
        if item != :_parent
            setproperty!(ids, item, getfield(tmp, item))
        end
        # if typeof(getfield(tmp, item)) <: Union{IDS,IDSvector}
        #     setfield!(getfield(tmp, item), :_parent, WeakRef(ids))
        # end
    end
    assign_expressions(ids)
end

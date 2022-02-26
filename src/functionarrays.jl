struct GlobalTime <: AbstractFloat end
const τ = GlobalTime()

#= ============================ =#
#  IDS and IDSvector structures  #
#= ============================ =#

abstract type IDS end
abstract type IDSvectorElement <: IDS end
abstract type IDSvectorStaticElement <: IDSvectorElement end
abstract type IDSvectorTimeElement <: IDSvectorElement end
mutable struct IDSvector{T} <: AbstractVector{T}
    _value::Vector{T}
    _parent::WeakRef
    function IDSvector(ids::Vector{T}) where {T<:IDSvectorElement}
        return new{T}(ids, WeakRef(missing))
    end
end

#= ============ =#
#  IMAS DD read  #
#= ============ =#

import JSON
import Memoize

"""
    imas_dd_ids_names(dirs::Vector{String}=["data_structures", "data_structures_extra"])

Return list of IDS names
"""
function imas_dd_ids_names(dirs::Vector{String} = ["data_structures", "data_structures_extra"])
    filenames = String[]
    for dir in dirs
        dir_filenames = readdir(joinpath(dirname(dirname(@__FILE__)), dir))
        append!(filenames, [filename for filename in dir_filenames if !startswith(filename, "_") && endswith(filename, ".json")])
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
    imas_info(ids::IDS, field::Symbol)

Return information of a filed of an IDS
"""
function imas_info(ids::Union{IDS,IDSvector}, field::Symbol)
    location = "$(f2u(ids)).$(field)"
    return imas_info(location)
end


function imas_info(ids::Type, field::Symbol)
    if ids <: IMAS.IDSvectorElement
        location = "$(IMAS._f2u(ids))[:].$(field)"
    else
        location = "$(IMAS._f2u(ids)).$(field)"
    end
    return imas_info(location)
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
#  Exceptions  #
#= === =#

struct IMASmissingDataException <: Exception
    ids::IDS
    field::Symbol
end

Base.showerror(io::IO, e::IMASmissingDataException) = print(io, "$(f2i(e.ids)).$(e.field) is missing")

abstract type IMASexpressionError <: Exception end

struct IMASexpressionRecursion <: IMASexpressionError
    call_stack::Vector{String}
end

function Base.showerror(io::IO, e::IMASexpressionRecursion)
    culprits = join(e.call_stack, "\n    * ")
    print(io, "These expressions are calling themselves recursively:\n    * $(culprits)\nAssign a numerical value to one of them to break the cycle.")
end

struct IMASbadExpression <: IMASexpressionError
    ids::IDS
    field::Symbol
    reason::String
end

Base.showerror(io::IO, e::IMASbadExpression) = print(io, "Bad expression $(f2i(e.ids)).$(e.field)\n$(e.reason)")

#= === =#
#  IDS  #
#= === =#

function Base.getproperty(ids::IDS, field::Symbol)
    value = getfield(ids, field)
    if value === missing
        # raise a nice error for missing values
        throw(IMASmissingDataException(ids, field))
    elseif typeof(value) <: Union{IDS,IDSvector}
        # nothing to do for data structures
        return value
    elseif typeof(value) <: Function
        # interpolate functions on given coordinates
        x = coordinates(ids, field)[:values]
        try
            return exec_expression_with_ancestor_args(ids, field, value, x)
        catch e
            if typeof(e) <: IMASexpressionRecursion
                rethrow(e)
            else
                throw(IMASbadExpression(ids, field, sprint(showerror, e)))
            end
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

function Base.setproperty!(ids::IDS, field::Symbol, v; skip_non_coordinates = false)
    _setproperty!(ids, field, v)
end

function Base.setproperty!(ids::IDS, field::Symbol, v::StepRangeLen; skip_non_coordinates = false)
    v = collect(v)
    setproperty!(ids, field, v; skip_non_coordinates = skip_non_coordinates)
end

function Base.setproperty!(ids::IDS, field::Symbol, v::AbstractArray; skip_non_coordinates = false)
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

function Base.setproperty!(ids::IDSvector{T}, field::Symbol, v::T; skip_non_coordinates = false) where {T<:IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    _setproperty!(ids, field, v)
end

function Base.setproperty!(ids::IDS, field::Symbol, v::Union{IDS,IDSvector}; skip_non_coordinates = false)
    setfield!(v, :_parent, WeakRef(ids))
    _setproperty!(ids, field, v)
end

function Base.deepcopy(ids::Union{IDS,IDSvector})
    ids1 = Base.deepcopy_internal(ids, Base.IdDict())
    ids1._parent = WeakRef(missing)
    return ids1
end

"""
    ids_ancestors(ids::IDS)::Dict{Symbol,Union{Missing,IDS,Int}}

Return dictionary with pointers to ancestors to an IDS
"""
function ids_ancestors(ids::IDS)::Dict{Symbol,Union{Missing,IDS,Int}}
    ancestors = Dict()
    # initialize ancestors to missing
    path = f2p(ids)
    pushfirst!(path, "dd")
    for (k, p) in enumerate(path)
        if typeof(p) <: String
            ancestors[Symbol(p)] = missing
        else
            ancestors[Symbol(path[k-1] * "_index")] = missing
        end
    end
    # traverse ancestors and assign pointers
    h = ids
    while typeof(h) <: Union{IDS,IDSvector}
        path = f2p(h)
        if typeof(path[end]) <: String
            ancestors[Symbol(path[end])] = h
        elseif path[end] != 0
            ancestors[Symbol(path[end-1])] = h
            ancestors[Symbol(path[end-1] * "_index")] = path[end]
        end
        h = h._parent.value
    end
    # for anc in keys(ancestors)
    #     println("$anc => $(typeof(ancestors[anc]))")
    # end
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
        if !(structure_name in expression_call_stack)
            push!(expression_call_stack, structure_name)
        else
            throw(IMASexpressionRecursion(copy(expression_call_stack)))
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

function Base.size(ids::IDSvector{T}) where {T<:IDSvectorElement}
    size(ids._value)
end

function Base.length(ids::IDSvector{T}) where {T<:IDSvectorElement}
    length(ids._value)
end

function Base.getindex(ids::IDSvector{T}) where {T<:IDSvectorTimeElement}
    return getindex(ids, global_time(ids))
end

function Base.getindex(ids::IDSvector{T}, time0::GlobalTime) where {T<:IDSvectorTimeElement}
    return getindex(ids, global_time(ids))
end

function Base.getindex(ids::IDSvector{T}, time0::AbstractFloat) where {T<:IDSvectorTimeElement}
    if length(ids) == 0
        ids[1]
    end
    time = time_parent(ids).time
    i = argmin(abs.(time .- time0))
    ids._value[i]
end

function Base.getindex(ids::IDSvector{T}, i::Int) where {T<:IDSvectorElement}
    try
        ids._value[i]
    catch
        error("Attempt to access $(length(ids))-element $(typeof(ids)) at index [$i]. Need to `resize!(ids, $i)`")
    end
end

function Base.setindex!(ids::IDSvector{T}, v::T) where {T<:IDSvectorTimeElement}
    setindex!(ids, v, global_time(ids))
end

function Base.setindex!(ids::IDSvector{T}, v::T, time0::GlobalTime) where {T<:IDSvectorTimeElement}
    setindex!(ids, v, global_time(ids))
end

function Base.setindex!(ids::IDSvector{T}, v::T, time0::AbstractFloat) where {T<:IDSvectorTimeElement}
    time = time_parent(ids).time
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
        v.time = time0
    end
    setfield!(v, :_parent, WeakRef(ids))
end

function Base.setindex!(ids::IDSvector{T}, v::T, i::Int) where {T<:IDSvectorElement}
    ids._value[i] = v
    setfield!(v, :_parent, WeakRef(ids))
end

function Base.push!(ids::IDSvector{T}, v::T) where {T<:IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    push!(ids._value, v)
end

function Base.pushfirst!(ids::IDSvector{T}, v::T) where {T<:IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    pushfirst!(ids._value, v)
end

function Base.pop!(ids::IDSvector{T}) where {T<:IDSvectorElement}
    pop!(ids._value)
end

function iterate(ids::IDSvector{T}) where {T<:IDSvectorElement}
    return ids[1], 2
end

function iterate(ids::IDSvector{T}, state) where {T<:IDSvectorElement}
    if isempty(state)
        nothing
    else
        ids[state], state + 1
    end
end

function Base.insert!(ids::IDSvector{T}, i, v::T) where {T<:IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    insert!(ids._value, i, v)
end

function Base.deleteat!(ids::IDSvector{T}, i::Int) where {T<:IDSvectorElement}
    return deleteat!(ids._value, i)
end

function Base.empty!(ids::IDS)
    tmp = typeof(ids)()
    for item in fieldnames(typeof(ids))
        if item != :_parent
            setproperty!(ids, item, getfield(tmp, item))
        end
    end
    assign_expressions(ids)
    return ids
end

function Base.empty!(ids::IDS, field::Symbol)
    struct_name = f2u(ids)
    if "$(struct_name).$(field)" in keys(expressions)
        return setproperty!(ids, field, expressions["$(struct_name).$(field)"])
    elseif typeof(getfield(ids, field)) <: Union{IDS,IDSvector}
        return empty!(getfield(ids, field))
    else
        setproperty!(ids, field, missing)
        return missing
    end
end

function Base.resize!(ids::IDSvector{T}) where {T<:IDSvectorTimeElement}
    return resize!(ids, global_time(ids))
end

function Base.resize!(ids::IDSvector{T}, time0::GlobalTime) where {T<:IDSvectorTimeElement}
    return resize!(ids, global_time(ids))
end

function Base.resize!(ids::IDSvector{T}, time0::AbstractFloat) where {T<:IDSvectorTimeElement}
    _time = time_parent(ids)
    time = _time.time
    if length(ids) > length(time)
        error(
            "Length [$(length(ids))] of $(p2i(f2p(ids)[1:end-1])) does not match length [$(length(time))] of $(f2i(_time)).time. Cannot resize based on global_time.",
        )
    elseif (length(ids) == 0) || (time0 > maximum(time))
        k = length(ids) + 1
        resize!(ids, k)
        if !(time0 in time)
            push!(time, time0)
        end
        if hasfield(typeof(ids[k]), :time)
            ids[k].time = time0
        end
    elseif time0 == maximum(time)
        empty!(ids[end])
        if hasfield(typeof(ids[end]), :time)
            ids[end].time = time0
        end
    elseif time0 < maximum(time)
        error("Cannot resize structure at time $time0 for a time array structure already ranging between $(time[1]) and $(time[end])")
    end
    return ids[end]
end

function Base.resize!(ids::IDSvector{T}, n::Int) where {T<:IDSvectorElement}
    if n == 0
        return empty!(ids)
    elseif n > length(ids)
        for k = length(ids):n-1
            push!(ids, eltype(ids)())
        end
    else
        if n < length(ids)
            for k = n:length(ids)-1
                pop!(ids)
            end
        end
    end
    empty!(ids[end])
    return ids[end]
end

function Base.fill!(target_ids::T, source_ids::T) where {T<:IDS}
    for field in fieldnames(typeof(target_ids))
        if field == :_parent
            continue
        end
        setproperty!(target_ids, field, deepcopy(getfield(source_ids, field)))
    end
    return target_ids
end

function _set_conditions(ids::IDS, conditions::Pair{String}...)
    for (path, value) in conditions
        h = ids
        for p in i2p(path)
            if typeof(p) <: Int
                if p > length(h)
                    resize!(h, p)
                end
                h = h[p]
            else
                p = Symbol(p)
                if ismissing(h, p)
                    setproperty!(h, p, value)
                end
                h = getproperty(h, p)
            end
        end
    end
    return ids
end

function _match(ids::IDSvector{T}, conditions) where {T<:IDSvectorElement}
    matches = Dict()
    for (k, item) in enumerate(ids)
        match = true
        for (path, value) in conditions
            h = item
            for p in i2p(path)
                if typeof(p) <: Int
                    if p > length(h)
                        match = false
                        break
                    end
                    h = h[p]
                else
                    p = Symbol(p)
                    if ismissing(h, p)
                        match = false
                        break
                    end
                    h = getproperty(h, p)
                end
            end
            if h != value
                match = false
                break
            end
        end
        if match
            matches[k] = item
        end
    end
    return matches
end

"""
    Base.resize!(ids::IDSvector{T}, conditions...) where {T <: IDSvectorElement}

Resize if a set of conditions are not met
If an entry matching the condition is found, then the content of the matching IDS is emptied, and the IDS is populated with the conditions
Returns selected IDS
"""
function Base.resize!(ids::IDSvector{T}, condition::Pair{String}, conditions::Pair{String}...; allow_multiple_matches = false) where {T<:IDSvectorElement}
    conditions = vcat(condition, collect(conditions))
    if length(ids) == 0
        return _set_conditions(resize!(ids, 1), conditions...)
    end
    matches = _match(ids, conditions)
    if length(matches) == 1
        return _set_conditions(empty!(collect(values(matches))[1]), conditions...)
    elseif length(matches) > 1
        if allow_multiple_matches
            for (kk, k) in reverse(collect(enumerate(keys(matches))))
                if kk == 1
                    return _set_conditions(empty!(matches[k]), conditions...)
                else
                    deleteat!(ids, k)
                end
            end
        else
            error("Multiple entries $([k for k in keys(matches)]) match resize! conditions: $conditions")
        end
    else
        return _set_conditions(resize!(ids, length(ids) + 1), conditions...)
    end
end

"""
    Base.deleteat!(ids::IDSvector{T}, conditions...) where {T <: IDSvectorElement}

If an entry matching the condition is found, then the content of the matching IDS is emptied
"""
function Base.deleteat!(ids::IDSvector{T}, condition::Pair{String}, conditions::Pair{String}...) where {T<:IDSvectorElement}
    conditions = vcat(condition, collect(conditions))
    if length(ids) == 0
        return ids
    end
    matches = _match(ids, conditions)
    for k in reverse(collect(keys(matches)))
        deleteat!(ids, k)
    end
    return ids
end
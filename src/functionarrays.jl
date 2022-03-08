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
    info(location::String)

Return information of a node in the IMAS data structure
"""
function info(location::String)
    if location == "dd.global_time"
        return Dict("units" => "s", "documentation" => "Generic global time")
    end
    location = replace(location, r"\[[0-9]+\]$" => "[:]")
    location = replace(location, r"\[:\]$" => "")
    return imas_dd_ids(split(location, ".")[1])[location]
end

"""
    info(ids::IDS, field::Symbol)

Return information of a filed of an IDS
"""
function info(ids::Union{IDS,IDSvector}, field::Symbol)
    location = "$(f2u(ids)).$(field)"
    return info(location)
end

function info(ids::Type, field::Symbol)
    if ids <: IDSvectorElement
        location = "$(_f2u(ids))[:].$(field)"
    else
        location = "$(_f2u(ids)).$(field)"
    end
    return info(location)
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

"""
    units(ids::IDS, field::Symbol)

Return string with units
"""
function units(location::String)
    return get(info(location), "units", "")
end

"""
    units(ids::IDS, field::Symbol)

Return string with units
"""
function units(ids::IDS, field::Symbol)
    location = "$(f2u(ids)).$(field)"
    return units(location)
end

"""
    coordinates(ids::IDS, field::Symbol)

Return two lists, one of coordinate names and the other with their values in the data structure
Coordinate value is `nothing` when the data does not have a coordinate
Coordinate value is `missing` if the coordinate is missing in the data structure
"""
function coordinates(ids::IDS, field::Symbol)
    txt = info(ids, field)
    # handle scalar quantities (which do not have coordinate)
    if !("coordinates" in keys(txt))
        return Dict(:names => [], :values => [])
    end
    coord_names = deepcopy(txt["coordinates"])
    coord_values = []
    for coord in coord_names
        if occursin("...", coord)
            push!(coord_values, nothing)
        else
            try
                h = goto(ids, u2fs(p2i(i2p(coord)[1:end-1])); f2 = f2fs)
                coord_leaf = Symbol(i2p(coord)[end])
                h = getfield(h, coord_leaf)
                # add value to the coord_values
                push!(coord_values, h)
            catch e
                if typeof(e) <: IMASdetachedHead
                    push!(coord_values, missing)
                else
                    rethrow()
                end
            end
        end
    end
    return Dict(:names => coord_names, :values => coord_values)
end

#= ========== =#
#  Exceptions  #
#= ========== =#

struct IMASdetachedHead <: Exception
    source::String
    destination::String
end

Base.showerror(io::IO, e::IMASdetachedHead) = print(io, "Could not reach `$(e.destination)` from `$(e.source)`")

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

struct AccessLog
    read::Set{String}
    expr::Set{String}
    write::Set{String}
end

const access_log = AccessLog(Set(String[]), Set(String[]), Set(String[]))

function empty_access_log!()
    empty!(access_log.read)
    empty!(access_log.expr)
    empty!(access_log.write)
end

function Base.show(io::IO, ::MIME"text/plain", access_log::AccessLog)
    for field in [:read, :expr, :write]
        log = getfield(access_log, field)
        for k in sort(collect(log))
            println(io, "$field: $k")
        end
    end
end

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
            value = exec_expression_with_ancestor_args(ids, field, value, x)
            push!(access_log.expr, "$(f2u(ids)).$(field)")
            return value
        catch e
            if typeof(e) <: IMASexpressionRecursion
                rethrow()
            else
                throw(IMASbadExpression(ids, field, sprint(showerror, e)))
            end
        end
    elseif field == :_parent
        return value
    else
        push!(access_log.read, "$(f2u(ids)).$(field)")
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
        tmp = setfield!(ids, field, v)
        if !(typeof(v) <: Union{IDS,IDSvector}) && (field != :_parent)
            push!(access_log.write, "$(f2u(ids)).$(field)")
        end
        return tmp
    catch e
        if typeof(e) <: TypeError
            error("$(typeof(v)) is the wrong type for $(f2u(ids)).$(field), it should be $(fieldtype(typeof(ids), field)))")
        else
            rethrow()
        end
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
    #@info ">"^length(expression_call_stack) * " " * structure_name
    # keep track of recursion
    if !(structure_name in expression_call_stack)
        push!(expression_call_stack, structure_name)
    else
        throw(IMASexpressionRecursion(copy(expression_call_stack)))
    end
    try
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

function Base.insert!(ids::IDSvector{T}, i, v::T) where {T<:IDSvectorElement}
    setfield!(v, :_parent, WeakRef(ids))
    insert!(ids._value, i, v)
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

#= ===== =#
#  UTILS  #
#= ===== =#
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

#= ==== =#
#  keys  #
#= ==== =#
"""
    _common_base_string(s1::String, s2::String)::Vector{String}

given two strings it returns a tuple of 3 strings that is the common initial part, and then the remaining parts
"""
function _common_base_string(s1::String, s2::String)::Tuple{String,String,String}
    index = nothing
    for k = 1:min(length(s1), length(s2))
        sub = SubString(s2, 1, k)
        if startswith(s1, sub)
            index = k
        end
    end
    if index === nothing
        return "", s1, s2
    else
        return string(SubString(s1, 1, index)), string(SubString(s1, index + 1, length(s1))), string(SubString(s2, index + 1, length(s2)))
    end
end

"""
    Base.keys(ids::IDS)

Returns list of fields with data in a IDS
"""
function Base.keys(ids::IDS)
    kkk = Symbol[]
    for k in fieldnames(typeof(ids))
        # hide the _parent field
        if k === :_parent
            continue
        end
        v = getfield(ids, k)
        # empty entries
        if v === missing
            continue
            # empty structures/arrays of structures (recursive)
        elseif typeof(v) <: Union{IDS,IDSvector}
            if length(keys(v)) > 0
                push!(kkk, k)
            end
            # entries with data
        else
            push!(kkk, k)
        end
    end
    return kkk
end

function Base.keys(ids::IDSvector)
    return collect(1:length(ids))
end


#= ====== =#
#  empty!  #
#= ====== =#

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

#= ======= =#
#  resize!  #
#= ======= =#

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

#= ========= =#
#  deleteat!  #
#= ========= =#

function Base.deleteat!(ids::IDSvector{T}, i::Int) where {T<:IDSvectorElement}
    return deleteat!(ids._value, i)
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

#= ========= =#
#  ismissing  #
#= ========= =#

"""
    ismissing(ids::IDS, leaf)::Bool

returns true/false if field is missing in IDS
"""
function Base.ismissing(ids::IDS, field::Symbol)::Bool
    try
        getproperty(ids, field)
        return false
    catch
        return true
    end
end

function Base.ismissing(ids::IDS, field::Vector)::Bool
    if length(field) == 1
        return ismissing(ids, field[1])
    end
    return ismissing(getfield(ids, field[1]), field[2:end])
end

function Base.ismissing(ids::IDSvector, field::Int)::Bool
    return length(ids) < field
end

function Base.ismissing(ids::IDSvector, field::Vector)::Bool
    if length(field) == 1
        return ismissing(ids, field[1])
    end
    if field[1] <= length(ids)
        return ismissing(ids[field[1]], field[2:end])
    else
        return true
    end
end

#= ==== =#
#  diff  #
#= ==== =#

_diff_function(v1, v2, tol) = maximum(abs.(v1 .- v2) ./ (tol + sum(abs.(v1) .+ abs.(v2)) / 2.0 / length(v1)))

"""
    diff(ids1::IDS, ids2::IDS, path = Any[]; tol = 1E-2, plot_function = nothing, recursive = true, differences = Dict(), verbose=true)

Compares two IDSs and returns dictionary with differences
"""
function Base.diff(ids1::IDS, ids2::IDS, path = Any[]; tol = 1E-2, plot_function = nothing, recursive = true, differences = Dict(), verbose = true)
    @assert typeof(ids1) === typeof(ids2) "IDS types are different:  $(typeof(ids1))  and  $(typeof(ids2))"
    for field in sort(collect(fieldnames(typeof(ids2))))
        if field === :_parent
            continue
        end
        v1 = missing
        v2 = missing
        try
            v1 = getproperty(ids1, field)
        catch e
            if typeof(e) <: Union{IMASmissingDataException,IMASexpressionError}
                v1 = missing
            else
                rethrow()
            end
        end
        try
            v2 = getproperty(ids2, field)
        catch e
            if typeof(e) <: Union{IMASmissingDataException,IMASexpressionError}
                v2 = missing
            else
                rethrow()
            end
        end

        pathname = p2i(vcat(path, field))
        if v1 === v2 === missing
            continue
        elseif typeof(v1) != typeof(v2)
            differences[pathname] = "types  $(typeof(v1)) --  $(typeof(v2))"
        elseif typeof(v1) <: IDS
            if recursive
                diff(v1, v2, vcat(path, field); tol, plot_function, recursive, differences)
            end
        elseif typeof(v1) <: IDSvector
            if recursive
                if length(v1) != length(v2)
                    differences[pathname] = "length  $(length(v1)) --  $(length(v2))"
                else
                    for k = 1:length(v1)
                        diff(v1[k], v2[k], vcat(path, field, k); tol, plot_function, recursive, differences)
                    end
                end
            end
        elseif typeof(v1) <: AbstractArray
            if length(v1) != length(v2)
                differences[pathname] = "length  $(length(v1)) --  $(length(v2))"
            elseif _diff_function(v1, v2, tol) > tol
                if (typeof(v1) <: AbstractVector) && (plot_function !== nothing)
                    display(plot_function(1:length(v1), [v1, v2], title = pathname, label = ""))
                end
                differences[pathname] = @sprintf("value %g", _diff_function(v1, v2, tol))
            end
        elseif typeof(v1) <: Number
            if _diff_function(v1, v2, tol) > tol
                differences[pathname] = "value  $v1 --  $v2"
            end
        elseif typeof(v1) <: Union{String,Symbol}
            if v1 != v2
                differences[pathname] = "value  $v1 --  $v2"
            end
        else
            println((pathname, typeof(v1), typeof(v2)))
            asdadas # raise error to force use to handle this explicitly
            differences[pathname] = "value  $v1 --  $v2"
        end
        if verbose && (pathname in keys(differences))
            printstyled(pathname; bold = true)
            printstyled(" ➡ "; color = :red)
            println(differences[pathname])
        end
    end
    return differences
end

#= ========== =#
#  navigation  #
#= ========== =#

"""
    top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)

Return top-level IDS in the DD hierarchy.
Considers IDS as maximum top level if IDS_is_absolute_top=true
"""
function top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool = true)
    if IDS_is_absolute_top & (typeof(ids) <: dd)
        error("Cannot call top(x::dd,IDS_is_absolute_top=true). Use `IDS_is_absolute_top=false`.")
    elseif !(typeof(ids._parent.value) <: Union{IDS,IDSvector})
        return ids
    elseif IDS_is_absolute_top & (typeof(ids._parent.value) <: dd)
        return ids
    else
        return top(ids._parent.value; IDS_is_absolute_top = IDS_is_absolute_top)
    end
end

"""
    top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)

Return top-level IDS in the DD hierarchy.
Considers IDS as maximum top level if IDS_is_absolute_top=true
"""

function top_ids(ids::Union{IDS,IDSvector})
    ids = top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top = true)
    if length(f2p(ids)) == 1
        return ids
    else
        return missing
    end
end

function top_dd(ids::Union{IDS,IDSvector})
    ids = top(ids; IDS_is_absolute_top = false)
    if typeof(ids) <: dd
        return ids
    else
        return missing
    end
end

"""
    parent(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)

Return parent IDS/IDSvector in the hierarchy
If IDS_is_absolute_top then returns `missing` instead of dd()
"""
function parent(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool = true)
    if !(typeof(ids._parent.value) <: Union{IDS,IDSvector})
        return missing
    elseif IDS_is_absolute_top & (typeof(ids._parent.value) <: dd)
        return missing
    else
        return ids._parent.value
    end
end

"""
    children(ids::Union{IDS,IDSvector})::Vector{Symbol}

Return children of a IDS/IDSvector
"""
function children(ids::Union{IDS,IDSvector})::Vector{Symbol}
    return [k for k in fieldnames(typeof(ids)) if k != :_parent]
end

"""
    goto(ids::IDS, location::String)

Reach location in a given IDS

# Arguments
- `f2::Function=f2i`: function used to process the IDS path to be compared to `location`
"""
function goto(ids::Union{IDS,IDSvector}, location::String; f2::Function = f2i)
    # find common ancestor
    cs, s1, s2 = _common_base_string(f2(ids), location)
    cs0 = replace(cs, r"\.$" => "")
    # go upstream until common acestor
    h = ids
    while f2(h) != cs0
        if h._parent.value === missing
            break
        end
        h = h._parent.value
    end
    # then dive into the location branch
    for k in i2p(s2)
        try
            h = getfield(h, Symbol(k))
        catch
            throw(IMASdetachedHead("$(f2(ids))", "$(location)"))
        end
    end
    return h
end

#= ========= =#
#  save/load  #
#= ========= =#

"""
    dict2imas(dct, ids::IDS=dd() ;verbose::Bool=false, path::Vector{String}=String[])::IDS

Populate IMAS data structure `ids` based on data contained in Julia dictionary `dct`.

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
- `skip_non_coordinates::Bool=false`: only assign coordinates to the data structure
"""
function dict2imas(dct, ids::T; verbose::Bool = false, path::Vector{String} = String[], skip_non_coordinates::Bool = false) where {T<:IDS}
    # recursively traverse `dtc` structure
    level = length(path)
    for (k, v) in dct
        if !hasfield(typeof(ids), Symbol(k))
            if !skip_non_coordinates
                @warn("$(f2i(ids)).$(k) was skipped in IMAS.jl data dictionary", maxlog = 1)
            end
            continue
        end
        # Struct
        if typeof(v) <: Dict
            if verbose
                println(("｜"^level) * string(k))
            end
            ff = getfield(ids, Symbol(k))
            dict2imas(v, ff; path = vcat(path, [string(k)]), verbose = verbose, skip_non_coordinates = skip_non_coordinates)

            # Array of struct
        elseif (typeof(v) <: Array) && (length(v) > 0) && (typeof(v[1]) <: Dict)
            ff = getfield(ids, Symbol(k))
            if verbose
                println(("｜"^level) * string(k))
            end
            if length(ff) < length(v)
                resize!(ff, length(v))
            end
            for i = 1:length(v)
                if verbose
                    println(("｜"^(level + 1)) * string(i))
                end
                dict2imas(v[i], ff[i]; path = vcat(path, [string(k), "[$i]"]), verbose = verbose, skip_non_coordinates = skip_non_coordinates)
            end

            # Leaf
        else
            if verbose
                print(("｜"^level) * string(k) * " → ")
            end
            target_type = Core.Compiler.typesubtract(struct_field_type(typeof(ids), Symbol(k)), Union{Missing,Function}, 1)
            if target_type <: AbstractArray
                if tp_ndims(target_type) == 2
                    v = transpose(reduce(hcat, v))
                end
                if tp_eltype(target_type) <: Real
                    v = convert(Array{Float64,tp_ndims(target_type)}, v)
                else
                    v = convert(Array{tp_eltype(target_type),tp_ndims(target_type)}, v)
                end
            end
            setproperty!(ids, Symbol(k), v; skip_non_coordinates = skip_non_coordinates)
            if verbose
                println(typeof(v))
            end
        end
    end

    return ids
end

tp_ndims(::Type{AbstractArray{T,N} where {T,N}}) = N
tp_ndims(v::UnionAll) = ndims(v.body)
tp_eltype(::Type{AbstractArray{T,N} where {T,N}}) = T
tp_eltype(v::UnionAll) = v.var.ub

"""
    json2imas(filename::String; verbose::Bool=false)::IDS

Load from a file with give `filename` the IMAS data structure saved in JSON format 

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
"""
function json2imas(filename::String; verbose::Bool = false)::IDS
    ids_data = dd()
    json_data = JSON.parsefile(filename)
    dict2imas(json_data, ids_data; verbose = verbose, skip_non_coordinates = true)
    dict2imas(json_data, ids_data; verbose = verbose, skip_non_coordinates = false)
    return ids_data
end

"""
    imas2dict(ids::Union{IDS,IDSvector})

Populate Julia structure of dictionaries and vectors with data from IMAS data structure `ids`
"""
function imas2dict(ids::Union{IDS,IDSvector})
    if typeof(ids) <: IDSvector
        dct = Any[]
    else
        dct = Dict()
    end
    return imas2dict(ids, dct)
end

function imas2dict(ids::Union{IDS,IDSvector}, dct::Union{Dict,Vector})
    items = sort(keys(ids))
    for item in items
        if typeof(ids) <: IDSvector # arrays of structures
            push!(dct, Dict())
            imas2dict(ids[item], dct[item])
        else
            value = missing
            try
                value = getproperty(ids, item)
            catch e
                if typeof(e) <: IMASexpressionError
                    #@debug(sprint(showerror, e))
                    continue
                else
                    throw
                end
            end
            if typeof(value) <: Union{IDS,IDSvector} # structures
                if typeof(value) <: IDS
                    dct[item] = Dict()
                else
                    dct[item] = Any[]
                end
                imas2dict(value, dct[item])
            else # field
                dct[item] = value
            end
        end
    end
    return dct
end

"""
    imas2json(ids::Union{IDS,IDSvector}, filename::String; kw...)

Save the IMAS data structure to a JSON file with give `filename`
"""
function imas2json(ids::Union{IDS,IDSvector}, filename::String; kw...)
    open(filename, "w") do io
        JSON.print(io, imas2dict(ids); kw...)
    end
end

#= ==== =#
#  show  #
#= ==== =#

function Base.show(io::IO, ids::Union{IDS,IDSvector}, depth::Int)
    items = keys(ids)
    for (k, item) in enumerate(sort(items))
        printstyled(io, "$('｜'^depth)"; color = :yellow)
        # arrays of structurs
        if typeof(ids) <: IDSvector
            printstyled(io, "[$(item)]\n"; bold = true, color = :green)
            show(io, ids[item], depth + 1)
        else
            value = getfield(ids, item)
            if typeof(value) <: Union{IDS,IDSvector}
                # structures
                if (typeof(ids) <: dd)
                    printstyled(io, "$(uppercase(string(item)))\n"; bold = true)
                else
                    printstyled(io, "$(string(item))\n"; bold = true)
                end
                show(io, value, depth + 1)
            else
                # field
                printstyled(io, "$(item)")
                printstyled(io, " ➡ "; color = :red)
                if typeof(value) <: Function
                    printstyled(io, "Function\n"; color = :blue)
                elseif typeof(value) <: String
                    printstyled(io, "\"$(value)\"\n"; color = :magenta)
                elseif typeof(value) <: Integer
                    printstyled(io, "$(value)\n"; color = :yellow)
                elseif typeof(value) <: AbstractFloat
                    printstyled(io, @sprintf("%g\n", value); color = :red)
                elseif typeof(value) <: AbstractArray
                    if length(value) < 5
                        if eltype(value) <: AbstractFloat
                            printstyled(io, "[$(join([@sprintf("%g",v) for v in value],","))]\n"; color = :green)
                        else
                            printstyled(io, "$value\n"; color = :green)
                        end
                    else
                        printstyled(io, "$(Base.summary(value))\n"; color = :green)
                    end
                else
                    printstyled(io, "$(Base.summary(value))\n"; color = :blue)
                end
            end
        end
        if (typeof(ids) <: dd) & (k < length(items))
            println(io, "")
        end
    end
end

# show function for the Jupyter notebook
function Base.show(io::IO, ::MIME"text/plain", ids::Union{IDS,IDSvector})
    return show(io, ids, 0)
end

# show function for inline prints
function Base.show(io::IO, ids::IDS)
    fnames = []
    for item in keys(ids)
        push!(fnames, item)
    end
    return println(io, "$(f2i(ids)){$(join(collect(map(x->"$x",fnames)),", "))}")
end

# show function for inline prints
function Base.show(io::IO, ids::IDSvector)
    fnames = []
    for item in keys(ids)
        push!(fnames, item)
    end
    if length(ids) < 2
        return println(io, "$(p2i(f2p(ids)[1:end-1]))[$(join(collect(map(x->"$x",fnames)),", "))]")
    else
        return println(io, "$(p2i(f2p(ids)[1:end-1]))[1...$(length(ids))]")
    end
end

#= ==== =#
#  time  #
#= ==== =#
include("time.jl")
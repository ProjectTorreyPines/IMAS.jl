import JSON

include("dd.jl")

include("utils.jl")

import Base:resize!

"""
    resize!(collection::IDSvector{T}, n::Int) where {T<:IDS}

Change size of a IDS array of structures
"""
function resize!(collection::IDSvector{T}, n::Int) where {T <: IDS}
    if n > length(collection)
        for k in length(collection):n - 1
            obj = eltype(collection)()
            setfield!(obj, :_parent, WeakRef(collection))
            push!(collection.value, obj)
            # println("push $(length(collection))")
        end
    elseif n < length(collection)
        for k in n:length(collection) - 1
            pop!(collection.value)
            # println("pop $(length(collection))")
        end
    end
    return collection
end

"""
    dict2imas(dct, ids::IDS=dd() ;verbose::Bool=false, path::Vector{String}=String[])::IDS

Populate IMAS data structure `ids` based on data contained in Julia dictionary `dct`.

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
- `skip_non_coordinates::Bool=false`: only assign coordinates to the data structure
"""
function dict2imas(dct, ids::T ;verbose::Bool=false, path::Vector{String}=String[], skip_non_coordinates::Bool=false) where {T <: IDS}
    # recursively traverse `dtc` structure
    level = length(path)
    for (k, v) in dct
        # Struct
        if typeof(v) <: Dict
            if verbose println(("｜"^level) * string(k)) end
            ff = getfield(ids, Symbol(k))
            dict2imas(v, ff; path=vcat(path, [string(k)]), verbose=verbose, skip_non_coordinates=skip_non_coordinates)

        # Array of struct
        elseif (typeof(v) <: Array) && (length(v) > 0) && (typeof(v[1]) <: Dict)
            ff = getfield(ids, Symbol(k))
            if verbose println(("｜"^level) * string(k)) end
            resize!(ff, length(v))
            for i in 1:length(v)
                if verbose println(("｜"^(level + 1)) * string(i)) end
                dict2imas(v[i], ff[i]; path=vcat(path, [string(k),"[$i]"]), verbose=verbose, skip_non_coordinates=skip_non_coordinates)
            end

        # Leaf
        else
            if verbose print(("｜"^level) * string(k) * " → ") end
            target_type = Core.Compiler.typesubtract(struct_field_type(typeof(ids), Symbol(k)), Union{Missing,Function}, 1)
            if target_type <: AbstractArray
                if ndims(target_type) == 2
                    v = reduce(hcat, v)
                end
                v = convert(Array{eltype(target_type),ndims(target_type)}, v)
            end
            setproperty!(ids, Symbol(k), v; skip_non_coordinates=skip_non_coordinates)
            if verbose println(typeof(v)) end
        end
    end

    return ids
end

Base.ndims(::Type{Vector{T} where T <: Real}) = 1
Base.ndims(::Type{Matrix{T} where T <: Real}) = 2
Base.eltype(::Type{Vector{T} where T <: Real}) = Real
Base.eltype(::Type{Matrix{T} where T <: Real}) = Real

"""
    json2imas(filename::String; verbose::Bool=false)::IDS

Load from a file with give `filename` the OMAS data structure saved in JSON format 

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
"""
function json2imas(filename::String; verbose::Bool=false)::IDS
    ids_data = dd()
    json_data = JSON.parsefile(filename)
    dict2imas(json_data, ids_data; verbose=verbose, skip_non_coordinates=true)
    dict2imas(json_data, ids_data; verbose=verbose, skip_non_coordinates=false)
    return ids_data
end

"""
    top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)

Return top-level IDS in the DD hierarchy.
Considers IDS as maximum top level if IDS_is_absolute_top=true
"""
function top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)
    if IDS_is_absolute_top & (typeof(ids) <: dd)
        error("Cannot call top(x::IMAS.dd,IDS_is_absolute_top=true). Use `IDS_is_absolute_top=false`.")
    elseif ids._parent.value === missing
        return ids
    elseif IDS_is_absolute_top & (typeof(ids._parent.value) <: dd)
        return ids
    else
        return top(ids._parent.value;IDS_is_absolute_top=IDS_is_absolute_top)
    end
end

"""
    parent(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)

Return parent IDS/IDSvector in the hierarchy
If IDS_is_absolute_top then returns `missing` instead of IMAS.dd()
"""
function parent(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)
    if ids._parent.value === missing
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
    assign_expressions(ids::Union{IDS,IDSvector})

Assign expressions to a IDS/IDSvector
NOTE: This is done not recursively
"""
function assign_expressions(ids::Union{IDS,IDSvector})
    struct_name = f2u(ids)
    for item in children(ids)
        if typeof(getfield(ids, item)) <: Union{IDS,IDSvector}
            continue
        elseif "$(struct_name).$(item)" in keys(expressions)
            setproperty!(ids, item, expressions["$(struct_name).$(item)"])
        end
    end
end

#= ===================== =#
#  IDS related functions  #
#= ===================== =#

"""
    coordinates(ids::IDS, field::Symbol)

Returns two lists, one of coordinate names and the other with their values in the data structure
Coordinate value is `nothing` when the data does not have a coordinate
Coordinate value is `missing` if the coordinate is missing in the data structure
"""
function coordinates(ids::IDS, field::Symbol)
    info = imas_info("$(f2u(ids)).$(field)")
    # handle scalar quantities (which do not have coordinate)
    if ! ("coordinates" in keys(info))
        return Dict(:names => [], :values => [])
    end
    coord_names = deepcopy(info["coordinates"])
    coord_values = []
    for (k, coord) in enumerate(coord_names)
        if occursin("...", coord)
            push!(coord_values, nothing)
        else
            # find common ancestor
            s1 = f2fs(ids)
            s2 = u2fs(p2i(i2p(coord)[1:end - 1]))
            cs, s1, s2 = IMAS.common_base_string(s1, s2)
            # go upstream until common acestor
            h = ids
            while f2fs(h) != cs
                h = h._parent.value
            end
            # then dive into the coordinates branch
            for k in i2p(s2)
                h = getfield(h, Symbol(k))
            end
            coord_leaf = i2p(coord)[end]
            h = getfield(h, Symbol(coord_leaf))
            # add value to the coord_values
            push!(coord_values, h)
        end
    end
    return Dict(:names => coord_names, :values => coord_values)
end

"""
    f2u(ids::Union{IDS, IDSvector, DataType, Symbol, String})

Returns universal IMAS location of a given IDS
"""
function f2u(ids::IDS)
    return _f2u(typeof(ids))
end

function f2u(ids::IDSvector)
    return _f2u(eltype(ids)) * "[:]"
end

function f2u(ids::IDSvectorElement)
    return _f2u(typeof(ids)) * "[:]"
end

function _f2u(ids::DataType)
    return _f2u(Base.typename(ids).name)
end

function _f2u(ids::Symbol)
    return _f2u(string(ids))
end

function _f2u(ids::String)
    if in(':', ids) | in('.', ids)
        error("`$ids` is not a qualified IDS type")
    end
    tmp = replace(ids, "___" => "[:].")
    tmp = replace(tmp, "__" => ".")
    return replace(tmp, r"\.$" => "")
end

"""
    return IDS type as a string
"""
function f2fs(ids::IDS)::String
    return _f2fs(typeof(ids))
end

function f2fs(ids::IDSvector)::String
    return _f2fs(eltype(ids))
end

function f2fs(ids::IDSvectorElement)::String
    return _f2fs(typeof(ids)) * "___"
end

function _f2fs(ids::DataType)::String
    return string(Base.typename(ids).name)
end


"""
    f2p(ids::Union{IDS,IDSvector})::Vector{Union{String,Int}}

Returns IMAS location of a given IDS

NOTE: indexes of arrays of structures that cannot be determined are set to 0
"""
function f2p(ids::Union{IDS,IDSvector})::Vector{Union{String,Int}}
    return f2p(ids, missing, nothing, Int[])
end

function f2p(ids::Union{IDS,IDSvector},
             child::Union{Missing,IDS,IDSvector},
             path::Union{Nothing,Vector},
             index::Vector{Int})
    # initialize path and index
    if path === nothing
        if typeof(ids) <: IDS
            if typeof(ids._parent.value) <: IDSvector
                name = string(Base.typename(typeof(ids)).name) * "___"
            else
                name = string(Base.typename(typeof(ids)).name)
            end
        elseif typeof(ids) <: IDSvector
            name = string(Base.typename(eltype(ids)).name) * "___"
        end
        path = replace(name, "___" => "__:__")
        path = Vector{Any}(Vector{String}(split(path, "__")))
        path = [k == ":" ? 0 : k for k in path if length(k) > 0]
    end

    # collect integers for arrays of structures
    if typeof(ids) <: IDSvector
        ix = findfirst([k === child for k in ids.value])
        if ix === nothing
            push!(index, 0)
        else
            push!(index, ix)
        end
    end

    # traverse IDSs upstream or return result once top is reached
    if ids._parent.value === missing
        index = reverse(index)
        path = reverse([(typeof(k) <: Int) & (length(index) > 0) ? pop!(index) : k for k in reverse(path)])
        return path
    else
        return f2p(ids._parent.value, ids, path, index)
    end
end

"""
    i2p(imas_location::String)::Vector{Union{String,Int}}

Split IMAS location in its elements
"""
function i2p(imas_location::String)::Vector{Union{String,Int}}
    path = []
    if length(imas_location) == 0
        return path
    end
    for s in split(imas_location, '.')
        if in('[', s)
            s, n = split(s, '[')
            n = strip(n, ']')
            s = string(s)
            push!(path, s)
            if n == ":"
                push!(path, ":")
            else
                push!(path, parse(Int, n))
        end
        else
            s = string(s)
            push!(path, s)
        end
    end
return path
end


"""
    p2i(path::Vector{Union{String,Int}})::String

Combine list of IMAS location elements into a string
"""
function p2i(path::Vector{Union{String,Int}})::String
    str = String[]
    for item in path
        if item == ":"
            push!(str, "[:]")
        elseif typeof(item) <: Int
            push!(str, "[$(string(item))]")
        else
            push!(str, item)
        end
    end
    return replace(join(str, "."), ".[" => "[")
end

function p2i(path::Vector)::String
    return p2i(Vector{Union{String,Int}}(path))
end


"""
    u2fs(imas_location::String)::String

return IDS/IDSvector type as a string starting from a universal IMAS location string
"""
function u2fs(imas_location::String)::String
    tmp = replace(imas_location, "." => "__")
    return replace(tmp, "[:]" => "_")
end

"""
    return IDS/IDSvector type starting from a universal IMAS location string

u2f(imas_location::String)
"""
function u2f(imas_location::String)
    return eval(Meta.parse(u2fs(imas_location)))
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

function Base.show(io::IO, ids::Union{IDS,IDSvector}, depth::Int)
    items = keys(ids)
    for (k, item) in enumerate(items)
        printstyled("$('｜'^depth)"; color=:yellow)
        # arrays of structurs
        if typeof(ids) <: IDSvector
            printstyled("[$(item)]\n"; bold=true, color=:green)
            show(io, ids[item], depth + 1)
        # structures
        elseif typeof(getfield(ids, item)) <: Union{IDS,IDSvector}
            if (typeof(ids) <: dd)
                printstyled("$(uppercase(string(item)))\n"; bold=true)
            else
                printstyled("$(string(item))\n"; bold=true)
            end
            show(io, getfield(ids, item), depth + 1)
        # field
        else
            value = getfield(ids, item)
            printstyled("$(item)")
            color = :red
            printstyled(" ➡ "; color=color)
            printstyled("$(Base.summary(getfield(ids, item)))\n"; color=:blue)
        end
        if (typeof(ids) <: dd) & (k < length(items))
            println()
end
    end
end

function Base.show(io::IO, ids::Union{IDS,IDSvector})
    return show(io, ids, 0)
end
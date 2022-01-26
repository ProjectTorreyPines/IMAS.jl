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
        if ! hasfield(typeof(ids), Symbol(k))
            if ! skip_non_coordinates
                println("$(f2i(ids)).$(k) was skipped in IMAS.jl data dictionary")
            end
            continue
        end
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
                if tp_ndims(target_type) == 2
                    v = transpose(reduce(hcat, v))
                end
                if tp_eltype(target_type) <: Real
                    v = convert(Array{Float64,tp_ndims(target_type)}, v)
                else
                    v = convert(Array{tp_eltype(target_type),tp_ndims(target_type)}, v)
                end
            end
            setproperty!(ids, Symbol(k), v; skip_non_coordinates=skip_non_coordinates)
            if verbose println(typeof(v)) end
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
function json2imas(filename::String; verbose::Bool=false)::IDS
    ids_data = dd()
    json_data = JSON.parsefile(filename)
    dict2imas(json_data, ids_data; verbose=verbose, skip_non_coordinates=true)
    dict2imas(json_data, ids_data; verbose=verbose, skip_non_coordinates=false)
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
        if typeof(ids) <: IDSvector
            # arrays of structures
            push!(dct, Dict())
            imas2dict(ids[item],dct[item])
        else
            value = getproperty(ids, item)
            # structures
            if typeof(value) <: Union{IDS,IDSvector}
                if typeof(value) <: IDS
                    dct[item] = Dict()
                else
                    dct[item] = Any[]
                end
                imas2dict(value, dct[item])
            # field
            else
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
    goto(ids::IDS, location::String)

Reach location in a given IDS

# Arguments
- `f2::Function=f2i`: function used to process the IDS path to be compared to `location`
"""
function goto(ids::IDS, location::String; f2::Function=f2i)
    # find common ancestor
    cs, s1, s2 = IMAS.common_base_string(f2(ids), location)
    cs0 = replace(cs, r"\.$" => "")
    # go upstream until common acestor
    h = ids
    while f2(h) != cs0
        h = h._parent.value
        if h === missing
            error("Could not reach `$(location)` from `$(f2(ids))`")
        end
    end
    # then dive into the location branch
    for k in i2p(s2)
        h = getfield(h, Symbol(k))
    end
    return h
end

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
            h = goto(ids, u2fs(p2i(i2p(coord)[1:end - 1])); f2=f2fs)
            coord_leaf = Symbol(i2p(coord)[end])
            h = getfield(h, coord_leaf)
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
    f2fs(ids::IDS)::String

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
    f2i(ids::Union{IDS,IDSvector})::String

return IMAS location of a given IDS
"""
function f2i(ids::Union{IDS,IDSvector})::String
    return p2i(f2p(ids))
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
    i2u(imas_location::String)::String

return universal IMAS location from IMAS location
ie. replaces indexes of arrays of structures with [:]
"""
function i2u(imas_location::String)::String
    return replace(imas_location, r"\[[0-9]+\]" => "[:]")
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
    for (k, item) in enumerate(sort(items))
        printstyled(io, "$('｜'^depth)"; color=:yellow)
        # arrays of structurs
        if typeof(ids) <: IDSvector
            printstyled(io, "[$(item)]\n"; bold=true, color=:green)
            show(io, ids[item], depth + 1)
        else
            value = getfield(ids, item)
            # structures
            if typeof(value) <: Union{IDS,IDSvector}
                if (typeof(ids) <: dd)
                    printstyled(io, "$(uppercase(string(item)))\n"; bold=true)
                else
                    printstyled(io, "$(string(item))\n"; bold=true)
                end
                show(io, value, depth + 1)
            # field
            else
                printstyled(io, "$(item)")
                printstyled(io, " ➡ "; color=:red)
                if typeof(value) <: Function
                    printstyled(io, "Function\n"; color=:green)
                elseif typeof(value) <: String
                    printstyled(io, "\"$(value)\"\n"; color=:magenta)
                elseif typeof(value) <: Number
                    printstyled(io, "$(value)\n"; color=:magenta)
                else
                    printstyled(io, "$(Base.summary(value))\n"; color=:blue)
                end
            end
        end
        if (typeof(ids) <: dd) & (k < length(items))
            println(io, "")
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", ids::Union{IDS,IDSvector})
    return show(io, ids, 0)
end

function Base.show(io::IO, ids::Union{IDS,IDSvector})
    fnames = []
    for item in keys(ids)
        push!(fnames, item)
    end
    return println(io, "$(f2i(ids)){$(join(collect(map(String,fnames)),", "))}")
end
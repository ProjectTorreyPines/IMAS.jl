import JSON

include("dd.jl")

include("utils.jl")

import Base:resize!

"""
    resize!(collection::FDSvector{T}, n::Int) where {T<:FDS}

Change size of a FDS array of structures
"""
function resize!(collection::FDSvector{T}, n::Int) where {T <: FDS}
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
    dict2fuse(dct, fds::FDS=dd() ;verbose::Bool=false, path::Vector{String}=String[])::FDS

Populate FUSE data structure `fds` based on data contained in Julia dictionary `dct`.

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
"""
function dict2fuse(dct, fds::T ;verbose::Bool=false, path::Vector{String}=String[]) where {T <: FDS}
    # recursively traverse `dtc` structure
    level = length(path)
    for (k, v) in dct
        # Struct
        if typeof(v) <: Dict
            if verbose println(("｜"^level) * string(k)) end
            ff = getfield(fds, Symbol(k))
            dict2fuse(v, ff; path=vcat(path, [string(k)]), verbose=verbose)

        # Array of struct
        elseif (typeof(v) <: Array) && (length(v) > 0) && (typeof(v[1]) <: Dict)
            ff = getfield(fds, Symbol(k))
            if verbose println(("｜"^level) * string(k)) end
            resize!(ff, length(v))
            for i in 1:length(v)
                if verbose println(("｜"^(level + 1)) * string(i)) end
                dict2fuse(v[i], ff[i]; path=vcat(path, [string(k),"[$i]"]), verbose=verbose)
            end

        # Leaf
        else
            if verbose print(("｜"^level) * string(k) * " → ") end
            target_type = typeintersect(conversion_types, struct_field_type(typeof(fds), Symbol(k)))
            if target_type <: AbstractArray
                if ndims(target_type) == 2
                    v = reduce(hcat, v)
                end
                v = convert(Array{eltype(target_type),ndims(target_type)}, v)
            end
            setproperty!(fds, Symbol(k), v)
            if verbose println(typeof(v)) end
        end
    end

    return fds
end

"""
    json2fuse(filename::String; verbose::Bool=false)::FDS

Load from a file with give `filename` the OMAS data structure saved in JSON format 

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
"""
function json2fuse(filename::String; verbose::Bool=false)::FDS
    fds_data = dd()
    json_data = JSON.parsefile(filename)
    dict2fuse(json_data, fds_data; verbose=verbose)
    return fds_data
end

"""
    top(fds::Union{FDS, FDSvector})

Return top-level FDS in the DD hierarchy
"""
function top(fds::Union{FDS,FDSvector}; stop_at_ids::Bool=true)
    if stop_at_ids & (typeof(fds) <: dd)
        error("Cannot call top(x::FUSE.dd,stop_at_ids=true). Use `stop_at_ids=false`.")
    elseif fds._parent.value === missing
        return fds
    elseif stop_at_ids & (typeof(fds._parent.value) <: dd)
        return fds
    else
        return top(fds._parent.value;stop_at_ids=stop_at_ids)
    end
end

#= ===================== =#
#  FDS related functions  #
#= ===================== =#

"""
    coordinates(fds::FDS, field::Symbol)

Returns two lists, one of coordinate names and the other with their values in the data structure
Coordinate value is `nothing` when the data does not have a coordinate
Coordinate value is `missing` if the coordinate is missing in the data structure
"""
function coordinates(fds::FDS, field::Symbol)
    coord_names = deepcopy(imas_info("$(f2u(fds)).$(field)")["coordinates"])
    coord_values = []
    for (k, coord) in enumerate(coord_names)
        if contains(coord, "...")
            push!(coord_values, nothing)
        else
            # find common ancestor
            s1 = f2fs(fds)
            s2 = u2fs(p2i(i2p(coord)[1:end - 1]))
            cs, s1, s2 = FUSE.common_base_string(s1, s2)
            common_fds = replace(string(cs), r"___$" => "")
            
            # go upstream until common acestor
            h = fds
            while f2fs(h) != common_fds
                h = h._parent.value
            end
            # then dive into the coordinates branch
            for k in i2p(s2)
                h = getproperty(h, Symbol(k))
            end
            coord_leaf = i2p(coord)[end]
            h = getproperty(h, Symbol(coord_leaf))

            if h === missing
                push!(coord_values, missing)
            else
                push!(coord_values, h)
            end
        end
    end
    return Dict(:names => coord_names, :values => coord_values)
end

"""
    f2u(fds::Union{FDS, FDSvector, DataType, Symbol, String})

Returns universal IMAS location of a given FDS
"""
function f2u(fds::Union{FDS,FDSvector})
    return f2u(typeof(fds))
end

function f2u(fds::DataType)
    return f2u(Base.typename(fds).name)
end

function f2u(fds::Symbol)
    return f2u(string(fds))
end

function f2u(fds::String)
    if in(':', fds) | in('.', fds)
        error("`$fds` is not a qualified FDS type")
    end
    tmp = replace(fds, "___" => "[:].")
    tmp = replace(tmp, "__" => ".")
    imas_location = replace(tmp, r"\.$" => "")
    return imas_location
end

"""
    f2p(fds::Union{FDS,FDSvector})::Vector{Union{String,Int}}

Returns IMAS location of a given FDS

NOTE: indexes of arrays of structures that cannot be determined are set to 0
"""
function f2p(fds::Union{FDS,FDSvector})::Vector{Union{String,Int}}
    return f2p(fds, missing, nothing, Int[])
end

function f2p(fds::Union{FDS,FDSvector},
             child::Union{Missing,FDS,FDSvector},
             path::Union{Nothing,Vector},
             index::Vector{Int})
    # initialize path and index
    if path === nothing
        if typeof(fds) <: FDS
            if typeof(fds._parent.value) <: FDSvector
                name = string(Base.typename(typeof(fds)).name) * "___"
            else
                name = string(Base.typename(typeof(fds)).name)
            end
        elseif typeof(fds) <: FDSvector
            name = string(Base.typename(eltype(fds)).name) * "___"
        end
        path = replace(name, "___" => "__:__")
        path = Vector{Any}(Vector{String}(split(path, "__")))
        path = [k == ":" ? 0 : k for k in path if length(k) > 0]
    end

    # collect integers for arrays of structures
    if typeof(fds) <: FDSvector
        ix = findfirst([k === child for k in fds.value])
        if ix === nothing
            push!(index, 0)
        else
            push!(index, ix)
        end
    end

    # traverse FDSs upstream or return result once top is reached
    if fds._parent.value === missing
        index = reverse(index)
        path = reverse([(typeof(k) <: Int) & (length(index) > 0) ? pop!(index) : k for k in reverse(path)])
        return path
    else
        return f2p(fds._parent.value, fds, path, index)
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
        if contains(s, '[')
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
    return FDS type as a string from IMAS location string
"""
function u2fs(imas_location::String)::String
    tmp = replace(imas_location, "." => "__")
    return replace(tmp, "[:]" => "_")
end

"""
    return FDS/FDSvector from IMAS location string
"""
function u2f(imas_location::String)
    return eval(Meta.parse(u2fs(imas_location)))
end

"""
    return FDS type as a string
"""
function f2fs(fds::Union{FDS,FDSvector})::String
    if typeof(fds) <: FDS
        if typeof(fds._parent.value) <: FDSvector
            return string(Base.typename(typeof(fds)).name)
        else
            return string(Base.typename(typeof(fds)).name)
        end
    elseif typeof(fds) <: FDSvector
        return string(Base.typename(eltype(fds)).name)
    end
end

import Base:keys

"""
    keys(fds::Union{FDS, FDSvector})

Returns list of fields with data in a FDS/FDSvector
"""
function Base.keys(fds::FDS)
    kkk = Symbol[]
    for k in fieldnames(typeof(fds))
        # hide the _parent field
        if k === :_parent
            continue
        end
        v = getfield(fds, k)
        # empty entries
        if v === missing
            continue
        # empty structures/arrays of structures (recursive)
        elseif typeof(v) <: Union{FDS,FDSvector}
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

function Base.keys(fds::FDSvector)
    return collect(1:length(fds))
end

import Base:show

function Base.show(io::IO, fds::Union{FDS,FDSvector}, depth::Int)
    items = keys(fds)
    for (k, item) in enumerate(items)
        printstyled("$('｜'^depth)"; color=:yellow)
        # arrays of structurs
        if typeof(fds) <: FDSvector
            printstyled("[$(item)]\n"; bold=true, color=:green)
            show(io, fds[item], depth + 1)
        # structures
        elseif typeof(getfield(fds, item)) <: Union{FDS,FDSvector}
            if (typeof(fds) <: dd)
                printstyled("$(uppercase(string(item)))\n"; bold=true)
            else
                printstyled("$(string(item))\n"; bold=true)
            end
            show(io, getfield(fds, item), depth + 1)
        # field
        else
            printstyled("$(item)")
            printstyled(" ➡ "; color=:red)
            printstyled("$(Base.summary(getfield(fds, item)))\n"; color=:blue)
        end
        if (typeof(fds) <: dd) & (k < length(items))
            println()
end
    end
end

function Base.show(io::IO, fds::Union{FDS,FDSvector})
    return show(io, fds, 0)
end
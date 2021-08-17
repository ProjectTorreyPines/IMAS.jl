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

function Base.getproperty(fds::FDS, field::Symbol)
    return getfield(fds, field)
end

"""
    coordinates(fds::FDS, field::Symbol)

Returns two lists, one of coordinate names and the other with their values in the data structure
Coordinate value is `nothing` when the data does not have a coordinate
Coordinate value is `missing` if the coordinate is missing in the data structure
"""
function coordinates(fds::FDS, field::Symbol)
    coord_names = deepcopy(imas_info("$(f2i(fds)).$(field)")["coordinates"])
    coord_values = []
    for (k, coord) in enumerate(coord_names)
        if contains(coord, "...")
            push!(coord_values, nothing)
        else
            h = top(fds)
            coord = replace(coord, ":" => "1")
            for k in i2p(coord)
                if (typeof(k) <: Int) & (typeof(h) <: FDSvector)
                    if length(keys(h)) <= k
                        h = h[k]
                        # println('*',k)
                    else
                        # println('-',k)
                    end
                elseif (typeof(k) <: Symbol) & (typeof(h) <: FDS)
                    if hasfield(typeof(h), k)
                        h = getfield(h, k)
                        # println('*',k)
                    else
                        # println('-',k)
                    end
                else
                    # println('-',k)
                end
            end
            if (h === fds) | (h === missing)
                push!(coord_values, missing)
            else
                push!(coord_values, h)
            end
        end
    end
    return coord_names, coord_values
end

function Base.setproperty!(fds::FDS, field::Symbol, v)
    if typeof(v) <: WeakRef
        return setfield!(fds, field, v)
    end

    if typeof(v) <: FDS
        v._parent = WeakRef(fds)
    end

    coords_names, coord_values = coordinates(fds, field)
    for (c_name, c_value) in zip(coords_names, coord_values)
        if c_value === missing
            error("Assign data to `$c_name` before assigning `$(f2i(fds)).$(field)`")
        end
    end

    # if ! (typeof(v) <: AbstractFunctionArray)
    #     target_type = typeintersect(convertsion_types, struct_field_type(typeof(fds), field))
    #     if target_type <: AbstractFunctionArray
    #         imas_info("$(f2i(fds)).$(field)")["coordinates"]
    #         domain = top(fds)
    #         v = NumericalFunctionVector(1:length(v), v)
    #     end
    # end

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
    v._parent = WeakRef(x)
end

import Base: push!, pop!

function push!(x::FDSvector{T}, v) where {T <: FDS}
    v._parent = WeakRef(x)
    push!(x.value, v)
end

function pop!(x::Vector{T}) where {T <: FDS}
    pop!(x.value)
end

#= ============ =#
#  FDSfunctions  #
#= ============ =#

"""
    f2i(fds::Union{FDS, FDSvector, DataType, Symbol, String})

Returns IMAS location of a given FDS
"""
function f2i(fds::Union{FDS,FDSvector})
    return f2i(typeof(fds))
end

function f2i(fds::DataType)
    return f2i(Base.typename(fds).name)
end

function f2i(fds::Symbol)
    return f2i(string(fds))
end

function f2i(fds::String)
    tmp = replace(fds, "___" => "[:].")
    tmp = replace(tmp, "__" => ".")
    imas_location = replace(tmp, r"_$" => "")
    return imas_location
end

"""
    i2p(imas_location::String)

Split IMAS location in its elements
"""
function i2p(imas_location::String)
    path = Any[]
    for s in split(imas_location, '.')
        if contains(s, '[')
            s, n = split(s, '[')
            n = strip(n, ']')
            push!(path, Symbol(s))
            if n == ":"
                push!(path, n)
            else
                push!(path, parse(Int, n))
            end
        else
            push!(path, Symbol(s))
        end
    end
    return path
end

"""
    p2i(path::Any[])

Combine list of IMAS location elements into a string
"""
function p2i(path::Vector{Any})
    str = String[]
    for k in path
        if typeof(k) <: Symbol
            push!(str, string(k))
        elseif typeof(k) <: Int
            push!(str, "[$(string(k))]")
        elseif (k == ":") | (k == ':') | (typeof(k) === Colon) 
            push!(str, "[:]")
        elseif typeof(k) <: String
            push!(str, k)
        end
    end
    return replace(join(str, "."), ".[" => "[")
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
        # arrays of structurs
        if typeof(fds) <: FDSvector
            printstyled("$(' '^depth)[$(item)]\n"; bold=true, color=:green)
            show(io, fds[item], depth + 1)
        # structures
        elseif typeof(getfield(fds, item)) <: Union{FDS,FDSvector}
            if (typeof(fds) <: dd)
                printstyled("$(' '^depth)$(uppercase(string(item)))\n"; bold=true)
            else
                printstyled("$(' '^depth)$(string(item))\n"; bold=true)
            end
            show(io, getfield(fds, item), depth + 1)
        # field
        else
            printstyled("$(' '^depth)$(item)")
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

#= ===================== =#
#  AbstractFunctionArray  #
#= ===================== =#

using Interpolations

abstract type AbstractFunctionArray{T,N} <: AbstractArray{T,N} end
const AFA = AbstractFunctionArray{T,N} where {T,N}

#= ====================== =#
#  NumericalFunctionArray  #
#= ====================== =#

struct NumericalFunctionVector{T} <: AbstractFunctionArray{T,1}
    domain::AbstractVector
    value::AbstractVector{T}
end

Base.broadcastable(x::NumericalFunctionVector) = Base.broadcastable(x.value)

function Base.getindex(x::NumericalFunctionVector, i::Int64)
    x.value[i]
end

Base.size(p::NumericalFunctionVector) = size(p.domain)

function Base.setindex!(x::NumericalFunctionVector, v, i::Int64)
    x.value[i] = v
end

function (x::NumericalFunctionVector)(y)
    LinearInterpolation(x.domain, x.value, )(y)
end

#= ======================= =#
#  AnalyticalFunctionArray  #
#= ======================= =#

struct AnalyticalFunctionVector <: AbstractFunctionArray{Float64,1}
    domain::AbstractVector
    func::Function
end

Base.broadcastable(x::AnalyticalFunctionVector) = Base.broadcastable(x.func(x.domain))

function Base.getindex(x::AnalyticalFunctionVector, i::Int64)
    x.func(x.domain[i])
end

Base.size(p::AnalyticalFunctionVector) = size(p.domain)

function Base.setindex!(x::AnalyticalFunctionVector, v, i::Int64)
    error("Cannot setindex! of a AnalyticalFunctionVector")
end

function (x::AnalyticalFunctionVector)(y)
    x.func(y)
end

#= ======== =#
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

function f2p(ids::Union{IDS,IDSvector}, child::Union{Missing,IDS,IDSvector}, path::Union{Nothing,Vector}, index::Vector{Int})
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
        ix = findfirst([k === child for k in ids._value])
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
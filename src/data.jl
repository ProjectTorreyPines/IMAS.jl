import JSON

include("dd.jl")

import Base: resize!

"""
    resize!(collection::Vector{T}, n::Int) where {T<:FUSE.FDS}

Change size of a FDS array of structures
"""
function resize!(collection::Vector{T}, n::Int) where {T<:FUSE.FDS}
    if n > length(collection)
        for k in length(collection):n - 1
            push!(collection, eltype(collection)())
            # println("push $(length(collection))")
        end
    elseif n < length(collection)
        for k in n:length(collection) - 1
            pop!(collection)
            # println("pop $(length(collection))")
        end
    end
    return collection
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
    dict2fuse(dct, fds::FDS=dd() ;verbose::Bool=false, path::Vector{String}=String[])::FDS

Populate FUSE data structure `fds` based on data contained in Julia dictionary `dct`.

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
"""
function dict2fuse(dct, fds::FDS=dd() ;verbose::Bool=false, path::Vector{String}=String[])::FDS
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
            for i in 1:length(v)
                if verbose println(("｜"^(level + 1)) * string(i)) end
                dict2fuse(v[i], ff[i]; path=vcat(path, [string(k),"[$i]"]), verbose=verbose)
            end

        # Leaf
        else
            if verbose print(("｜"^level) * string(k) * " → ") end
            target_type = typeintersect(supported_types, struct_field_type(typeof(fds), Symbol(k)))
            if target_type <: Array
                if ndims(target_type) == 2
                    v = hcat(v...)'
                end
            end
            v = convert(target_type, v)
            setfield!(fds, Symbol(k), v)
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
    return FUSE.dict2fuse(JSON.parsefile(filename); verbose=verbose)
end
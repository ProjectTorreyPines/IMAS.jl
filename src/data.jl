import JSON

include("dd.jl")

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
            target_type = typeintersect(convertsion_types, struct_field_type(typeof(fds), Symbol(k)))
            if target_type <: Array
                if ndims(target_type) == 2
                    v = reduce(hcat, v)
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
    fds_data = dd()
    json_data = JSON.parsefile(filename)
    dict2fuse(json_data, fds_data; verbose=verbose)
    return fds_data
end

"""
    top(fds::Union{FDS, FDSvector})

Return top-level FDS in the DD hierarchy
"""
function top(fds::Union{FDS,FDSvector})
    if fds._parent.value === missing
        return fds
    else
        return top(fds._parent.value)
    end
end

mutable struct DDigestField
    tag::String
    ids::IMAS.IDS
    field::Symbol
    value::Any
end

function digest(ids::Union{IMAS.IDS,IMAS.IDSvector}; reverse::Bool=true, join_with::String=".")
    out = digest!(Dict{String,Any}(), ids)

    # generate all possible tags
    tmp = Dict{String,Vector{String}}()
    for key in keys(out)
        skey = split(key, ".")
        if skey[end] == "value"
            skey = skey[1:end-1]
        end
        if reverse
            tmp[key] = [join(skey[end:-1:end-k+1], join_with) for k in eachindex(skey)]
        else
            tmp[key] = [join(skey[end-k+1:end], join_with) for k in eachindex(skey)]
        end
    end

    # find shortest tag
    result = Dict{Symbol,Any}()
    l = 1
    while !isempty(out)
        for key1 in keys(out)
            collision = false
            for key2 in keys(out)
                if key1 != key2 && tmp[key1][l] == tmp[key2][l]
                    collision = true
                end
            end
            if !collision
                tag = tmp[key1][l]
                result[Symbol(tag)] = pop!(out, key1)
                result[Symbol(tag)].tag = tag
            end
        end
        l += 1
    end

    return result
end

function Base.show(io::IO, node::DDigestField)
    print(io, node.value)
    u = units(node.ids, node.field)
    if !(isempty(u) || u == "-")
        print(io, " [$u]")
    end
end

function digest!(out::Dict{String,Any}, ids::Union{IMAS.IDS,IMAS.IDSvector})
    for field in keys(ids)
        if typeof(ids) <: IMAS.IDS
            value = getraw(ids, field)
        elseif typeof(ids) <: IMAS.IDSvector
            value = ids[field]
        end
        if value === missing
            continue
        elseif field == :source
            continue
        elseif typeof(value) <: Union{IMAS.IDS,IMAS.IDSvector}
            digest!(out, value)
        else
            fname = string(IMAS.f2i(ids)) * ".$field"
            value = getproperty(ids, field, missing)
            if value !== missing
                out[fname] = DDigestField(fname, ids, field, get_time_array(ids, field))
            end
        end
    end
    return out
end

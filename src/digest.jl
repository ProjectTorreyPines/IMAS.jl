struct DDigestField
    tag::String
    ids::IMAS.IDS
    field::Symbol
    value::Any
end

function digest(ids::Union{IMAS.IDS,IMAS.IDSvector})
    out = digest!(Dict{String,String}(), ids)

    # generate all possible tags
    tmp = Dict{String,Vector{String}}()
    for key in keys(out)
        skey = split(key, ".")
        if skey[end] == "value"
            skey = skey[1:end-1]
        end
        tmp[key] = [join(skey[end:-1:end-k+1], "__") for k in eachindex(skey)]
    end

    # find shortest tag
    result_strings = Dict{String,String}()
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
                result_strings[tmp[key1][l]] = pop!(out, key1)
            end
        end
        l += 1
    end

    # evaluate locations
    result_values = Dict{String,Any}()
    for (key, location) in result_strings
        parent_location = p2i(i2p(location)[2:end-1])
        parent_ids = goto(ids, parent_location; f2=f2fs)
        field = Symbol(i2p(location)[end])

        value = getproperty(parent_ids, field, missing)
        if value !== missing
            result_values[key] = DDigestField(key, parent_ids, field, get_time_array(parent_ids, field))
        end
    end

    return result_values
end

function Base.show(io::IO, node::DDigestField)
    print(io, node.value)
    u = units(node.ids, node.field)
    if !(isempty(u) || u == "-")
        print(io, " [$u]")
    end
end

function digest!(out::Dict{String,String}, ids::Union{IMAS.IDS,IMAS.IDSvector})
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
            out[fname] = fname
        end
    end
    return out
end

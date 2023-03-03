const DDigest = Dict{Symbol,String}()

function add_ddigest(location::String)
    ids_type = u2f(location)
    return add_ddigest(ids_type)
end

function add_ddigest(@nospecialize(ids_type::Type{<:IDS}))
    location = fs2u(ids_type)
    out = Dict{Symbol,String}()
    for (fname, ftype) in zip(fieldnames_(ids_type), fieldtypes_(ids_type))
        floc = "$location.$fname"
        if fname == :source
            continue
        elseif ftype <: IDS
            tmp = add_ddigest(floc)
            for k in keys(tmp)
                if k == :value
                    out[fname] = tmp[k]
                else
                    out[Symbol("$(fname)__$k")] = tmp[k]
                end
            end
        elseif ftype <: IDSvector
            #pass
        else
            out[fname] = "dd." * floc
        end
    end
    return out
end

# empty!(IMAS.DDigest)
# merge!(IMAS.DDigest, add_ddigest("summary.global_quantities"))
# merge!(IMAS.DDigest, add_ddigest("equilibrium.time_slice[:].global_quantities"))
# merge!(IMAS.DDigest, add_ddigest("core_profiles.global_quantities"))

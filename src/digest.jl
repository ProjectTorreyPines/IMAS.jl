abstract type AbstractIDSdict end

Base.@kwdef mutable struct DDigest <: AbstractIDSdict
    _data::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

function Base.getproperty(igest::AbstractIDSdict, field::Symbol)
    data = getfield(igest,:_data)
    if field ∉ keys(data)
        data[field] = typeof(igest)()
    else
        return data[field]
    end
end

function Base.setproperty!(igest::AbstractIDSdict, field::Symbol, value::Any)
    data = getfield(igest, :_data)
    if typeof(value) <: IDS
        data[field] = typeof(igest)()
        for (subfield,subvalue) in value
            setproperty!(data[field], subfield, subvalue)
        end
    elseif 
        data[field] = value
    end
end

function Base.show(io::IO, igest::AbstractIDSdict)
    return show(io, getfield(igest, :_data))
end

function Base.keys(igest::AbstractIDSdict)
    return keys(getfield(igest,:_data))
end

#==========================#

function digest(dd::IMAS.dd)
    igest = DDigest()
    if !isempty(dd.equilibrium.time_slice)
        for (field,value) in dd.equilibrium.time_slice[].global_quantities
            setproperty!(igest.equilibrium, field, value)
        end
    end

    if !isempty(dd.core_profiles.global_quantities)
        for (field,value) in dd.core_profiles.global_quantities
            setproperty!(igest.equilibrium, field, value)
        end
    end

    return igest
end

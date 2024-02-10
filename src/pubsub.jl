import Jedis
import FXP

mutable struct FxpDD
    enabled::Bool
    identifier::String
    client::Jedis.Client
end

function fxp_connect(dd::IMAS.dd; host::String="localhost", port::Int=55000, password::String="redispw")
    aux = getfield(dd, :_aux)
    client = Jedis.Client(; host, port, password, keepalive_enable=true)
    identifier = string(objectid(dd))
    b = Base.text_colors[:bold]
    n = Base.disable_text_style[:bold]
    @info("fxp @ $(b)HOST$(n):$host, $(b)PORT$(n):$port, $(b)SESSION$(n):$identifier")
    aux[:fxp] = FxpDD(true, identifier, client)
end

function fxp_enable(dd::IMAS.dd, enable_switch::Bool=true)
    aux = getfield(dd, :_aux)
    if :fxp ∈ keys(aux)
        aux[:fxp].enabled = enable_switch
    else
        return false
    end
    return enable_switch
end

function is_fxping(dd::IMAS.dd)
    aux = getfield(dd, :_aux)
    if :fxp ∉ keys(aux)
        return false
    else
        client = fxp_client(dd)
        return aux[:fxp].enabled && !Jedis.isclosed(client)
    end
end

function fxp_client(dd)
    aux = getfield(dd, :_aux)
    if :fxp ∈ keys(aux)
        return aux[:fxp].client
    else
        return nothing
    end
end

function fxp_identifier(dd::IMAS.dd)
    aux = getfield(dd, :_aux)
    if :fxp ∈ keys(aux)
        return aux[:fxp].identifier
    else
        return string(objectid(dd))
    end
end

function dd_service(dd::IMAS.dd, service_name::String)
    return "$(fxp_identifier(dd))__$(service_name)"
end

function Base.deepcopy_internal(x::IMAS.dd{T}, dict::IdDict) where {T<:Real}
    if haskey(dict, x)
        return dict[x]
    end

    # Use the default `deepcopy_internal` to handle the copying process.
    # This creates a shallow copy of `x`, so all mutable fields will
    # still be copied correctly later.
    new_obj = IMAS.dd{T}()

    # Register the new object in `dict` to handle cyclic references.
    dict[x] = new_obj

    # Explicitly perform a deep copy on all fields using the default behavior.
    # Since `copy` has already been used, we don't need to copy primitive fields,
    # but we do need to ensure deep copying of mutable fields.
    fields = fieldnames(typeof(x))
    for field in fields
        field_val = getfield(x, field)
        copied_val = Base.deepcopy_internal(field_val, dict)
        setfield!(new_obj, field, copied_val)
    end

    if :fxp ∈ keys(getfield(x, :_aux))
        getfield(new_obj, :_aux)[:fxp] = getfield(x, :_aux)[:fxp]
    end

    return new_obj
end

function fxp_push!(dd::IMAS.dd, service_name::String; data...)
    client = fxp_client(dd)
    dd_service_name = dd_service(dd, service_name)
    return push!(client, dd_service_name; data...)
end

function fxp_pop!(dd::IMAS.dd, service_name::String; timeout::Float64)
    client = fxp_client(dd)
    dd_service_name = dd_service(dd, service_name)
    return pop!(client, dd_service_name; timeout)
end

function fxp_has_service_provider(dd::IMAS.dd, service_name::String)
    aux = getfield(dd, :_aux)
    if :fxp ∈ keys(aux)
        client = fxp_client(dd)
        return FXP.has_service_provider(service_name; client)
    else
        return false
    end
end

function fxp_negotiate_service(dd::IMAS.dd, service_name::String)
    client = fxp_client(dd)
    return FXP.negotiate_service(fxp_identifier(dd), service_name; client)
end

function fxp_request_service(dd::IMAS.dd, service_name::String)
    if is_fxping(dd) && fxp_has_service_provider(dd, service_name)
        fxp_negotiate_service(dd, service_name)
        return true
    end
    return false
end
import Jedis
import JSON

mutable struct StreamDD
    enabled::Bool
    identifier::String
    client::Jedis.Client
end

function stream(dd::IMAS.dd, enable_switch::Bool=true; host::String="localhost", port::Int=55000, password::String="redispw")
    aux = getfield(dd, :_aux)
    if :stream ∈ keys(aux)
        aux[:stream].enabled = enable_switch
    elseif enable_switch
        client = Jedis.Client(; host, port, password, keepalive_enable=true)
        identifier = string(objectid(dd))
        b = Base.text_colors[:bold]
        n = Base.disable_text_style[:bold]
        @info("dd stream @ $(b)HOST$(n):$host, $(b)PORT$(n):$port, $(b)CHANNEL$(n):$identifier")
        aux[:stream] = StreamDD(enable_switch, identifier, client)
    end
    return enable_switch
end

function is_streaming(dd::IMAS.dd)
    aux = getfield(dd, :_aux)
    if :stream ∉ keys(aux)
        return false
    else
        client = stream_client(dd)
        return aux[:stream].enabled && !Jedis.isclosed(client)
    end
end

function stream_client(dd)
    aux = getfield(dd, :_aux)
    if :stream ∈ keys(aux)
        return aux[:stream].client
    else
        return nothing
    end
end

function stream_identifier(dd::IMAS.dd)
    aux = getfield(dd, :_aux)
    if :stream ∈ keys(aux)
        return aux[:stream].identifier
    else
        return string(objectid(dd))
    end
end

function stream_name(dd::IMAS.dd, name::String)
    return "$(stream_identifier(dd))__$(name)"
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

    if :stream ∈ keys(getfield(x, :_aux))
        getfield(new_obj, :_aux)[:stream] = getfield(x, :_aux)[:stream]
    end

    return new_obj
end

########
# PUSH #
########
function Base.push!(client::Jedis.Client, stream::String; data...)
    return Jedis.execute(["XADD", stream, "*", "data", JSON.sprint(data)], client)
end

function stream_push!(dd::IMAS.dd, name::String; data...)
    client = stream_client(dd)
    dd_stream_name = stream_name(dd, name)
    return push!(client, dd_stream_name; data...)
end

#######
# POP #
#######
function Base.pop!(client::Jedis.Client, stream::String; timeout::Float64, error_on_timeout::Bool=true)
    ms_timeout = Int(round(timeout * 1000))
    out = Jedis.execute(String["XREAD", "BLOCK", string(ms_timeout), "COUNT", "1", "STREAMS", stream, "0-0"], client)
    if isempty(out)
        if error_on_timeout
            error("Reading data from stream `$stream` has timed out")
        else
            return nothing
        end
    end
    payload = out[1][2][1][2]
    @assert payload[1] == "data"
    data = JSON.parse(payload[2])
    delkey = out[1][2][1][1]
    Jedis.execute(["XDEL", stream, delkey], client)
    return Dict(Symbol(k) => v for (k,v) in data)
end

function stream_pop!(dd::IMAS.dd, name::String; timeout::Float64)
    client = stream_client(dd)
    dd_stream_name = stream_name(dd, name)
    return pop!(client, dd_stream_name; timeout)
end

####################
# REGISTER SERVICE #
####################
function stream_register_service(dd_stream_name::String; client::Jedis.Client)
    if dd_stream_name ∈ client.subscriptions
        return nothing
    else
        serviced = []
        @async Jedis.subscribe(dd_stream_name; client) do msg
            push!(serviced, msg)
        end
        Jedis.wait_until_subscribed(client)
        return serviced
    end
end

function stream_register_service(dd::IMAS.dd, name::String)
    client = pubsub_client(dd)
    dd_stream_name = stream_name(dd, name)
    return stream_register_service(dd_stream_name; client)
end

###########################
# CHECK SERVICE PROVIDERS #
###########################
function stream_has_service_provider(dd::IMAS.dd, name::String)
    subs = subscribers(dd, name)
    @assert length(subs) <= 1 "Too many service providers: $(subs)"
    return length(subs) == 1
end

function subscribers(dd::IMAS.dd, name::String)
    aux = getfield(dd, :_aux)
    if :stream ∈ keys(aux)
        client = stream_client(dd)
        dd_stream_name = stream_name(dd, name)
        return Jedis.execute("PUBSUB CHANNELS $(dd_stream_name)*", client)
    else
        return []
    end
end

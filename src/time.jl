"""
    set_timedep_value!(time_ids::IDS, value_ids::IDS, value_symbol::Symbol, time0::Real, value0::Real; time_symbol::Symbol=:time)

Set the value of a time-depentent quantity
"""
function set_timedep_value!(time_ids::IDS, value_ids::IDS, value_symbol::Symbol, time0::Real, value0::Real; time_symbol::Symbol = :time)
    if is_missing(value_ids, value_symbol) || is_missing(time_ids, time_symbol)
        setproperty!(time_ids, time_symbol, [time0])
        setproperty!(value_ids, value_symbol, [value0])
    else
        time = getproperty(time_ids, time_symbol)
        value = getproperty(value_ids, value_symbol)
        if any(time .== time0)
            time_index = argmin(abs.(time .- time0))
            value[time_index] = value0
        elseif time0 > maximum(time)
            push!(time, time0)
            push!(value, value0)
        elseif time0 < minimum(time)
            pushfirst!(time, time0)
            pushfirst!(value, value0)
        else
            for (k, t) in enumerate(time)
                if t < time0
                    continue
                else
                    setproperty!(time_ids, time_symbol, vcat(time[1:k-1], time0, time[k:end]))
                    setproperty!(value_ids, value_symbol, vcat(value[1:k-1], value0, value[k:end]))
                    break
                end
            end
        end
    end
    return getproperty(time_ids, time_symbol), getproperty(value_ids, value_symbol)
end

"""
    look for time array information
"""
function time_array(x::Union{IDS,IDSvector{T}}) where {T<:IDSvectorElement}
    time_array = nothing
    missing_time_locations = []

    # traverse IDS hierarchy upstream looking for a time array
    h = x
    while h._parent.value !== missing
        field_names_types = NamedTuple{fieldnames(typeof(h))}(fieldtypes(typeof(h)))
        if :time in keys(field_names_types) && typeintersect(field_names_types[:time], AbstractVector) !== Union{}
            if is_missing(h, :time)
                push!(missing_time_locations,h)
            elseif length(h.time) == 0
                push!(missing_time_locations,h)
            else
                time_array = h.time
                break
            end
        end
        h = h._parent.value
    end
    if time_array === nothing
        if length(missing_time_locations)==0
            throw("Could not find time array information for $(p2i(f2p(x)[1:end-1]))[$time]")
        else
            time_array = missing_time_locations[end].time = Float64[]
        end
    end
    time_array
end

function global_time(ids::Union{IDS,IDSvector})::Real
    dd = top_dd(ids)
    if dd === missing
        throw("Could not reach top level dd where global time is defined")
    end
    return dd.global_time
end

function set_time_array(ids, location,value)
    time = time_array(ids)
    time0 = global_time(ids)
    i = nothing
    if (length(time) == 0) || (time0 < minimum(time))
        pushfirst!(time, time0)
        if location !== :time
            setproperty!(ids, location, [value])
        end
        i = 1
    else
        i = argmin(abs.(time .- time0))
        # perfect match --> overwrite
        if minimum(abs.(time .- time0)) == 0
            if location !== :time
                getproperty(ids, location)[i] = value
            end
        else
            insert!(time, i+1, time0)
            if location !== :time
                insert!(getproperty(ids, location), i+1, value)
            end
        end
    end
    i = argmin(abs.(time .- time0))
    return getproperty(ids, location)[i]
end

function get_time_array(ids, location)
    time = time_array(ids)
    time0 = global_time(ids)
    i = argmin(abs.(time .- time0))
    return getproperty(ids, location)[i]
end

macro timedep(ex)
    quote
        local expr = $(Meta.QuoteNode(ex))
        if expr.head == :(=)
            local value = $(esc(ex.args[2]))
            local ids = $(esc(ex.args[1].args[1]))
            local location = $(esc(ex.args[1].args[2]))
            local tmp = set_time_array(ids, location, value)
        else
            local ids = $(esc(ex.args[1]))
            local location = $(esc(ex.args[2]))
            local tmp = get_time_array(ids, location)
        end
        tmp
    end
end

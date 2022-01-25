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

function time_array(x::IDSvector{T}; raise_errors::Bool=true) where {T<:IDSvectorElement}
    time_array = []
    h = x
    while h._parent.value !== missing
        h = h._parent.value
        if hasfield(typeof(h), :time)
            if is_missing(h, :time)
                if raise_errors
                    throw("$(f2i(h)).time is not set and $(p2i(f2p(x)[1:end-1]))[$time] could not be determined")
                else
                    h.time = Real[]
                end
            elseif length(h.time) == 0
                if raise_errors
                    throw("$(f2i(h)).time is empty and $(p2i(f2p(x)[1:end-1]))[$time] could not be determined")
                end
            end
            time_array = h.time
            break
        end
    end
    if time_array === nothing
        if raise_errors
            throw("Could not find time array information for $(p2i(f2p(x)[1:end-1]))[$time]")
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
"""
    look for time array information
"""
function time_array(ids::Union{IDS,IDSvector{T}}) where {T<:IDSvectorElement}
    time_array = nothing
    missing_time_locations = []

    # traverse IDS hierarchy upstream looking for a time array
    h = ids
    while h._parent.value !== missing
        field_names_types = NamedTuple{fieldnames(typeof(h))}(fieldtypes(typeof(h)))
        if :time in keys(field_names_types) && typeintersect(field_names_types[:time], AbstractVector) !== Union{}
            if is_missing(h, :time)
                push!(missing_time_locations, h)
            elseif length(h.time) == 0
                push!(missing_time_locations, h)
            else
                time_array = h.time
                break
            end
        end
        h = h._parent.value
    end
    if time_array === nothing
        if length(missing_time_locations) == 0
            error("Could not find time array information for $(p2i(f2p(ids)[1:end-1]))[$time]")
        else
            time_array = missing_time_locations[end].time = Float64[]
        end
    end
    time_array
end

function global_time(ids::Union{IDS,IDSvector})::Real
    dd = top_dd(ids)
    if dd === missing
        error("Could not reach top level dd where global time is defined")
    end
    return dd.global_time
end

"""
    set_time_array(ids, location,value)

set data to a time-dependent array at the dd.global_time
"""
function set_time_array(ids::Union{IDS,IDSvector{T}}, location::Symbol, value) where {T<:IDSvectorElement}
    time = time_array(ids)
    time0 = global_time(ids)
    # no time information
    if length(time) == 0
        push!(time, time0)
        if location !== :time
            setproperty!(ids, location, [value])
        end
    else
        i = argmin(abs.(time .- time0))
        # perfect match --> overwrite
        if minimum(abs.(time .- time0)) == 0
            if location !== :time
                if is_missing(ids, location) || (length(getproperty(ids, location)) == 0)
                    setproperty!(ids, location, vcat([NaN for k = 1:i-1], value))
                else
                    last_value = getproperty(ids, location)
                    if length(last_value)<i
                        reps = i - length(last_value) - 1
                        append!(last_value, vcat([last_value[end] for k = 1:reps], value))
                    else
                        last_value[i] = value
                    end
                end
            end
        # append
        elseif time0 > maximum(time)
            push!(time, time0)
            if location !== :time
                if is_missing(ids, location) || (length(getproperty(ids, location)) == 0)
                    setproperty!(ids, location, vcat([NaN for k = 1:length(time)-1], value))
                else
                    last_value = getproperty(ids, location)
                    reps = length(time) - length(last_value) - 1
                    append!(last_value, vcat([last_value[end] for k = 1:reps], value))
                end
            end
        else
            error("Could not add time array information for $(f2i(ids)).$location[$time]")
        end
    end
    i = argmin(abs.(time .- time0))
    return getproperty(ids, location)[i]
end

"""
    get_time_array(ids, location)

get data from a time-dependent array at the dd.global_time
"""
function get_time_array(ids::Union{IDS,IDSvector{T}}, location::Symbol) where {T<:IDSvectorElement}
    time = time_array(ids)
    time0 = global_time(ids)
    i = argmin(abs.(time .- time0))
    return getproperty(ids, location)[i]
end

function _timedep(ex)
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

"""
    timedep

Macro for getting/setting data of a time-dependent array at the dd.global_time
"""
macro timedep(ex)
    return _timedep(ex)
end

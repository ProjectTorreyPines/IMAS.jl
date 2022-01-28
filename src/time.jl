"""
    time_locations(ids::Union{IDS,IDSvector{T}}) where {T<:IDSvectorElement}

Return all possible time locations with info about whether they are arrays or not
"""
function time_locations(ids::Union{IDS,IDSvector{T}}) where {T<:IDSvectorElement}
    time_locations = []
    time_array = []

    # traverse IDS hierarchy upstream looking for a :time
    h = ids
    while h._parent.value !== missing
        field_names_types = NamedTuple{fieldnames(typeof(h))}(fieldtypes(typeof(h)))
        if :time in keys(field_names_types)
            pushfirst!(time_locations, h)
            pushfirst!(time_array, typeintersect(field_names_types[:time], AbstractVector) !== Union{})
        end
        h = h._parent.value
    end

    return time_locations, time_array
end

"""
    time_array(ids::Union{IDS,IDSvector{T}}) where {T<:IDSvectorElement}

Look for time array information
"""
function time_array(ids::Union{IDS,IDSvector{T}}) where {T<:IDSvectorElement}
    locs, tarr = time_locations(ids)
    h = [locs[k] for k in 1:length(locs) if tarr[k]][end]
    if is_missing(h, :time)
        h.time = Float64[]
    end
    return h.time
end

"""
    global_time(ids::Union{IDS,IDSvector})::Real

get the dd.global_time of a given IDS
"""
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

"""
    ddtime

Macro for getting/setting data of a time-dependent array at the dd.global_time
"""
macro ddtime(ex)
    return _ddtime(ex)
end

function _ddtime(ex)
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

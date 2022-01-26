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
    insert_time_index(time_array, time; overwrite = false)

return insertion point of time to keep time_array sorted

NOTE: Negative times are treated as "special" and are appended at the end
"""
function insert_time_index(time_array, time; overwrite = false)
    if time < 0
        t0 = 1E6 + time
    end
    for (k, t) in enumerate(time_array)
        if t == time
            if overwrite
                return k
            else
                throw("time slice $time already present: time=$time_array")
            end
        elseif t > t0
            return k
        end
    end
    return length(time_array) + 1
end

"""
    insert_time_slice!(dd::IMAS.dd, args...; kw...)

Insert a time slice of an IDS in the data dictionary keeping time arrays sorted
"""
function insert_time_slice!(dd::IMAS.dd, args...; kw...)
    return insert_time_slice!(getproperty(dd, Symbol(IMAS.f2p(args[1])[1])), args...; kw...)
end

function insert_time_slice!(eq::IMAS.equilibrium, eqt::IMAS.equilibrium__time_slice; time = nothing, R0 = nothing, overwrite = false)
    if time === nothing
        time = eqt.time
    end
    if R0 === nothing
        R0 = eq.vacuum_toroidal_field.r0
    end
    if IMAS.is_missing(eq, :time)
        eq.time = Real[]
        eq.vacuum_toroidal_field.b0 = Real[]
    end
    # find time_index of insertion
    time_index = insert_time_index(eq.time, time; overwrite = overwrite)

    # set R0
    eq.vacuum_toroidal_field.r0 = R0
    # set time
    eqt.time = time

    if overwrite && (length(eq.time) >= time_index)
        # update [eq.time]
        eq.time[time_index] = time
        # update [b0]
        eq.vacuum_toroidal_field.b0[time_index] = eqt.profiles_1d.f[end] / R0
        # update [eq.time_slice]
        eq.time_slice[time_index] = eqt
    else
        # insert [eq.time]
        insert!(eq.time, time_index, time)
        # insert [b0]
        insert!(eq.vacuum_toroidal_field.b0, time_index, eqt.profiles_1d.f[end] / R0)
        # insert [eq.time_slice]
        insert!(eq.time_slice, time_index, eqt)
    end
    return eq
end

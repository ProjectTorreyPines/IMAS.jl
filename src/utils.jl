"""
    set_field_time_array(ids::IDS, field::Symbol, time_index::Integer, value::Real)

Set the value at a given time_index of a time dependent array in a given IDS
"""
function set_field_time_array(ids::IDS, field::Symbol, time_index::Integer, value::Real)
    array = missing
    try
        array = getproperty(ids, field)
    catch
        # pass
    end
    if (array === missing)
        if time_index == 1
            return setproperty!(ids, field, [value])
        else
            error("$(f2i(ids)).$(field) is empty: cannot set data at index $(time_index)")
        end
    end
    if time_index <= length(array)
        array[time_index] = value
    elseif time_index == (length(array) + 1)
        push!(array, value)
    else
        error("$(f2i(ids)).$(field) has length $(length(array)): cannot set data at index $(time_index)")
    end
    return setproperty!(ids, field, array)
end

"""
    get_time_index(ids_time_slice::IDS, time::Real)::Integer

Return index of a given time and resize array of structures accordingly if necessary
"""
function get_time_index(ids_time_slice::IDSvector, time::Real)::Integer
    ids = top(ids_time_slice)
    if is_missing(ids, :time)
        time_index = 1
        set_field_time_array(ids, :time, time_index, time)
        resize!(ids_time_slice, time_index)
    elseif time in ids.time
        time_index = findall(x -> x == time, ids.time)[1]
    elseif maximum(ids.time) < time
        time_index = length(ids.time) + 1
        set_field_time_array(ids, :time, time_index, time)
        resize!(ids_time_slice, time_index)
    else
        error("Cannot append time slice of $(f2i(ids_time_slice)) at $(time) seconds")
    end
    return time_index
end


"""
    common_base_string(s1::String, s2::String)::Vector{String}

given two strings it returns a tuple of 3 strings that is the common initial part, and then the remaining parts
"""
function common_base_string(s1::String, s2::String)::Tuple{String,String,String}
    index = nothing
    for k in 1:min(length(s1), length(s2))
        sub = SubString(s2, 1, k)
        if startswith(s1, sub)
            index = k
        end
    end
    if index === nothing
        return "", s1, s2
    else
        return string(SubString(s1, 1, index)), string(SubString(s1, index + 1, length(s1))), string(SubString(s2, index + 1, length(s2)))
    end
end

"""
    is_missing(ids::IDS, leaf)::Bool

returns true/false if field is missing in IDS
"""
function is_missing(ids::IDS, field)::Bool
    try
        getproperty(ids, field)
        return false
    catch
        return true
    end
end

"""
    iscallable(f)

returns true if argument is callable
"""
function iscallable(f)
    return !isempty(methods(f))
end

"""
    norm(x::Vector{T} where T<:Real)::Vector{Real}

Normalize a vector so that the first item in the array is 0 and the last one is 1
This is handy where psi_norm should be used (and IMAS does not define a psi_norm array)
"""
function norm(x::Vector{T} where T<:Real)::Vector{Real}
    return (x-x[1])/(x[end]-x[1])
end
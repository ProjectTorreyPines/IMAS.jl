"""
    get_from(dd::IMAS.dd, what::Symbol, from_where::Symbol)

Gets `what` from `where` from `dd`, for example

    get_from(dd, :ip, :equilibrium)

IMAS stores the same physical quantities in different IDSs.
This function can used to abstract where information should come from,
and is generally handy when coupling different codes/modules/actors.
"""
function get_from(dd::IMAS.dd, what::Symbol, from_where::Symbol)
    return get_from(dd, Val{what}, from_where)
end

# ip [A]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:ip}}, from_where::Symbol)::T where {T<:Real}
    if from_where == :equilibrium
        return dd.equilibrium.time_slice[].global_quantities.ip
    elseif from_where == :core_profiles
        return @ddtime(dd.core_profiles.global_quantities.ip)
    elseif from_where == :pulse_schedule
        return @ddtime(dd.pulse_schedule.flux_control.i_plasma.reference.data)
    else
        error("$what from $from_where doesn't exist yet")
    end
end

# beta_normal [-]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:beta_normal}}, from_where::Symbol)::T where {T<:Real}
    if from_where == :equilibrium
        return dd.equilibrium.time_slice[].global_quantities.beta_normal
    elseif from_where == :core_profiles
        return @ddtime(dd.core_profiles.global_quantities.beta_tor_norm)
    else
        error("$what from $from_where doesn't exist yet")
    end
end

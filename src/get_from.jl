#===
IMAS stores the same physical quantities in different IDSs.
These sets of functions can used to abstract where information should come from,
and is generally handy when coupling different codes/modules/actors.
===#

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
function get_from(dd::IMAS.dd{T}, what::Type{Val{:βn}}, from_where::Symbol)::T where {T<:Real}
    if from_where == :equilibrium
        return dd.equilibrium.time_slice[].global_quantities.beta_normal
    elseif from_where == :core_profiles
        return @ddtime(dd.core_profiles.global_quantities.beta_tor_norm)
    else
        error("$what from $from_where doesn't exist yet")
    end
end

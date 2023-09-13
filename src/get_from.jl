#===
IMAS stores the same physical quantities in different IDSs.
These sets of functions can used to abstract where information should come from,
and is generally handy when coupling different codes/modules/actors.
===#

# ip [A]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:ip}}, from_where::Symbol; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :equilibrium
        return dd.equilibrium.time_slice[time0].global_quantities.ip
    elseif from_where == :core_profiles
        return IMAS.get_time_array(dd.core_profiles.global_quantities, :ip, time0, :linear)
    elseif from_where == :pulse_schedule
        return IMAS.get_time_array(dd.pulse_schedule.flux_control.i_plasma.reference, :data, time0, :linear)
    else
        error("$what from $from_where doesn't exist yet")
    end
end

# beta_normal [-]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:βn}}, from_where::Symbol; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :equilibrium
        return dd.equilibrium.time_slice[time0].global_quantities.beta_normal
    elseif from_where == :core_profiles
        return IMAS.get_time_array(dd.core_profiles.global_quantities, :beta_tor_norm, time0, :linear)
    else
        error("$what from $from_where doesn't exist yet")
    end
end

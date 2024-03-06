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
        return IMAS.get_time_array(dd.pulse_schedule.flux_control.i_plasma, :reference, time0, :linear)
    else
        error("`get_from(dd, $what, $from_where)` doesn't exist yet")
    end
end

# vacuum_r0_b0 [m], [T]
function get_from(dd::IMAS.dd, what::Type{Val{:vacuum_r0_b0}}, from_where::Symbol; time0::Float64=dd.global_time)
    if from_where == :equilibrium
        eqt = dd.equilibrium.time_slice[time0]
        return (r0=eqt.global_quantities.vacuum_toroidal_field.r0, b0=eqt.global_quantities.vacuum_toroidal_field.b0)
    elseif from_where == :pulse_schedule
        return (r0=dd.pulse_schedule.tf.r0, b0=IMAS.get_time_array(dd.pulse_schedule.tf.b_field_tor_vacuum, :reference, time0, :linear))
    else
        error("`get_from(dd, $what, $from_where)` doesn't exist yet")
    end
end

# vloop [V]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:vloop}}, from_where::Symbol; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :equilibrium
        return vloop(dd.equilibrium; time0)
    elseif from_where == :core_profiles
        return vloop(dd.core_profiles.profiles_1d[time0], dd.equilibrium.time_slice[time0])
    elseif from_where == :pulse_schedule
        return IMAS.get_time_array(dd.pulse_schedule.flux_control.loop_voltage, :reference, time0, :linear)
    elseif from_where == :controllers__ip
        return vloop(dd.controllers; time0)
    else
        error("`get_from(dd, $what, $from_where)` doesn't exist yet")
    end
end

# beta_normal [-]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:βn}}, from_where::Symbol; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :equilibrium
        return dd.equilibrium.time_slice[time0].global_quantities.beta_normal
    elseif from_where == :core_profiles
        return IMAS.get_time_array(dd.core_profiles.global_quantities, :beta_tor_norm, time0, :linear)
    else
        error("`get_from(dd, $what, $from_where)` doesn't exist yet")
    end
end

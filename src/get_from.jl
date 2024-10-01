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
        if time0 >= dd.core_profiles.time[end]
            return Ip(dd.core_profiles.profiles_1d[end])
        else
            return IMAS.get_time_array(dd.core_profiles.global_quantities, :ip, time0, :linear)
        end
    elseif from_where == :pulse_schedule
        return IMAS.get_time_array(dd.pulse_schedule.flux_control.i_plasma, :reference, time0, :linear)
    end
    return error("`get_from(dd, $what, Val{:$from_where})` doesn't exist yet")
end

# vacuum_r0_b0 [m], [T]
function get_from(dd::IMAS.dd, what::Type{Val{:vacuum_r0_b0}}, from_where::Symbol; time0::Float64=dd.global_time)
    if from_where == :equilibrium
        eqt = dd.equilibrium.time_slice[time0]
        return (r0=eqt.global_quantities.vacuum_toroidal_field.r0, b0=eqt.global_quantities.vacuum_toroidal_field.b0)
    elseif from_where == :pulse_schedule
        return (r0=dd.pulse_schedule.tf.r0, b0=IMAS.get_time_array(dd.pulse_schedule.tf.b_field_tor_vacuum, :reference, time0, :linear))
    end
    return error("`get_from(dd, $what, Val{:$from_where})` doesn't exist yet")
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
    end
    return error("`get_from(dd, $what, Val{:$from_where})` doesn't exist yet")
end

# beta_normal [-]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:βn}}, from_where::Symbol; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :equilibrium
        return dd.equilibrium.time_slice[time0].global_quantities.beta_normal
    elseif from_where == :core_profiles
        if time0 >= dd.core_profiles.time[end]
            return beta_tor_norm(dd.equilibrium, dd.core_profiles.profiles_1d[end])
        else
            return IMAS.get_time_array(dd.core_profiles.global_quantities, :beta_tor_norm, time0, :linear)
        end
    end
    return error("`get_from(dd, $what, Val{:$from_where})` doesn't exist yet")
end

# ne_ped w/ pedestal fit [m^-3]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:ne_ped}}, from_where::Symbol, rho_ped::Nothing; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :core_profiles
        cp1d = dd.core_profiles.profiles_1d[time0]
        _, w_ped = IMAS.pedestal_finder(cp1d.electrons.density_thermal, cp1d.grid.psi_norm)
        rho_ped = 1.0 - w_ped
    else
        rho_ped = NaN
    end
    ne_ped = get_from(dd, what, from_where, rho_ped; time0)
    return ne_ped
end

# ne_ped [m^-3]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:ne_ped}}, from_where::Symbol, rho_ped::Float64; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :summary
        return IMAS.get_time_array(dd.summary.local.pedestal.n_e, :value, time0, :linear)
    elseif from_where == :core_profiles
        cp1d = dd.core_profiles.profiles_1d[time0]
        return IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(rho_ped)
    elseif from_where == :pulse_schedule
        if !ismissing(dd.pulse_schedule.density_control.n_e_pedestal, :reference)
            return IMAS.get_time_array(dd.pulse_schedule.density_control.n_e_pedestal, :reference, time0, :linear)
        elseif !ismissing(dd.pulse_schedule.density_control.n_e_pedestal_greenwald_fraction, :reference)
            ne_ped_frac = IMAS.get_time_array(dd.pulse_schedule.density_control.n_e_pedestal_greenwald_fraction, :reference, time0, :linear)
            ne_gw = IMAS.greenwald_density(dd.pulse_schedule; time0) # must evaluate the greenwald density from pulse_schedule!
            return ne_ped_frac * ne_gw
        elseif !ismissing(dd.pulse_schedule.density_control.n_e_greenwald_fraction, :reference)
            ne_frac = IMAS.get_time_array(dd.pulse_schedule.density_control.n_e_greenwald_fraction, :reference, time0, :linear)
            ne_gw = IMAS.greenwald_density(dd.pulse_schedule; time0) # must evaluate the greenwald density from pulse_schedule!            
            return ne_frac * ne_gw
        end
        error("`get_from(dd, $what, Val{:$from_where})` does not have data")
    end
    return error("`get_from(dd, $what, Val{:$from_where})` doesn't exist yet")
end

# zeff_ped w/ pedestal fit [-]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:zeff_ped}}, from_where::Symbol, rho_ped::Nothing; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :core_profiles
        cp1d = dd.core_profiles.profiles_1d[time0]
        _, w_ped = IMAS.pedestal_finder(cp1d.electrons.density_thermal, cp1d.grid.psi_norm)
        rho_ped = 1.0 - w_ped
    else
        rho_ped = NaN
    end
    return get_from(dd, what, from_where, rho_ped; time0)
end

# zeff_ped [-]
function get_from(dd::IMAS.dd{T}, what::Type{Val{:zeff_ped}}, from_where::Symbol, rho_ped::Float64; time0::Float64=dd.global_time)::T where {T<:Real}
    if from_where == :summary
        return IMAS.get_time_array(dd.summary.local.pedestal.zeff, :value, time0, :linear)
    elseif from_where == :core_profiles
        cp1d = dd.core_profiles.profiles_1d[time0]
        return IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.zeff).(rho_ped)
    elseif from_where == :pulse_schedule
        if !ismissing(dd.pulse_schedule.density_control.zeff_pedestal, :reference)
            return IMAS.get_time_array(dd.pulse_schedule.density_control.zeff_pedestal, :reference, time0, :linear)
        elseif !ismissing(dd.pulse_schedule.density_control.zeff, :reference)
            return IMAS.get_time_array(dd.pulse_schedule.density_control.zeff_pedestal, :reference, time0, :linear)
        end
        error("`get_from(dd, $what, Val{:$from_where})` does not have data")
    end
    return error("`get_from(dd, $what, Val{:$from_where})` doesn't exist yet")
end

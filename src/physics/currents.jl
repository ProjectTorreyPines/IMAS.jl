document[Symbol("Physics currents")] = Symbol[]

function J_tor(eqt1d::IMAS.equilibrium__time_slice___profiles_1d)
    return @. (-(eqt1d.dpressure_dpsi + eqt1d.f_df_dpsi * eqt1d.gm1 / mks.μ_0) * (2π)) / eqt1d.gm9
end

"""
    j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice{T}, cp1d::IMAS.core_profiles__profiles_1d{T}, ip::T, j_ohmic_shape::AbstractVector{T}=cp1d.conductivity_parallel) where {T<:Real}

Sets `j_ohmic` parallel current density to what it would be at steady-state, based on `conductivity_parallel`, `j_non_inductive` and a target total `ip`

Requires constant loop voltage: `Vl = 2π * η * <J_oh⋅B> / (F * <R⁻²>) = constant`
"""
function j_ohmic_steady_state(
    eqt::IMAS.equilibrium__time_slice{T},
    cp1d::IMAS.core_profiles__profiles_1d{T},
    ip::T,
    j_ohmic_shape::AbstractVector{T}=cp1d.conductivity_parallel
) where {T<:Real}
    rho_tor_norm = cp1d.grid.rho_tor_norm
    rho_eq = eqt.profiles_1d.rho_tor_norm
    F = cubic_interp1d(rho_eq, eqt.profiles_1d.f).(rho_tor_norm)
    gm1 = cubic_interp1d(rho_eq, eqt.profiles_1d.gm1).(rho_tor_norm) # <R⁻²>
    j_oh_par_norm = j_ohmic_shape .* F .* gm1  # arbitrary normalized Joh,par = <Joh.B> / B0 (constant)

    j_non_inductive = cp1d.j_non_inductive

    I_ohmic_tor = ip - Ip_non_inductive(cp1d, eqt)
    I_tor_2_par = trapz(cp1d.grid.area, Jpar_2_Jtor(rho_tor_norm, fill(I_ohmic_tor, size(rho_tor_norm)), false, eqt)) / (I_ohmic_tor * cp1d.grid.area[end])
    I_ohmic_par_guess = I_ohmic_tor .* I_tor_2_par

    j_oh_par_norm .*= sign(I_ohmic_par_guess)
    I_ohmic_par_guess = abs(I_ohmic_par_guess)

    function cost(x)
        j_ohmic = x[1] .* j_oh_par_norm
        j_total = j_ohmic .+ j_non_inductive
        j_tor = Jpar_2_Jtor(rho_tor_norm, j_total, true, eqt)
        return abs(ip - trapz(cp1d.grid.area, j_tor)) .^ 2
    end

    res = Optim.optimize(cost, [I_ohmic_par_guess], Optim.Newton(); autodiff=:forward)
    I_ohmic_par = res.minimizer

    return I_ohmic_par .* j_oh_par_norm
end

@compat public j_ohmic_steady_state
push!(document[Symbol("Physics currents")], :j_ohmic_steady_state)

"""
    JtoR_2_JparB(rho_tor_norm::Vector{<:Real}, JtoR::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)

Given `<Jt/R>` returns `<J⋅B>`

Transformation obeys `<J⋅B> = (1/f)*(<B^2>/<1/R^2>)*(<Jt/R> + dp/dpsi*(1 - f^2*<1/R^2>/<B^2>))`

Includes_bootstrap set to true if input current includes bootstrap

NOTE: `Jtor ≂̸ JtoR`

    <Jt/R> = <Jt/R>/<1/R> * <1/R> = Jtor * <1/R> = Jtor * gm9

NOTE: `Jpar ≂̸ JparB`

    JparB = Jpar * B0
"""
function JtoR_2_JparB(rho_tor_norm::Vector{<:Real}, JtoR::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    itp_B2 = cubic_interp1d(rho_eq, eqt.profiles_1d.gm5)
    itp_invR2 = cubic_interp1d(rho_eq, eqt.profiles_1d.gm1)
    itp_f = cubic_interp1d(rho_eq, eqt.profiles_1d.f)
    itp_dpdpsi = cubic_interp1d(rho_eq, eqt.profiles_1d.dpressure_dpsi)
    return Jtor_2_JparB.(rho_tor_norm, JtoR, includes_bootstrap, Ref(itp_B2), Ref(itp_invR2), Ref(itp_f), Ref(itp_dpdpsi))
end

function Jtor_2_JparB(rho_tor_norm::Real, JtoR::Real, includes_bootstrap::Bool, B2_itp, invR2_itp, f_itp, dpdpsi_itp)
    fsa_B2 = B2_itp(rho_tor_norm)
    fsa_invR2 = invR2_itp(rho_tor_norm)
    f = f_itp(rho_tor_norm)
    dpdpsi = dpdpsi_itp(rho_tor_norm)
    if includes_bootstrap
        # add diamagnetic term to get included with bootstrap currrent
        JtoR += dpdpsi * (1.0 - fsa_invR2 * f^2 / fsa_B2) * 2pi
    end
    return fsa_B2 * JtoR / (f * fsa_invR2)
end

@compat public JtoR_2_JparB
push!(document[Symbol("Physics currents")], :JtoR_2_JparB)

"""
    JparB_2_JtoR(rho_tor_norm::Vector{<:Real}, JparB::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)

Given `<J⋅B>` returns `<Jt/R>`

Transformation obeys `<J⋅B> = (1/f)*(<B^2>/<1/R^2>)*(<Jt/R> + dp/dpsi*(1 - f^2*<1/R^2>/<B^2>))`

Includes_bootstrap set to true if input current includes bootstrap

NOTE: `Jtor ≂̸ JtoR`

    <Jt/R> = <Jt/R>/<1/R> * <1/R> = Jtor * <1/R> = Jtor * gm9

NOTE: `Jpar ≂̸ JparB`

    JparB = Jpar * B0
"""
function JparB_2_JtoR(rho_tor_norm::Vector{<:Real}, JparB::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    B2_itp = cubic_interp1d(rho_eq, eqt.profiles_1d.gm5)
    invR2_itp = cubic_interp1d(rho_eq, eqt.profiles_1d.gm1)
    f_itp = cubic_interp1d(rho_eq, eqt.profiles_1d.f)
    dpdpsi_itp = cubic_interp1d(rho_eq, eqt.profiles_1d.dpressure_dpsi)
    return JparB_2_JtoR.(rho_tor_norm, JparB, includes_bootstrap, Ref(B2_itp), Ref(invR2_itp), Ref(f_itp), Ref(dpdpsi_itp))
end

function JparB_2_JtoR(rho_tor_norm::Real, JparB::Real, includes_bootstrap::Bool, B2_itp, invR2_itp, f_itp, dpdpsi_itp)
    fsa_B2 = B2_itp(rho_tor_norm)
    fsa_invR2 = invR2_itp(rho_tor_norm)
    f = f_itp(rho_tor_norm)
    dpdpsi = dpdpsi_itp(rho_tor_norm)
    JtoR = f * fsa_invR2 * JparB / fsa_B2
    if includes_bootstrap
        # subtract diamagnetic term to get included with bootstrap currrent
        JtoR -= dpdpsi * (1.0 - fsa_invR2 * f^2 / fsa_B2) * 2pi
    end
    return JtoR
end

@compat public JparB_2_JtoR
push!(document[Symbol("Physics currents")], :JparB_2_JtoR)

"""
    Jtor_2_Jpar(rho_tor_norm::Vector{<:Real}, Jtor::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)

Given `Jtor` returns `Jpar`

Includes_bootstrap set to true if input current includes bootstrap

NOTE: `Jtor ≂̸ JtoR`

    <Jt/R> = <Jt/R>/<1/R> * <1/R> = Jtor * <1/R> = Jtor * gm9

NOTE: `Jpar ≂̸ JparB`

    JparB = Jpar * B0
"""
function Jtor_2_Jpar(rho_tor_norm::Vector{<:Real}, Jtor::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    JtoR = Jtor .* cubic_interp1d(rho_eq, eqt.profiles_1d.gm9).(rho_tor_norm)
    JparB = JtoR_2_JparB(rho_tor_norm, JtoR, includes_bootstrap, eqt)
    eq = top_ids(eqt)
    B0 = get_time_array(eq.vacuum_toroidal_field, :b0, eqt.time)
    Jpar = JparB ./ B0
    return Jpar
end

@compat public Jtor_2_Jpar
push!(document[Symbol("Physics currents")], :Jtor_2_Jpar)

"""
    Jpar_2_Jtor(rho_tor_norm::Vector{<:Real}, Jpar::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)

Given `Jpar` returns `Jtor`

Includes_bootstrap set to true if input current includes bootstrap

NOTE: `Jtor ≂̸ JtoR`

    <Jt/R> = <Jt/R>/<1/R> * <1/R> = Jtor * <1/R> = Jtor * gm9

NOTE: `Jpar ≂̸ JparB`

    JparB = Jpar * B0
"""
function Jpar_2_Jtor(rho_tor_norm::Vector{<:Real}, Jpar::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    eq = top_ids(eqt)
    B0 = get_time_array(eq.vacuum_toroidal_field, :b0, eqt.time, :constant)
    JparB = Jpar .* B0
    JtoR = JparB_2_JtoR(rho_tor_norm, JparB, includes_bootstrap, eqt)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    Jtor = JtoR ./ cubic_interp1d(rho_eq, eqt.profiles_1d.gm9).(rho_tor_norm)
    return Jtor
end

@compat public Jpar_2_Jtor
push!(document[Symbol("Physics currents")], :Jpar_2_Jtor)

"""
    vloop(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}

`Vloop = 2π * η * <J_oh⋅B> / (F * <R⁻²>)`: method emphasizes the resistive nature of the plasma
"""
function vloop(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}; method::Symbol=:area) where {T<:Real}
    @assert method in (:area, :edge, :mean)
    rho_tor_norm = cp1d.grid.rho_tor_norm
    rho_eq = eqt.profiles_1d.rho_tor_norm
    F = cubic_interp1d(rho_eq, eqt.profiles_1d.f).(rho_tor_norm)
    gm1 = cubic_interp1d(rho_eq, eqt.profiles_1d.gm1).(rho_tor_norm) # <R⁻²>
    _, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0
    Vls = 2π .* cp1d.j_ohmic .* B0 ./ (cp1d.conductivity_parallel .* F .* gm1)
    if method == :area
        return trapz(cp1d.grid.area, Vls) / cp1d.grid.area[end]
    elseif method == :edge
        return Vls[end]
    elseif method == :mean
        return sum(Vls) / length(Vls)
    end
end

"""
    vloop(eq::IMAS.equilibrium{T}, n::Int=1; time0::Float64=global_time(eq)) where {T<:Real}

`Vloop = dψ/dt`: method emphasizes the inductive nature of the loop voltage. Assumes COCOS 11.

n is the number of time_slices priors use to take the difference in psi_boundary
"""
function vloop(eq::IMAS.equilibrium{T}, n::Int=1; time0::Float64=global_time(eq)) where {T<:Real}
    @assert length(eq.time) > n "vloop from equilibrium can only be calculated in presence of at least two time slices"
    index = nearest_causal_time(eq.time, time0).index
    return (eq.time_slice[index].global_quantities.psi_boundary - eq.time_slice[index-n].global_quantities.psi_boundary) / (eq.time[index] - eq.time[index-n])
end

"""
    vloop(eq::IMAS.equilibrium{T}, n::Int; time0::Float64=global_time(eq)) where {T<:Real}

`Vloop = dψ/dt`: method emphasizes the inductive nature of the loop voltage. Assumes COCOS 11.

δt is the time difference used to take the difference in psi_boundary
"""
function vloop(eq::IMAS.equilibrium{T}, δt::Float64; time0::Float64=global_time(eq)) where {T<:Real}
    @assert δt > 0
    index0 = nearest_causal_time(eq.time, time0).index
    index1 = nearest_causal_time(eq.time, time0 - δt).index
    @assert index0 != index1 "error calculating vloop: use a larger δt"
    return (eq.time_slice[index0].global_quantities.psi_boundary - eq.time_slice[index1].global_quantities.psi_boundary) / (eq.time[index0] - eq.time[index1])
end

"""
    vloop(ct::IMAS.controllers{T}; time0::Float64=global_time(ct)) where {T<:Real}

Returns `vloop` at `time0` from controller named `ip`
"""
function vloop(ct::IMAS.controllers{T}; time0::Float64=global_time(ct)) where {T<:Real}
    vl = vloop_time(ct)
    return get_time_array(vl.time, vl.data, [time0], :linear)[1]
end

@compat public vloop
push!(document[Symbol("Physics currents")], :vloop)

"""
    vloop_time(ct::IMAS.controllers{T}) where {T<:Real}

Returns named tuple with `time` and `data` with `vloop` from controller named `ip`
"""
function vloop_time(ct::IMAS.controllers{T}) where {T<:Real}
    ctrl = controller(ct, "ip")
    return (time=ctrl.outputs.time, data=ctrl.outputs.data[1, :])
end

"""
    vloop_time(eq::IMAS.equilibrium{T}) where {T<:Real}

Returns named tuple with `time` and `data` with `vloop` from equilibrium
"""
function vloop_time(eq::IMAS.equilibrium{T}, n::Int) where {T<:Real}
    time = eq.time[n+1:end]
    data = [vloop(eq, n; time0) for time0 in time]
    return (time=time, data=data)
end

function vloop_time(eq::IMAS.equilibrium{T}, δt::Float64) where {T<:Real}
    n = findfirst(>(δt), eq.time .- eq.time[1])
    time = eq.time[n:end]
    data = [vloop(eq, δt; time0) for time0 in time]
    return (time=time, data=data)
end

"""
    vloop_time(cp::IMAS.core_profiles{T}, eq::IMAS.equilibrium) where {T<:Real}

Returns named tuple with `time` and `data` with `vloop` from core_profiles
"""
function vloop_time(cp::IMAS.core_profiles{T}, eq::IMAS.equilibrium; method::Symbol=:area) where {T<:Real}
    time = cp.time
    data = [vloop(cp.profiles_1d[time0], eq.time_slice[time0]; method) for time0 in time]
    return (time=time, data=data)
end

@compat public vloop_time
push!(document[Symbol("Physics currents")], :vloop_time)

"""
    Ip_non_inductive(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}

Integrated toroidal non-inductive current
"""
function Ip_non_inductive(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}
    return trapz(cp1d.grid.area, Jpar_2_Jtor(cp1d.grid.rho_tor_norm, cp1d.j_non_inductive, true, eqt))
end

@compat public Ip_non_inductive
push!(document[Symbol("Physics currents")], :Ip_non_inductive)

"""
    Ip_bootstrap(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}

Integrated toroidal bootstrap current
"""
function Ip_bootstrap(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}
    return trapz(cp1d.grid.area, Jpar_2_Jtor(cp1d.grid.rho_tor_norm, cp1d.j_bootstrap, true, eqt))
end

@compat public Ip_bootstrap
push!(document[Symbol("Physics currents")], :Ip_bootstrap)

"""
    Ip_ohmic(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}

Integrated toroidal ohmic current
"""
function Ip_ohmic(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}
    return trapz(cp1d.grid.area, Jpar_2_Jtor(cp1d.grid.rho_tor_norm, cp1d.j_ohmic, false, eqt))
end

@compat public Ip_ohmic
push!(document[Symbol("Physics currents")], :Ip_ohmic)

"""
    Ip(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}

Integrated toroidal total current (based on core_profiles)
"""
function Ip(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}
    return trapz(cp1d.grid.area, Jpar_2_Jtor(cp1d.grid.rho_tor_norm, cp1d.j_total, true, eqt))
end

"""
    Ip(eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}

Integrated toroidal total current (based on equilibrium)
"""
function Ip(eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}
    return Ip(eqt.profiles_1d)
end

function Ip(eqt1d::IMAS.equilibrium__time_slice___profiles_1d{T}) where {T<:Real}
    # equivalent to: trapz(eqt1d.volume, eqt1d.j_tor .* eqt1d.gm9) / (2π)
    return trapz(eqt1d.area, eqt1d.j_tor)
end

@compat public Ip
push!(document[Symbol("Physics currents")], :Ip)

"""
    plasma_lumped_resistance(dd::IMAS.dd)

Returns equivalent plasma lumped resistance in ohms
"""
function plasma_lumped_resistance(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    P_ohm = dd.core_sources.source[:ohmic].profiles_1d[].electrons.power_inside[end]
    I_ni = Ip_non_inductive(cp1d, eqt)
    I_p = Ip(cp1d, eqt)
    I_ohm = I_p - I_ni
    R_p = P_ohm / (I_p * I_ohm)
    return R_p
end

@compat public plasma_lumped_resistance
push!(document[Symbol("Physics currents")], :plasma_lumped_resistance)

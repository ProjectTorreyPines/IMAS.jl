document[Symbol("Physics currents")] = Symbol[]

function J_tor(eqt1d::IMAS.equilibrium__time_slice___profiles_1d)
    return @. (-(eqt1d.dpressure_dpsi + eqt1d.f_df_dpsi * eqt1d.gm1 / mks.μ_0) * (2π)) / eqt1d.gm9
end

"""
    j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice{T}, cp1d::IMAS.core_profiles__profiles_1d{T}, ip::T, j_ohmic_shape::AbstractVector{T}=cp1d.conductivity_parallel) where {T<:Real}

Sets `j_ohmic` parallel current density to what it would be at steady-state, based on `conductivity_parallel`, `j_non_inductive` and a target total `ip`

Requires constant loop voltage: `Vl = 2π * η * <J_oh⋅B> / (F * <R⁻²>) = constant`
"""
function j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice{T}, cp1d::IMAS.core_profiles__profiles_1d{T}, ip::T, j_ohmic_shape::AbstractVector{T}=cp1d.conductivity_parallel) where {T<:Real}
    rho_tor_norm = cp1d.grid.rho_tor_norm
    rho_eq = eqt.profiles_1d.rho_tor_norm
    F = interp1d(rho_eq, eqt.profiles_1d.f, :cubic).(rho_tor_norm)
    gm1 = interp1d(rho_eq, eqt.profiles_1d.gm1, :cubic).(rho_tor_norm) # <R⁻²>
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
    fsa_B2 = interp1d(rho_eq, eqt.profiles_1d.gm5, :cubic).(rho_tor_norm)
    fsa_invR2 = interp1d(rho_eq, eqt.profiles_1d.gm1, :cubic).(rho_tor_norm)
    f = interp1d(rho_eq, eqt.profiles_1d.f, :cubic).(rho_tor_norm)
    dpdpsi = interp1d(rho_eq, eqt.profiles_1d.dpressure_dpsi, :cubic).(rho_tor_norm)
    if includes_bootstrap
        # diamagnetic term to get included with bootstrap currrent
        JtoR_dia = dpdpsi .* (1.0 .- fsa_invR2 .* f .^ 2 ./ fsa_B2) .* 2pi
        return fsa_B2 .* (JtoR .+ JtoR_dia) ./ (f .* fsa_invR2)
    else
        return fsa_B2 .* JtoR ./ (f .* fsa_invR2)
    end
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
    fsa_B2 = interp1d(rho_eq, eqt.profiles_1d.gm5, :cubic).(rho_tor_norm)
    fsa_invR2 = interp1d(rho_eq, eqt.profiles_1d.gm1, :cubic).(rho_tor_norm)
    f = interp1d(rho_eq, eqt.profiles_1d.f, :cubic).(rho_tor_norm)
    dpdpsi = interp1d(rho_eq, eqt.profiles_1d.dpressure_dpsi, :cubic).(rho_tor_norm)
    if includes_bootstrap
        # diamagnetic term to get included with bootstrap currrent
        JtoR_dia = dpdpsi .* (1.0 .- fsa_invR2 .* f .^ 2 ./ fsa_B2) .* 2pi
        return f .* fsa_invR2 .* JparB ./ fsa_B2 .- JtoR_dia
    else
        return f .* fsa_invR2 .* JparB ./ fsa_B2
    end
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
    JtoR = Jtor .* interp1d(rho_eq, eqt.profiles_1d.gm9, :cubic).(rho_tor_norm)
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
    Jtor = JtoR ./ interp1d(rho_eq, eqt.profiles_1d.gm9, :cubic).(rho_tor_norm)
    return Jtor
end

@compat public Jpar_2_Jtor
push!(document[Symbol("Physics currents")], :Jpar_2_Jtor)

"""
    vloop(cp1d::IMAS.core_profiles__profiles_1d{T})::T where {T<:Real}

`Vloop = 2π * η * <J_oh⋅B> / (F * <R⁻²>)`: method emphasizes the resistive nature of the plasma
"""
function vloop(cp1d::IMAS.core_profiles__profiles_1d{T}, eqt::IMAS.equilibrium__time_slice{T}; method::Symbol=:area)::T where {T<:Real}
    rho_tor_norm = cp1d.grid.rho_tor_norm
    rho_eq = eqt.profiles_1d.rho_tor_norm
    F = interp1d(rho_eq, eqt.profiles_1d.f, :cubic).(rho_tor_norm)
    gm1 = interp1d(rho_eq, eqt.profiles_1d.gm1, :cubic).(rho_tor_norm) # <R⁻²>
    _, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0
    Vls = 2π .* cp1d.j_ohmic .* B0 ./ (cp1d.conductivity_parallel .* F .* gm1)
    if method === :area
        return trapz(cp1d.grid.area, Vls) / cp1d.grid.area[end]
    elseif method === :edge
        return Vls[end]
    elseif method === :mean
        return sum(Vls) / length(Vls)
    else
        throw(ArgumentError("method should be :area, :mean, or :edge"))
    end
end

"""
    vloop(eq::IMAS.equilibrium{T}; time0::Float64=global_time(eq))::T where {T<:Real}

`Vloop = dψ/dt`: method emphasizes the inductive nature of the loop voltage. Assumes COCOS 11.
"""
function vloop(eq::IMAS.equilibrium{T}; time0::Float64=global_time(eq))::T where {T<:Real}
    @assert length(eq.time) > 2 "vloop from equilibrium can only be calculated in presence of at least two time slices"
    index = nearest_causal_time(eq.time, time0).index
    return (eq.time_slice[index].global_quantities.psi_boundary - eq.time_slice[index-1].global_quantities.psi_boundary) / (eq.time[index] - eq.time[index-1])
end

"""
    vloop(ct::IMAS.controllers{T}; time0::Float64=global_time(ct))::T where {T<:Real}

Returns `vloop` at `time0` from controller named `ip`
"""
function vloop(ct::IMAS.controllers{T}; time0::Float64=global_time(ct))::T where {T<:Real}
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
    Ip(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}

Integrated toroidal total current (based on core_profiles)
"""
function Ip(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}
    return trapz(cp1d.grid.area, cp1d.j_tor)
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
    I_p = Ip(cp1d)
    I_ohm = I_p - I_ni
    R_p = P_ohm / (I_p * I_ohm)
    return R_p
end

@compat public plasma_lumped_resistance
push!(document[Symbol("Physics currents")], :plasma_lumped_resistance)

"""
    j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Sets j_ohmic to what it would be at steady-state, based on parallel conductivity and j_non_inductive
"""
function j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    j_non_inductive_tor = Jpar_2_Jtor(cp1d.grid.rho_tor_norm, cp1d.j_non_inductive, true, eqt)
    I_ohmic = eqt.global_quantities.ip - integrate(cp1d.grid.area, j_non_inductive_tor)
    return I_ohmic .* cp1d.conductivity_parallel ./ integrate(cp1d.grid.area, cp1d.conductivity_parallel)
end

"""
    j_ohmic_steady_state!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Sets j_ohmic as expression in core_profiles that evaluates to what it would be at steady-state, based on parallel conductivity and j_non_inductive
"""
function j_ohmic_steady_state!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    cp1d.j_ohmic = j_ohmic_steady_state(eqt, cp1d)
    empty!(cp1d, :j_total) # restore total as expression, to make things self-consistent
    return nothing
end

"""
    j_total_from_equilibrium!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Sets j_total as expression in core_profiles that evaluates to the total parallel current in the equilibrirum
"""
function j_total_from_equilibrium!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    cp1d.j_total = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.j_parallel, :cubic).(cp1d.grid.rho_tor_norm)
    empty!(cp1d, :j_ohmic) # restore ohmic as expression, to make things self-consistent
    return nothing
end

"""
    JtoR_2_JparB(rho_tor_norm, JtoR, includes_bootstrap, eqt::IMAS.equilibrium__time_slice)

Given <Jt/R> returns <J⋅B>
Transformation obeys <J⋅B> = (1/f)*(<B^2>/<1/R^2>)*(<Jt/R> + dp/dpsi*(1 - f^2*<1/R^2>/<B^2>))
Includes_bootstrap set to true if input current includes bootstrap
NOTE: Jtor ≂̸ JtoR
JtoR = = <Jt/R> = <Jt/R>/<1/R> * <1/R> = Jtor * <1/R> = Jtor * gm9
NOTE: Jpar ≂̸ JparB
JparB = Jpar * B0
"""
function JtoR_2_JparB(rho_tor_norm, JtoR, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    fsa_B2 = interp1d(rho_eq, eqt.profiles_1d.gm5).(rho_tor_norm)
    fsa_invR2 = interp1d(rho_eq, eqt.profiles_1d.gm1).(rho_tor_norm)
    f = interp1d(rho_eq, eqt.profiles_1d.f).(rho_tor_norm)
    dpdpsi = interp1d(rho_eq, eqt.profiles_1d.dpressure_dpsi).(rho_tor_norm)
    if includes_bootstrap
        # diamagnetic term to get included with bootstrap currrent
        JtoR_dia = dpdpsi .* (1.0 .- fsa_invR2 .* f .^ 2 ./ fsa_B2) .* 2pi
        return fsa_B2 .* (JtoR .+ JtoR_dia) ./ (f .* fsa_invR2)
    else
        return fsa_B2 * JtoR / (f * fsa_invR2)
    end
end

"""
    JparB_2_JtoR(rho_tor_norm, JparB, includes_bootstrap, eqt::IMAS.equilibrium__time_slice)

Given <J⋅B> returns <Jt/R>
Transformation obeys <J⋅B> = (1/f)*(<B^2>/<1/R^2>)*(<Jt/R> + dp/dpsi*(1 - f^2*<1/R^2>/<B^2>))
Includes_bootstrap set to true if input current includes bootstrap
NOTE: Jtor ≂̸ JtoR
JtoR = = <Jt/R> = <Jt/R>/<1/R> * <1/R> = Jtor * <1/R> = Jtor * gm9
NOTE: Jpar ≂̸ JparB
JparB = Jpar * B0
"""
function JparB_2_JtoR(rho_tor_norm, JparB, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    fsa_B2 = interp1d(rho_eq, eqt.profiles_1d.gm5).(rho_tor_norm)
    fsa_invR2 = interp1d(rho_eq, eqt.profiles_1d.gm1).(rho_tor_norm)
    f = interp1d(rho_eq, eqt.profiles_1d.f).(rho_tor_norm)
    dpdpsi = interp1d(rho_eq, eqt.profiles_1d.dpressure_dpsi).(rho_tor_norm)
    if includes_bootstrap
        # diamagnetic term to get included with bootstrap currrent
        JtoR_dia = dpdpsi .* (1.0 .- fsa_invR2 .* f .^ 2 ./ fsa_B2) .* 2pi
        return f .* fsa_invR2 .* JparB ./ fsa_B2 .- JtoR_dia
    else
        return f .* fsa_invR2 .* JparB ./ fsa_B2
    end
end

function Jpar_2_Jtor(rho_tor_norm, Jpar, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    eq = top_ids(eqt)
    B0 = interp1d(eq.time, eq.vacuum_toroidal_field.b0, :constant).(eqt.time)
    JparB = Jpar .* B0
    JtoR = JparB_2_JtoR(rho_tor_norm, JparB, includes_bootstrap, eqt)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    Jtor = JtoR ./ interp1d(rho_eq, eqt.profiles_1d.gm9).(rho_tor_norm)
    return Jtor
end

function Jtor_2_Jpar(rho_tor_norm, Jtor, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    JtoR = Jtor .* interp1d(rho_eq, eqt.profiles_1d.gm9).(rho_tor_norm)
    JparB = JtoR_2_JparB(rho_tor_norm, JtoR, includes_bootstrap, eqt)
    eq = top_ids(eqt)
    B0 = interp1d(eq.time, eq.vacuum_toroidal_field.b0, :constant).(eqt.time)
    Jpar = JparB ./ B0
    return Jpar
end
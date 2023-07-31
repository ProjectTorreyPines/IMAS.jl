"""
    j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Sets j_ohmic parallel current density to what it would be at steady-state, based on parallel conductivity and j_non_inductive and a target Ip
"""
function j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    j_oh_par_norm = cp1d.conductivity_parallel ./ integrate(cp1d.grid.area, cp1d.conductivity_parallel)

    j_non_inductive_tor = Jpar_2_Jtor(cp1d.grid.rho_tor_norm, cp1d.j_non_inductive, true, eqt)
    I_ohmic_tor = eqt.global_quantities.ip - integrate(cp1d.grid.area, j_non_inductive_tor)
    I_tor_2_par = integrate(cp1d.grid.area, Jpar_2_Jtor(cp1d.grid.rho_tor_norm, fill(I_ohmic_tor, size(cp1d.grid.rho_tor_norm)), false, eqt)) / (I_ohmic_tor * cp1d.grid.area[end])
    I_ohmic_par_guess = I_ohmic_tor .* I_tor_2_par

    j_oh_par_norm .*= sign(I_ohmic_par_guess)
    I_ohmic_par_guess = abs(I_ohmic_par_guess)

    function cost(x)
        j_ohmic = x[1] .* j_oh_par_norm
        j_total = j_ohmic .+ cp1d.j_non_inductive
        j_tor = Jpar_2_Jtor(cp1d.grid.rho_tor_norm, j_total, true, eqt)
        return abs(eqt.global_quantities.ip - IMAS.integrate(cp1d.grid.area, j_tor)) .^ 2
    end

    res = Optim.optimize(cost, [I_ohmic_par_guess], Optim.Newton(); autodiff=:forward)
    I_ohmic_par = res.minimizer

    return I_ohmic_par .* j_oh_par_norm
end

"""
    j_ohmic_steady_state!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Sets j_ohmic parallel current density to what it would be at steady-state, based on parallel conductivity and j_non_inductive and a target Ip
"""
function j_ohmic_steady_state!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    cp1d.j_ohmic = j_ohmic_steady_state(eqt, cp1d)
    # restore j_total and j_tor as expression, to make sure things are self-consistent
    empty!(cp1d, :j_total)
    empty!(cp1d, :j_tor)
    return nothing
end

"""
    j_total_from_equilibrium!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Sets j_total parallel current density as expression in core_profiles that evaluates to the total parallel current in the equilibrium
"""
function j_total_from_equilibrium!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    cp1d.j_total = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.j_parallel, :cubic).(cp1d.grid.rho_tor_norm)
    # restore j_ohmic and j_tor as expression, to make sure things are self-consistent
    empty!(cp1d, :j_ohmic)
    empty!(cp1d, :j_tor)
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

function Jpar_2_Jtor(rho_tor_norm::Vector{<:Real}, Jpar::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    eq = top_ids(eqt)
    B0 = get_time_array(eq.vacuum_toroidal_field, :b0, eqt.time, :constant)
    JparB = Jpar .* B0
    JtoR = JparB_2_JtoR(rho_tor_norm, JparB, includes_bootstrap, eqt)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    Jtor = JtoR ./ interp1d(rho_eq, eqt.profiles_1d.gm9).(rho_tor_norm)
    return Jtor
end

function Jtor_2_Jpar(rho_tor_norm::Vector{<:Real}, Jtor::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    JtoR = Jtor .* interp1d(rho_eq, eqt.profiles_1d.gm9).(rho_tor_norm)
    JparB = JtoR_2_JparB(rho_tor_norm, JtoR, includes_bootstrap, eqt)
    eq = top_ids(eqt)
    B0 = interp1d(eq.time, eq.vacuum_toroidal_field.b0, :constant).(eqt.time)
    Jpar = JparB ./ B0
    return Jpar
end
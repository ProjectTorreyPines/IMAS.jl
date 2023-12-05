"""
    j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice{T}, cp1d::IMAS.core_profiles__profiles_1d{T}, Ip::T) where {T<:Real}

Sets j_ohmic parallel current density to what it would be at steady-state, based on parallel conductivity and j_non_inductive and a target Ip
"""
function j_ohmic_steady_state(eqt::IMAS.equilibrium__time_slice{T}, cp1d::IMAS.core_profiles__profiles_1d{T}, Ip::T) where {T<:Real}
    j_oh_par_norm = cp1d.conductivity_parallel ./ integrate(cp1d.grid.area, cp1d.conductivity_parallel)

    rho_tor_norm = cp1d.grid.rho_tor_norm
    j_non_inductive = cp1d.j_non_inductive

    j_non_inductive_tor = Jpar_2_Jtor(rho_tor_norm, j_non_inductive, true, eqt)
    I_ohmic_tor = Ip - integrate(cp1d.grid.area, j_non_inductive_tor)
    I_tor_2_par = integrate(cp1d.grid.area, Jpar_2_Jtor(rho_tor_norm, fill(I_ohmic_tor, size(rho_tor_norm)), false, eqt)) / (I_ohmic_tor * cp1d.grid.area[end])
    I_ohmic_par_guess = I_ohmic_tor .* I_tor_2_par

    j_oh_par_norm .*= sign(I_ohmic_par_guess)
    I_ohmic_par_guess = abs(I_ohmic_par_guess)

    function cost(x)
        j_ohmic = x[1] .* j_oh_par_norm
        j_total = j_ohmic .+ j_non_inductive
        j_tor = Jpar_2_Jtor(rho_tor_norm, j_total, true, eqt)
        return abs(Ip - IMAS.integrate(cp1d.grid.area, j_tor)) .^ 2
    end

    res = Optim.optimize(cost, [I_ohmic_par_guess], Optim.Newton(); autodiff=:forward)
    I_ohmic_par = res.minimizer

    return I_ohmic_par .* j_oh_par_norm
end

"""
    j_ohmic_steady_state!(eqt::IMAS.equilibrium__time_slice{T}, cp1d::IMAS.core_profiles__profiles_1d{T}; Ip::T) where {T<:Real}

Sets j_ohmic parallel current density to what it would be at steady-state, based on parallel conductivity and j_non_inductive and a target Ip
"""
function j_ohmic_steady_state!(eqt::IMAS.equilibrium__time_slice{T}, cp1d::IMAS.core_profiles__profiles_1d{T}, Ip::T) where {T<:Real}
    cp1d.j_ohmic = j_ohmic_steady_state(eqt, cp1d, Ip)
    # empty j_total and j_tor to turn them into expressions and make sure things are self-consistent
    empty!(cp1d, :j_total)
    empty!(cp1d, :j_tor)
    return nothing
end

"""
    JtoR_2_JparB(rho_tor_norm::Vector{<:Real}, JtoR::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)

Given <Jt/R> returns <J⋅B>

Transformation obeys <J⋅B> = (1/f)*(<B^2>/<1/R^2>)*(<Jt/R> + dp/dpsi*(1 - f^2*<1/R^2>/<B^2>))

Includes_bootstrap set to true if input current includes bootstrap

NOTE: Jtor ≂̸ JtoR

    JtoR = <Jt/R> = <Jt/R>/<1/R> * <1/R> = Jtor * <1/R> = Jtor * gm9

NOTE: Jpar ≂̸ JparB

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
        return fsa_B2 * JtoR / (f * fsa_invR2)
    end
end

"""
    JparB_2_JtoR(rho_tor_norm::Vector{<:Real}, JparB::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)

Given <J⋅B> returns <Jt/R>

Transformation obeys <J⋅B> = (1/f)*(<B^2>/<1/R^2>)*(<Jt/R> + dp/dpsi*(1 - f^2*<1/R^2>/<B^2>))

Includes_bootstrap set to true if input current includes bootstrap

NOTE: Jtor ≂̸ JtoR

    JtoR = <Jt/R> = <Jt/R>/<1/R> * <1/R> = Jtor * <1/R> = Jtor * gm9

NOTE: Jpar ≂̸ JparB

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

function Jpar_2_Jtor(rho_tor_norm::Vector{<:Real}, Jpar::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    eq = top_ids(eqt)
    B0 = get_time_array(eq.vacuum_toroidal_field, :b0, eqt.time, :constant)
    JparB = Jpar .* B0
    JtoR = JparB_2_JtoR(rho_tor_norm, JparB, includes_bootstrap, eqt)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    Jtor = JtoR ./ interp1d(rho_eq, eqt.profiles_1d.gm9, :cubic).(rho_tor_norm)
    return Jtor
end

function Jtor_2_Jpar(rho_tor_norm::Vector{<:Real}, Jtor::Vector{<:Real}, includes_bootstrap::Bool, eqt::IMAS.equilibrium__time_slice)
    rho_eq = eqt.profiles_1d.rho_tor_norm
    JtoR = Jtor .* interp1d(rho_eq, eqt.profiles_1d.gm9, :cubic).(rho_tor_norm)
    JparB = JtoR_2_JparB(rho_tor_norm, JtoR, includes_bootstrap, eqt)
    eq = top_ids(eqt)
    B0 = get_time_array(eq.vacuum_toroidal_field, :b0, eqt.time)
    Jpar = JparB ./ B0
    return Jpar
end

"""
    vloop(cp1d::IMAS.core_profiles__profiles_1d{T})::T where {T<:Real}

Vloop = η*J: method emphasizes the resistive nature of the plasma.
"""
function vloop(cp1d::IMAS.core_profiles__profiles_1d{T})::T where {T<:Real}
    return integrate(cp1d.grid.area, cp1d.j_tor ./ cp1d.conductivity_parallel) / cp1d.grid.area[end]
end

"""
    vloop(eq::IMAS.equilibrium{T}; time0::Float64=global_time(eq))::T where {T<:Real}

`Vloop = dψ/dt` method emphasizes the inductive nature of the loop voltage.
"""
function vloop(eq::IMAS.equilibrium{T}; time0::Float64=global_time(eq))::T where {T<:Real}
    @assert length(eq.time) > 2 "vloop from equilibrium can only be calculated in presence of at least two time slices"
    index = causal_time_index(eq.time, time0)
    return (eq.time_slice[index].global_quantities.psi_boundary - eq.time_slice[index-1].global_quantities.psi_boundary) / (eq.time[index] - eq.time[index-1])
end

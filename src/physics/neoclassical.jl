document[Symbol("Physics neoclassical")] = Symbol[]

"""
    spitzer_conductivity(ne, Te, Zeff)

Calculates the Spitzer conductivity in [1/(Ω*m)]
"""
function spitzer_conductivity(ne, Te, Zeff)
    return @. 1.9012e4 * Te^1.5 / (Zeff * (0.58 + 0.74 / (0.76 + Zeff)) * lnLambda_e(ne, Te))
end

@compat public spitzer_conductivity
push!(document[Symbol("Physics neoclassical")], :spitzer_conductivity)

"""
    collision_frequencies(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}

Returns named tuple with `nue`, `nui`, `nu_exch` in [1/s]

NOTE: transpiled from the TGYRO `collision_rates` subroutine
"""
function collision_frequencies(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}
    # from TGYRO `collision_rates` subroutine
    mp = cgs.mp # g
    me = cgs.me # g
    e = cgs.e # statcoul
    k = cgs.k # erg/eV

    Te = cp1d.electrons.temperature # ev
    ne = cp1d.electrons.density_thermal / 1E6 # cm^-3

    loglam = 24.0 .- log.(sqrt.(ne) ./ Te)

    # 1/tau_ee (Belli 2008) in 1/s
    nue = @. sqrt(2) * pi * ne * e^4.0 * loglam / (sqrt(me) * (k * Te)^1.5)

    # 1/tau_ii (Belli 2008) in 1/s
    nui = zeros(T, length(Te))
    for ion in cp1d.ion
        if !ismissing(ion, :temperature) # ion temperature may be missing for purely fast-ions species
            Ti = ion.temperature
            ni = ion.density_thermal / 1E6
            Zi = avgZ(ion)
            mi = ion.element[1].a * mp
            nui += @. sqrt(2) * pi * ni * Zi * e^4.0 * loglam / (sqrt(mi) * (k * Ti)^1.5)
        end
    end

    # c_exch = 1.8e-19 is the formulary exch. coefficient
    c_exch = 2.0 * (4.0 / 3) * sqrt(2.0 * pi) * e^4 / k^1.5

    # nu_exch in 1/s
    nu_exch = zeros(T, length(Te))
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            Ti = ion.temperature
            ni = ion.density_thermal / 1E6
            Zi = avgZ(ion)
            mi = ion.element[1].a * mp
            nu_exch .+= @. c_exch * sqrt(me * mi) * Zi^2 * ni * loglam / (me * Ti + mi * Te)^1.5
        end
    end

    return (nue=nue, nui=nui, nu_exch=nu_exch)
end

@compat public collision_frequencies
push!(document[Symbol("Physics neoclassical")], :collision_frequencies)

"""
    Sauter_neo2021_bootstrap(dd::IMAS.DD; neo_2021::Bool, same_ne_ni::Bool)

Calculates bootstrap current from the data dictionary. Both keyword arguments are required.

See `Sauter_neo2021_bootstrap(eqt, cp1d; ...)` for details on the keyword arguments.
"""
function Sauter_neo2021_bootstrap(dd::IMAS.DD; neo_2021::Bool, same_ne_ni::Bool)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    return Sauter_neo2021_bootstrap(eqt, cp1d; neo_2021, same_ne_ni)
end

"""
    Sauter_neo2021_bootstrap(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; neo_2021::Bool=true, same_ne_ni::Bool=false)

Calculates bootstrap current

  - neo_2021: A.Redl, et al., Phys. Plasma 28, 022502 (2021) instead of O Sauter, et al., Phys. Plasmas 9, 5140 (2002); doi:10.1063/1.1517052 (https://crppwww.epfl.ch/~sauter/neoclassical)

  - same_ne_ni: assume same inverse scale length for electrons and ions
"""
function Sauter_neo2021_bootstrap(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; neo_2021::Bool=true, same_ne_ni::Bool=false)
    psi = cp1d.grid.psi

    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature
    Ti = cp1d.t_i_average

    pe = cp1d.electrons.pressure
    p = cp1d.pressure_thermal

    Zeff = cp1d.zeff

    rho = cp1d.grid.rho_tor_norm
    rho_eq = eqt.profiles_1d.rho_tor_norm
    fT = interp1d(rho_eq, eqt.profiles_1d.trapped_fraction).(rho)
    I_psi = interp1d(rho_eq, eqt.profiles_1d.f).(rho)

    ν_es = nuestar(eqt, cp1d; Zeff)
    ν_is = nuistar(eqt, cp1d)

    B0 = B0_geo(eqt)
    ip = eqt.global_quantities.ip

    return Sauter_neo2021_bootstrap(psi, ne, Te, Ti, pe, p, Zeff, fT, I_psi, ν_es, ν_is, ip, B0; neo_2021, same_ne_ni)
end

@compat public Sauter_neo2021_bootstrap
push!(document[Symbol("Physics neoclassical")], :Sauter_neo2021_bootstrap)

function Sauter_neo2021_bootstrap(
    psi::T,
    ne::T,
    Te::T,
    Ti::T,
    pe::T,
    p::T,
    Zeff::T,
    fT::T,
    I_psi::T,
    nuestar::T,
    nuistar::T,
    ip::Real,
    B0::Real;
    neo_2021::Bool,
    same_ne_ni::Bool
) where {T<:AbstractVector{<:Real}}
    psi = psi ./ 2π # COCOS 11 to COCOS 1 --> 1/2π
    dP_dpsi = gradient(psi, p)
    dTi_dpsi = gradient(psi, Ti)
    dTe_dpsi = gradient(psi, Te)
    dne_dpsi = gradient(psi, ne)

    R_pe = pe ./ p

    # ========== 1
    function F31(X, Zeff)
        # Equation 14a (also used in equation 16a)
        if neo_2021
            return @. (
                (1.0 + 0.15 / (Zeff^1.2 - 0.71)) * X
                -
                0.22 / (Zeff^1.2 - 0.71) * X^2
                + 0.01 / (Zeff^1.2 - 0.71) * X^3
                + 0.06 / (Zeff^1.2 - 0.71) * X^4
            )  # eq(10) from A.Redl, et al.
        else
            return @. (
                (1.0 + 1.4 / (Zeff + 1.0)) * X - 1.9 / (Zeff + 1) * X^2 + 0.3 / (Zeff + 1.0) * X^3 + 0.2 / (Zeff + 1.0) * X^4
            )
        end
        # F31 double checked against BScoeff.m, checked against neo_theory.f90
    end

    # ========== 2
    if neo_2021
        f31teff = @. fT / (
            1.0
            + (0.67 * (1.0 - 0.7 * fT) * sqrt(nuestar)) / (0.56 + 0.44 * Zeff)
            + (0.52 + 0.086 * sqrt(nuestar)) * (1.0 + 0.87 * fT) * nuestar / (1.0 + 1.13 * (Zeff - 1.0)^0.5)
        )  # eq(11) from A.Redl, et al.
    else
        f31teff = @. fT / (1.0 + (1.0 - 0.1 * fT) * sqrt(nuestar) + 0.5 * (1.0 - fT) * nuestar / Zeff)  # Equation 14b, checked
    end

    L_31 = F31(f31teff, Zeff)

    # ========== 3
    if neo_2021
        f32eiteff = @. fT / (
            1.0
            + ((0.87 * (1.0 + 0.39 * fT) * sqrt(nuestar)) / (1.0 + 2.95 * (Zeff - 1.0)^2))
            + 1.53 * (1.0 - 0.37 * fT) * nuestar * (2.0 + 0.375 * (Zeff - 1.0))
        )  # eq(16) from A.Redl, et al.
    else
        # L32 from equation 15
        f32eiteff = @. fT / (
            1.0 + (1.0 + 0.6 * fT) * sqrt(nuestar) + 0.85 * (1.0 - 0.37 * fT) * nuestar * (1.0 + Zeff)
        )
    end
    # eqn 15e double checked against BScoeff.m, checked against neo_theory.f90
    Y = f32eiteff

    # ========== 4
    if neo_2021
        F32ei = @. (
            -(0.4 + 1.93 * Zeff) / (Zeff * (0.8 + 0.6 * Zeff)) * (Y - Y^4)
            +
            5.5 / (1.5 + 2.0 * Zeff) * (Y^2 - Y^4 - 0.8 * (Y^3 - Y^4))
            -
            1.3 / (1.0 + 0.5 * Zeff) * Y^4
        )  # eq (15) from A.Redl, et al.
    else
        F32ei = @. (
            -(0.56 + 1.93 * Zeff) / (Zeff * (1.0 + 0.44 * Zeff)) * (Y - Y^4)
            +
            4.95 / (1.0 + 2.48 * Zeff) * (Y^2 - Y^4 - 0.55 * (Y^3 - Y^4))
            -
            1.2 / (1.0 + 0.5 * Zeff) * Y^4
        )  # Equation 15c # checked against neo_theory.f90
    end

    # ========== 5
    # Which Z should be used in F32ei? Neo uses the same "zeff" as elsewhere.
    if neo_2021
        f32eeteff = @. fT / (
            1.0
            + (0.23 * (1.0 - 0.96 * fT) * sqrt(nuestar)) / Zeff^0.5
            + (0.13 * (1.0 - 0.38 * fT) * nuestar / Zeff^2)
              *
              (sqrt(1.0 + 2.0 * (Zeff - 1.0)^0.5) + fT^2 * sqrt((0.075 + 0.25 * (Zeff - 1.0)^2) * nuestar))
        )  # eq(14) from A.Redl, et al.
    else
        f32eeteff = @. fT / (
            1.0 + 0.26 * (1.0 - fT) * sqrt(nuestar) + 0.18 * (1.0 - 0.37 * fT) * nuestar / sqrt(Zeff)
        )
    end
    # eqn 15d double checked against BScoeff.m, checked against neo_theory.f90

    # ========== 6
    X = f32eeteff
    if neo_2021
        F32ee = @. (
            (0.1 + 0.6 * Zeff) / (Zeff * (0.77 + 0.63 * (1.0 + (Zeff - 1.0)^1.1))) * (X - X^4)
            + 0.7 / (1.0 + 0.2 * Zeff) * (X^2 - X^4 - 1.2 * (X^3 - X^4))
            + 1.3 / (1.0 + 0.5 * Zeff) * X^4)
        # eq(13) from A.Redl, et al.
    else
        F32ee = @. (
            (0.05 + 0.62 * Zeff) / (Zeff * (1.0 + 0.44 * Zeff)) * (X - X^4)
            + 1.0 / (1.0 + 0.22 * Zeff) * (X^2 - X^4 - 1.2 * (X^3 - X^4))
            + 1.2 / (1.0 + 0.5 * Zeff) * X^4
        )  # Equation 15b, checked
    end

    L_32 = F32ee .+ F32ei

    # ========== 7
    if neo_2021
        f34teff = @. fT / (
            (1.0 + 0.25 * (1.0 - 0.7 * fT) * sqrt(nuestar) * (1.0 + 0.45 * (Zeff - 1.0)^0.5))
            +
            (0.61 * (1.0 - 0.41 * fT) * nuestar) / (Zeff^0.5)
        )  # eq(18) from from A.Redl, et al. ; which is actually f33teff
    else
        # L34 from equation 16
        # eqn 16b is very similar to 14b but there is an extra factor of 0.5 in front of the last appearance of fT
        f34teff = @. fT / (1.0 + (1.0 - 0.1 * fT) * sqrt(nuestar) + 0.5 * (1.0 - 0.5 * fT) * nuestar / Zeff)  # Equation 16b, checked
    end
    # eqn 16b double checked against BScoeff.m, checked against neo_theory.f90

    L_34 = F31(f34teff, Zeff)  # Equation 16a # double checked against BScoeff.m, checked against neo_theory.f90 (or eq(19) from from A.Redl, et al.)

    # ========== 8
    if neo_2021
        alpha0 = @. (
            -(0.62 + 0.055 * (Zeff - 1.0)) / (0.53 + 0.17 * (Zeff - 1.0)) * (1.0 - fT) / (1.0 - (0.31 - 0.065 * (Zeff - 1.0)) * fT - 0.25 * fT^2)
        )  # eq(20) from from A.Redl, et al.
    else
        # alpha from equation 17
        alpha0 = @. -1.17 * (1.0 - fT) / (1.0 - 0.22 * fT - 0.19 * fT^2)  # Checked, double checked against BScoeff.m
    end
    # Double checked against neo_theory.f90

    # ========== 9
    if neo_2021
        alpha = @. ((alpha0 + 0.7 * Zeff * fT^0.5 * sqrt(nuistar)) / (1.0 + 0.18 * sqrt(nuistar)) - 0.002 * nuistar^2 * fT^6) * (
            1.0 / (1.0 + 0.004 * nuistar^2 * fT^6)
        )  # eq (21) from A.Redl, et al.
    else
        alpha = @. ((alpha0 + 0.25 * (1.0 - fT^2) * sqrt(nuistar)) / (1.0 + 0.5 * sqrt(nuistar)) + 0.315 * nuistar^2 * fT^6) / (
            1.0 + 0.15 * nuistar^2 * fT^6
        )  # equation 17a with correction from the erratum (Sauter 2002) #checked
        # eqn 17 double checked against BScoeff.m, checked against neo_theory.f90
    end

    # Assemble the result ==========================================
    front = @. -I_psi * pe# * sign(q)
    bra1 = @. L_31 * dP_dpsi / pe  # First term in the brackets
    bra2 = @. L_32 * dTe_dpsi / Te  # Second term in the brackets
    bra3 = @. L_34 * alpha * (1.0 - R_pe) / R_pe * dTi_dpsi / Ti  # Last term in the brackets

    if same_ne_ni
        # This definition of the bootstrap current assumes that
        #  d ln(n_e)     d ln(n_i)
        #  ---------  =  ---------
        #    d psi         d psi
        # or put another way, d_psi(ln(n_e))=d_psi(ln(n_i)) same inverse scale length for electrons and ions

        # This is the second equation in Sauter"s conclusion. You are less likely to
        # get an ugly double peak if you have badly matched electron and ion density profiles, but the result doesn"t agree
        # as nicely with TRANSP. jB has more assumptions in it than j_boot, so it is intrinsically worse if you trust your
        # profile fits. If your profile fits aren"t great, it might actually improve matters by forcing some physics in.
        j_boot = @. (
            -I_psi
            * p
            * (L_31 * dne_dpsi / ne + R_pe * (L_31 + L_32) * dTe_dpsi / Te + (1 - R_pe) * (L_31 + alpha * L_34) * dTi_dpsi / Ti)
            #* sign(q)
        )
        # Double checked jB against jdotB_BS.m
    else
        # This is the second term in the first equation of the conclusion of Sauter 1999
        # with correction from Sauter 2002). It tends to be less smooth than when same_ne_ni is assumed.
        # The bootstrap current (times B) is the second term in equation for <j_par*B>,
        # the first term is ohmic current (times B)
        j_boot = @. front * (bra1 + bra2 + bra3)
    end

    return abs.(j_boot) .* sign(ip) ./ abs(B0)
end

function collisionless_bootstrap_coefficient(dd::IMAS.DD)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    return collisionless_bootstrap_coefficient(eqt, cp1d)
end

"""
    collisionless_bootstrap_coefficient(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Returns the collisional bootstrap coefficient Cbs defines as `jbootfract = Cbs * sqrt(ϵ) * βp`

See: Gi et al., Fus. Eng. Design 89 2709 (2014)

See: Wilson et al., Nucl. Fusion 32 257 (1992)
"""
function collisionless_bootstrap_coefficient(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    βp = eqt.global_quantities.beta_pol
    ϵ = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    jbootfract = Ip_bootstrap(cp1d, eqt) / eqt.global_quantities.ip
    return jbootfract / (sqrt(ϵ) * βp)
end

@compat public collisionless_bootstrap_coefficient
push!(document[Symbol("Physics neoclassical")], :collisionless_bootstrap_coefficient)

"""
    nuestar(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Calculate the electron collisionality, ν_*e, as a dimensionless measure of the
frequency of electron collisions relative to their characteristic transit frequency.
"""
function nuestar(eqt::IMAS.equilibrium__time_slice{T}, cp1d::IMAS.core_profiles__profiles_1d{T};
    Zeff::Vector{T}=cp1d.zeff) where {T<:Real}
    rho = cp1d.grid.rho_tor_norm
    Te = cp1d.electrons.temperature
    ne = cp1d.electrons.density_thermal

    R_eq = (eqt.profiles_1d.r_outboard .+ eqt.profiles_1d.r_inboard) .* 0.5
    @assert all(x -> x > 0.0, R_eq) "R_eq .> 0.0 but $R_eq"

    R_itp = linear_interp1d(eqt.profiles_1d.rho_tor_norm, R_eq)

    a_eq = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) .* 0.5
    @assert all(x -> x > 0.0, a_eq) "a_eq .> 0.0 but $a_eq"

    a_itp = linear_interp1d(eqt.profiles_1d.rho_tor_norm, a_eq)

    q_itp = linear_interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.q)

    return @. 6.921e-18 * abs(q_itp(rho)) * R_itp(rho) * ne * Zeff * lnLambda_e(ne, Te) / (Te^2 * (a_itp(rho) / R_itp(rho))^1.5)
end

@compat public nuestar
push!(document[Symbol("Physics neoclassical")], :nuestar)


"""
    nuistar(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Calculate the ion collisionality, ν_*i, as a dimensionless measure of the
frequency of ion collisions relative to their characteristic transit frequency.
"""
function nuistar(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    rho = cp1d.grid.rho_tor_norm

    R = (eqt.profiles_1d.r_outboard .+ eqt.profiles_1d.r_inboard) ./ 2.0
    R = interp1d(eqt.profiles_1d.rho_tor_norm, R).(rho)
    a = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    a = interp1d(eqt.profiles_1d.rho_tor_norm, a).(rho)

    q = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.q).(rho)
    ne = cp1d.electrons.density_thermal
    nis = hcat((ion.density_thermal for ion in cp1d.ion)...)
    ni = sum(nis; dims=2)[:, 1]
    Ti = cp1d.t_i_average

    # dominant ion (the one with the most particles)
    Zs = [ion.element[1].z_n for ion in cp1d.ion]
    Zdom = [Zs[dom[2]] for dom in argmax(nis; dims=2)[:, 1]]

    Zavg = ne ./ ni

    return @. 4.90e-18 * abs(q) * R * ni * Zdom^4 * lnLambda_i(ni, Ti, Zavg) / (Ti^2 * (a / R)^1.5)
end

@compat public nuistar
push!(document[Symbol("Physics neoclassical")], :nuistar)

"""
    lnLambda_e(ne, Te)

Compute the electron Coulomb logarithm, ln(Λ_e)
"""
function lnLambda_e(ne::Real, Te::Real)
    return 23.5 - log(sqrt(ne / 1e6) * Te^(-5.0 / 4.0)) - sqrt(1e-5 + (log(Te) - 2)^2 / 16.0)
end

"""
    lnLambda_i(ni, Ti, Zavg)

Compute the ion Coulomb logarithm, ln(Λ_i)
"""
function lnLambda_i(ni::Real, Ti::Real, Zavg::Real)
    return 30.0 - log(Zavg^3 * sqrt(ni) / (Ti^1.5))
end

"""
    neo_conductivity(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Calculates the neo-classical conductivity in 1/(Ohm*meter) based on the NEO 2021 modifcation
"""
function neo_conductivity(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    rho = cp1d.grid.rho_tor_norm
    Te = cp1d.electrons.temperature
    ne = cp1d.electrons.density_thermal
    Zeff = cp1d.zeff

    tf_itp = linear_interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.trapped_fraction)

    nue = nuestar(eqt, cp1d; Zeff)

    # neo 2021
    f33teff =
        @. tf_itp(rho) / (
            1 + 0.25 * (1 - 0.7 * tf_itp(rho)) * sqrt(nue) * (1 + 0.45 * (Zeff - 1)^0.5) +
            0.61 * (1 - 0.41 * tf_itp(rho)) * nue / Zeff^0.5
        )

    F33 = @. 1 - (1 + 0.21 / Zeff) * f33teff + 0.54 / Zeff * f33teff^2 - 0.33 / Zeff * f33teff^3

    conductivity_parallel = spitzer_conductivity(ne, Te, Zeff) .* F33

    return conductivity_parallel
end

@compat public neo_conductivity
push!(document[Symbol("Physics neoclassical")], :neo_conductivity)

"""
    sawtooth_conductivity(conductivity::AbstractVector{T}, q::AbstractVector{T}; q_sawtooth::Float64=1.0) where {T<:Real}

Jardin's model for stationary sawteeth to raise q>1
"""
function sawtooth_conductivity(conductivity::AbstractVector{T}, q::AbstractVector{T}; q_sawtooth::Float64=1.0) where {T<:Real}
    i_q = findlast(abs.(q) .< q_sawtooth)
    if i_q !== nothing
        conductivity[1:i_q] .= conductivity[i_q]
    end
    return conductivity
end

@compat public sawtooth_conductivity
push!(document[Symbol("Physics neoclassical")], :sawtooth_conductivity)

#= ================================================================================== =#
#  Species rotation from radial force balance
#
#  Every species s on a flux surface has a flow of the form (Hinton & Wong 1985)
#      V_s = ω_s(ψ) R² ∇φ + u_s(ψ) B ,     ω_s = ω_E - ω_*s
#  with the common ExB frequency ω_E = -dΦ/dψ (`rotation_frequency_tor_sonic`) and the
#  diamagnetic frequency ω_*s = (1 / (Z_s e n_s)) dp_s/dψ, where ψ is the poloidal flux
#  per radian. The flux-surface toroidal rotation frequency of species s is then
#      Ω_φ,s = ω_E + Δ_s ,     Δ_s = -ω_*s + u_s F ⟨R⁻²⟩
#  so all species are tied to one free rotation function per surface, and Δ_s depends only
#  on the n, T profiles (and on the poloidal-flow model that sets u_s).
#= ================================================================================== =#

"""
    poloidal_flux_per_radian(cp1d::IMAS.core_profiles__profiles_1d)

Poloidal flux on the `cp1d` grid normalized as in Hinton & Wong (ψ per radian, B_p = ∇φ × ∇ψ),
the normalization in which `ω_E = -dΦ/dψ` is an angular frequency.

IMAS stores the total poloidal flux in Wb (COCOS 11, `exp_Bp = 1`), so that |∇ψ| = 2π R B_p.
"""
function poloidal_flux_per_radian(cp1d::IMAS.core_profiles__profiles_1d)
    COCOS = cocos(internal_cocos)
    return cp1d.grid.psi .* (COCOS.sigma_Bp / (2π)^COCOS.exp_Bp)
end

@compat public poloidal_flux_per_radian
push!(document[Symbol("Physics neoclassical")], :poloidal_flux_per_radian)

"""
    rotation_taper(rho::Real, rho_correction::Real)

Factor `min(ρ/ρ_c, 1)²` that tapers gradient-driven rotation offsets to zero near the magnetic
axis, where gradients of fast-ion-subtracted thermal pressure are noisy. `rho_correction = 0` disables it.
"""
function rotation_taper(rho::Real, rho_correction::Real)
    return rho_correction > 0 ? min(rho / rho_correction, one(rho))^2 : one(rho)
end

"""
    diamagnetic_rotation_frequency(cp1d::IMAS.core_profiles__profiles_1d, ion::IMAS.core_profiles__profiles_1d___ion; rho_correction::Real=0.2, psi=poloidal_flux_per_radian(cp1d))

Diamagnetic angular frequency of an ion species `ω_*s = (1 / (Z_s e n_s)) dp_s/dψ` [rad/s], with ψ per radian.

Returns zero where the thermal density of the species vanishes.
"""
function diamagnetic_rotation_frequency(
    cp1d::IMAS.core_profiles__profiles_1d,
    ion::IMAS.core_profiles__profiles_1d___ion;
    rho_correction::Real=0.2,
    psi::AbstractVector{<:Real}=poloidal_flux_per_radian(cp1d))
    rho = cp1d.grid.rho_tor_norm
    ni = ion.density_thermal
    dp_dpsi = gradient(psi, ion.density_thermal .* ion.temperature .* mks.e)
    return @. ifelse(ni > 0, dp_dpsi / (max(ni, eps(one(ni))) * ion.z_ion * mks.e), zero(dp_dpsi)) * rotation_taper(rho, rho_correction)
end

@compat public diamagnetic_rotation_frequency
push!(document[Symbol("Physics neoclassical")], :diamagnetic_rotation_frequency)

"""
    thermal_rotation_frequency(cp1d::IMAS.core_profiles__profiles_1d, ion::IMAS.core_profiles__profiles_1d___ion; rho_correction::Real=0.2, psi=poloidal_flux_per_radian(cp1d))

Temperature-gradient angular frequency of an ion species `ω_T,s = (1 / (Z_s e)) dT_s/dψ` [rad/s], with ψ per radian.
"""
function thermal_rotation_frequency(
    cp1d::IMAS.core_profiles__profiles_1d,
    ion::IMAS.core_profiles__profiles_1d___ion;
    rho_correction::Real=0.2,
    psi::AbstractVector{<:Real}=poloidal_flux_per_radian(cp1d))
    rho = cp1d.grid.rho_tor_norm
    # temperature is in eV, so (1/(Z e)) d(e T)/dψ = (dT/dψ) / Z
    return gradient(psi, ion.temperature) ./ ion.z_ion .* rotation_taper.(rho, rho_correction)
end

@compat public thermal_rotation_frequency
push!(document[Symbol("Physics neoclassical")], :thermal_rotation_frequency)

"""
    main_ion_index(cp1d::IMAS.core_profiles__profiles_1d)

Index of the thermal ion species with the largest (grid-summed) thermal density
"""
function main_ion_index(cp1d::IMAS.core_profiles__profiles_1d)
    return argmax([hasdata(ion, :density_thermal) || hasdata(ion, :density) ? sum(ion.density_thermal) : 0.0 for ion in cp1d.ion])
end

@compat public main_ion_index
push!(document[Symbol("Physics neoclassical")], :main_ion_index)

"""
    poloidal_flow_coefficient(nustar::Real, epsilon::Real)

Neoclassical main-ion poloidal flow coefficient `k_i`, interpolated across collisionality regimes

    k_i = [(1.17 - 0.35 √ν*) / (1 + 0.7 √ν*) - 2.1 ν*² ε³] / (1 + ν*² ε³)

which tends to 1.17 in the banana regime and -2.1 in the Pfirsch-Schlüter regime
(Hinton & Hazeltine 1976; interpolation as in Kim, Diamond & Groebner, Phys. Fluids B 3, 2050 (1991)).
"""
function poloidal_flow_coefficient(nustar::Real, epsilon::Real)
    sn = sqrt(nustar)
    ne3 = nustar^2 * epsilon^3
    return ((1.17 - 0.35 * sn) / (1.0 + 0.7 * sn) - 2.1 * ne3) / (1.0 + ne3)
end

@compat public poloidal_flow_coefficient
push!(document[Symbol("Physics neoclassical")], :poloidal_flow_coefficient)

"""
    rotation_frequency_offset(
        cp1d::IMAS.core_profiles__profiles_1d,
        ion::IMAS.core_profiles__profiles_1d___ion;
        poloidal_flow::Symbol=:none,
        eqt::Union{Nothing,IMAS.equilibrium__time_slice}=nothing,
        rho_correction::Real=0.2,
        psi=poloidal_flux_per_radian(cp1d))

Offset `Δ_s` of the toroidal rotation frequency of species `ion` from the ExB frequency, `Ω_φ,s = ω_E + Δ_s` [rad/s]

`poloidal_flow` selects the model for the poloidal flow of each species:

  - `:none`: no poloidal flow for any species, `Δ_s = -ω_*s`

  - `:analytic` (requires `eqt`): main ion (see `main_ion_index`) with the neoclassical poloidal flow
    `Ω_φ,i = ω_E - ω_*i + k_i f ω_T,i`, where `k_i = poloidal_flow_coefficient(ν*_i, ε)` and `f = F² ⟨R⁻²⟩ / ⟨B²⟩`.
    Every other ion has its parallel flow locked to the main ion by friction (trace-impurity limit),
    `Ω_φ,z = Ω_φ,i - (1 - f)(ω_*z - ω_*i)`. NOTE: the impurity thermal force, which is of the same order
    as the `k_i` term, is neglected: for quantitative impurity rotation use species flows from NEO.
"""
function rotation_frequency_offset(
    cp1d::IMAS.core_profiles__profiles_1d,
    ion::IMAS.core_profiles__profiles_1d___ion;
    poloidal_flow::Symbol=:none,
    eqt::Union{Nothing,IMAS.equilibrium__time_slice}=nothing,
    rho_correction::Real=0.2,
    psi::AbstractVector{<:Real}=poloidal_flux_per_radian(cp1d))

    if poloidal_flow == :none
        return -diamagnetic_rotation_frequency(cp1d, ion; rho_correction, psi)

    elseif poloidal_flow == :analytic
        @assert eqt !== nothing "rotation_frequency_offset: `poloidal_flow=:analytic` requires the equilibrium time slice `eqt`"
        rho = cp1d.grid.rho_tor_norm
        eqt1d = eqt.profiles_1d
        itp(y) = interp1d(eqt1d.rho_tor_norm, y).(rho)
        f_geo = itp(eqt1d.f .^ 2 .* eqt1d.gm1 ./ eqt1d.gm5)
        epsilon = itp((eqt1d.r_outboard .- eqt1d.r_inboard) ./ (eqt1d.r_outboard .+ eqt1d.r_inboard))
        k_i = poloidal_flow_coefficient.(nuistar(eqt, cp1d), epsilon)

        main_ion = cp1d.ion[main_ion_index(cp1d)]
        ωstar_main = diamagnetic_rotation_frequency(cp1d, main_ion; rho_correction, psi)
        Δ_main = -ωstar_main .+ k_i .* f_geo .* thermal_rotation_frequency(cp1d, main_ion; rho_correction, psi)
        if ion === main_ion
            return Δ_main
        else
            ωstar_ion = diamagnetic_rotation_frequency(cp1d, ion; rho_correction, psi)
            return Δ_main .- (1.0 .- f_geo) .* (ωstar_ion .- ωstar_main)
        end

    else
        error("rotation_frequency_offset: `poloidal_flow` must be `:none` or `:analytic`, not `:$poloidal_flow`")
    end
end

@compat public rotation_frequency_offset
push!(document[Symbol("Physics neoclassical")], :rotation_frequency_offset)

"""
    ωtor2sonic(cp1d::IMAS.core_profiles__profiles_1d; ind::Int=0, poloidal_flow::Symbol=:none, eqt=nothing, rho_correction::Real=0.2)

ExB rotation frequency `ω_E = Ω_φ,ref - Δ_ref` [rad/s] from the toroidal rotation of a reference ion species

`ind` is the index of the reference species (0 = first ion with non-zero `rotation_frequency_tor`).
See `rotation_frequency_offset` for `poloidal_flow` and `rho_correction`.
"""
function ωtor2sonic(
    cp1d::IMAS.core_profiles__profiles_1d{T};
    ind::Int=0,
    poloidal_flow::Symbol=:none,
    eqt::Union{Nothing,IMAS.equilibrium__time_slice}=nothing,
    rho_correction::Real=0.2) where {T<:Real}
    if ind == 0
        ind = findfirst(x -> hasdata(x, :rotation_frequency_tor) && sum(abs.(x.rotation_frequency_tor)) != 0.0, cp1d.ion)
        if ind === nothing
            return zero(cp1d.grid.rho_tor_norm)
        end
    end
    ion = cp1d.ion[ind]
    return ion.rotation_frequency_tor .- rotation_frequency_offset(cp1d, ion; poloidal_flow, eqt, rho_correction)
end

@compat public ωtor2sonic
push!(document[Symbol("Physics neoclassical")], :ωtor2sonic)

"""
    ωtor2sonic!(cp1d::IMAS.core_profiles__profiles_1d; ind::Int=0, poloidal_flow::Symbol=:none, eqt=nothing, rho_correction::Real=0.2)

Updates `cp1d.rotation_frequency_tor_sonic` in place from the toroidal rotation of a reference ion species (see `ωtor2sonic`)
"""
function ωtor2sonic!(cp1d::IMAS.core_profiles__profiles_1d{T}; kw...) where {T<:Real}
    return cp1d.rotation_frequency_tor_sonic = ωtor2sonic(cp1d; kw...)
end

@compat public ωtor2sonic!
push!(document[Symbol("Physics neoclassical")], :ωtor2sonic!)

"""
    sonic2ωtor(cp1d::IMAS.core_profiles__profiles_1d, ion::IMAS.core_profiles__profiles_1d___ion; poloidal_flow::Symbol=:none, eqt=nothing, rho_correction::Real=0.2)

Toroidal rotation frequency of ion species `ion`, `Ω_φ,s = ω_E + Δ_s` [rad/s], from the ExB rotation `rotation_frequency_tor_sonic`

See `rotation_frequency_offset` for `poloidal_flow` and `rho_correction`.
"""
function sonic2ωtor(
    cp1d::IMAS.core_profiles__profiles_1d,
    ion::IMAS.core_profiles__profiles_1d___ion;
    poloidal_flow::Symbol=:none,
    eqt::Union{Nothing,IMAS.equilibrium__time_slice}=nothing,
    rho_correction::Real=0.2)
    return cp1d.rotation_frequency_tor_sonic .+ rotation_frequency_offset(cp1d, ion; poloidal_flow, eqt, rho_correction)
end

@compat public sonic2ωtor
push!(document[Symbol("Physics neoclassical")], :sonic2ωtor)

"""
    sonic2ωtor!(cp1d::IMAS.core_profiles__profiles_1d, ion::IMAS.core_profiles__profiles_1d___ion; kw...)

In-place conversion of ExB rotation to the toroidal rotation of `ion` (see `sonic2ωtor`)
"""
function sonic2ωtor!(cp1d::IMAS.core_profiles__profiles_1d, ion::IMAS.core_profiles__profiles_1d___ion; kw...)
    return ion.rotation_frequency_tor = sonic2ωtor(cp1d, ion; kw...)
end

"""
    sonic2ωtor!(cp1d::IMAS.core_profiles__profiles_1d; kw...)

All ions: Updates the toroidal rotation frequency of every ion species from the ExB rotation (see `sonic2ωtor`)
"""
function sonic2ωtor!(cp1d::IMAS.core_profiles__profiles_1d; kw...)
    for ion in cp1d.ion
        sonic2ωtor!(cp1d, ion; kw...)
    end
end

@compat public sonic2ωtor!
push!(document[Symbol("Physics neoclassical")], :sonic2ωtor!)

"""
    toroidal_angular_momentum_density(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice; psi=poloidal_flux_per_radian(cp1d), rho_correction::Real=0.2)

Toroidal angular momentum density of the thermal ions `Σ_s m_s n_s ⟨R²⟩ Ω_φ,s` [kg m⁻¹ s⁻¹],
with `Ω_φ,s = ω_E + Δ_s` from `rotation_frequency_offset` (no poloidal flow).

`psi` can be passed for profiles that do not carry their own `grid.psi` (e.g. a copy at a previous time step).
Electron and fast-ion contributions are neglected.
"""
function toroidal_angular_momentum_density(
    cp1d::IMAS.core_profiles__profiles_1d,
    eqt::IMAS.equilibrium__time_slice;
    psi::AbstractVector{<:Real}=poloidal_flux_per_radian(cp1d),
    rho_correction::Real=0.2)
    eqt1d = eqt.profiles_1d
    # ⟨R²⟩ (gm10) is filled by `flux_surfaces`; equilibria loaded from file may only carry ⟨R⟩ (gm8): ⟨R⟩² ≈ ⟨R²⟩ to O(ε²/2)
    R2 = hasdata(eqt1d, :gm10) ? eqt1d.gm10 : eqt1d.gm8 .^ 2
    R2_avg = interp1d(eqt1d.rho_tor_norm, R2).(cp1d.grid.rho_tor_norm)
    ωE = cp1d.rotation_frequency_tor_sonic
    L = zero(ωE)
    for ion in cp1d.ion
        hasdata(ion, :density_thermal) || hasdata(ion, :density) || continue
        Ω = ωE .+ rotation_frequency_offset(cp1d, ion; rho_correction, psi)
        L = L .+ (ion.element[1].a * mks.m_p) .* ion.density_thermal .* Ω
    end
    return L .* R2_avg
end

@compat public toroidal_angular_momentum_density
push!(document[Symbol("Physics neoclassical")], :toroidal_angular_momentum_density)

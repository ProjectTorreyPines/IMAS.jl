function spitzer_conductivity(ne, Te, Zeff)
    return 1.9012e4 .* Te .^ 1.5 ./ (Zeff .* 0.58 .+ 0.74 ./ (0.76 .+ Zeff) .* lnLambda_e(ne, Te))
end

function collision_frequencies(dd::IMAS.dd)
    # from TGYRO `collision_rates` subroutine
    cp1d = dd.core_profiles.profiles_1d[]

    Te = cp1d.electrons.temperature # ev
    ne = cp1d.electrons.density_thermal / 1E6 # cm^-3
    me = constants.m_e * 1E3 # g
    mp = constants.m_p * 1E3 # g
    e = gacode_units.e # statcoul
    k = gacode_units.k # erg/eV

    loglam = 24.0 .- log.(sqrt.(ne) ./ Te)

    # 1/tau_ee (Belli 2008) in 1/s
    nue = @. sqrt(2) * pi * ne * e^4.0 * loglam / (sqrt(me) * (k * Te)^1.5)

    # 1/tau_ii (Belli 2008) in 1/s
    nui = zeros(length(Te))
    for ion in cp1d.ion
        if !ismissing(ion, :temperature) # ion temperature may be missing for purely fast-ions species
            Ti = ion.temperature
            ni = ion.density / 1E6
            Zi = avgZ(ion.element[1].z_n, Ti)
            mi = ion.element[1].a * mp
            nui += @. sqrt(2) * pi * ni * Zi * e^4.0 * loglam / (sqrt(mi) * (k * Ti)^1.5)
        end
    end

    # c_exch = 1.8e-19 is the formulary exch. coefficient
    c_exch = 2.0 * (4.0 / 3) * sqrt(2.0 * pi) * e^4 / k^1.5

    # nu_exch in 1/s
    nu_exch = zeros(length(Te))
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            Ti = ion.temperature
            ni = ion.density / 1E6
            Zi = avgZ(ion.element[1].z_n, Ti)
            mi = ion.element[1].a * mp
            nu_exch .+= @. c_exch * sqrt(me * mi) * Zi^2 * ni * loglam / (me * Ti + mi * Te)^1.5
        end
    end

    return (nue=nue, nui=nui, nu_exch=nu_exch)
end

function Sauter_neo2021_bootstrap(dd::IMAS.dd; neo_2021::Bool=true, same_ne_ni::Bool=false)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    return Sauter_neo2021_bootstrap(eqt, cp1d; neo_2021, same_ne_ni)
end

"""
    Sauter_neo2021_bootstrap(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; neo_2021::Bool=true, same_ne_ni::Bool=false)

* neo_2021: A.Redl, et al., Phys. Plasma 28, 022502 (2021) instead of O Sauter, et al., Phys. Plasmas 9, 5140 (2002); doi:10.1063/1.1517052 (https://crppwww.epfl.ch/~sauter/neoclassical)

* same_ne_ni: assume same inverse scale length for electrons and ions
"""
function Sauter_neo2021_bootstrap(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; neo_2021::Bool=false, same_ne_ni::Bool=false)
    psi = cp1d.grid.psi

    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature
    Ti = cp1d.ion[1].temperature

    pe = cp1d.electrons.pressure
    p = cp1d.pressure_thermal

    Zeff = cp1d.zeff

    rho = cp1d.grid.rho_tor_norm
    rho_eq = eqt.profiles_1d.rho_tor_norm
    fT = interp1d(rho_eq, eqt.profiles_1d.trapped_fraction).(rho)
    I_psi = interp1d(rho_eq, eqt.profiles_1d.f).(rho)

    nuestar = IMAS.nuestar(eqt, cp1d)
    nuistar = IMAS.nuistar(eqt, cp1d)

    B0 = B0_geo(eqt)
    ip = eqt.global_quantities.ip

    return Sauter_neo2021_bootstrap(psi, ne, Te, Ti, pe, p, Zeff, fT, I_psi, nuestar, nuistar, ip, B0; neo_2021, same_ne_ni)
end

function Sauter_neo2021_bootstrap(psi::T, ne::T, Te::T, Ti::T, pe::T, p::T, Zeff::T, fT::T, I_psi::T, nuestar::T, nuistar::T, ip::Real, B0::Real; neo_2021::Bool=false, same_ne_ni::Bool=false) where {T<:AbstractVector{<:Real}}
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

function collisionless_bootstrap_coefficient(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    collisionless_bootstrap_coefficient(eqt, cp1d)
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
    jbootfract = IMAS.integrate(cp1d.grid.area, cp1d.j_bootstrap) / eqt.global_quantities.ip
    jbootfract / (sqrt(ϵ) * βp)
end

function nuestar(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    rho = cp1d.grid.rho_tor_norm
    Te = cp1d.electrons.temperature
    ne = cp1d.electrons.density
    Zeff = cp1d.zeff

    R = (eqt.profiles_1d.r_outboard .+ eqt.profiles_1d.r_inboard) ./ 2.0
    R = interp1d(eqt.profiles_1d.rho_tor_norm, R).(rho)
    a = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    a = interp1d(eqt.profiles_1d.rho_tor_norm, a).(rho)

    eps = a ./ R

    q = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.q).(rho)

    return @. 6.921e-18 * abs(q) * R * ne * Zeff * lnLambda_e(ne, Te) / (Te^2 * eps^1.5)
end

function nuistar(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    rho = cp1d.grid.rho_tor_norm

    R = (eqt.profiles_1d.r_outboard .+ eqt.profiles_1d.r_inboard) ./ 2.0
    R = interp1d(eqt.profiles_1d.rho_tor_norm, R).(rho)
    a = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    a = interp1d(eqt.profiles_1d.rho_tor_norm, a).(rho)

    eps = a ./ R

    q = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.q).(rho)
    ne = cp1d.electrons.density
    nis = hcat((ion.density for ion in cp1d.ion)...)
    ni = sum(nis; dims=2)[:, 1]
    Ti = cp1d.ion[1].temperature

    # dominant ion (the one with the most particles)
    Zs = [ion.z_ion for ion in cp1d.ion]
    Zdom = [Zs[dom[2]] for dom in argmax(nis; dims=2)[:, 1]]

    Zavg = ne ./ ni

    return @. 4.90e-18 * abs(q) * R * ni * Zdom^4 * lnLambda_i(ni, Ti, Zavg) / (Ti^2 * eps^1.5)
end

function lnLambda_e(ne, Te)
    return @. 23.5 - log(sqrt(ne / 1e6) * Te^(-5.0 / 4.0)) - sqrt(1e-5 + (log(Te) - 2)^2 / 16.0)
end

function lnLambda_i(ni, Ti, Zavg)
    return @. 30.0 - log(Zavg^3 * sqrt(ni) / (Ti^1.5))
end

"""
    neo_conductivity(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Calculates the neo-classical conductivity in 1/(Ohm*meter) based on the NEO 2021 modifcation and stores it in dd

More info see omfit_classes.utils_fusion.py neo_conductivity function
"""
function neo_conductivity(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    rho = cp1d.grid.rho_tor_norm
    Te = cp1d.electrons.temperature
    ne = cp1d.electrons.density
    Zeff = cp1d.zeff

    trapped_fraction = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.trapped_fraction).(rho)

    nue = nuestar(eqt, cp1d)

    # neo 2021
    f33teff =
        trapped_fraction ./ (
            1 .+ 0.25 .* (1 .- 0.7 .* trapped_fraction) .* sqrt.(nue) .* (1 .+ 0.45 .* (Zeff .- 1) .^ 0.5) .+
            0.61 .* (1 .- 0.41 .* trapped_fraction) .* nue ./ Zeff .^ 0.5
        )

    F33 = 1 .- (1 .+ 0.21 ./ Zeff) .* f33teff .+ 0.54 ./ Zeff .* f33teff .^ 2 .- 0.33 ./ Zeff .* f33teff .^ 3

    conductivity_parallel = spitzer_conductivity(ne, Te, Zeff) .* F33

    return conductivity_parallel
end

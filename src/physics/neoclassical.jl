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
    nue = @. sqrt(2) * pi * ne * e^4.0 * loglam / (sqrt(me) * (k * Te) ^ 1.5)

    # 1/tau_ii (Belli 2008) in 1/s
    nui = zeros(length(Te))
    for ion in cp1d.ion
        Ti = ion.temperature
        ni = ion.density / 1E6
        Zi = avgZ(ion.element[1].z_n, Ti)
        mi = ion.element[1].a * mp
        nui += @. sqrt(2) * pi * ni * Zi * e^4.0 * loglam / (sqrt(mi) * (k * Ti) ^ 1.5)
    end

    # c_exch = 1.8e-19 is the formulary exch. coefficient
    c_exch = 2.0 * (4.0 / 3) * sqrt(2.0 * pi) * e^4 / k^1.5

    # nu_exch in 1/s
    nu_exch = zeros(length(Te))
    for ion in cp1d.ion
        Ti = ion.temperature
        ni = ion.density / 1E6
        Zi = avgZ(ion.element[1].z_n, Ti)
        mi = ion.element[1].a * mp
        nu_exch .+= @. c_exch * sqrt(me * mi) * Zi^2 * ni * loglam / (me * Ti + mi * Te) ^ 1.5
    end

    return nue, nui, nu_exch
end

function Sauter_neo2021_bootstrap(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    return Sauter_neo2021_bootstrap(eqt, cp1d)
end

function Sauter_neo2021_bootstrap(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    rho = cp1d.grid.rho_tor_norm
    rho_eq = eqt.profiles_1d.rho_tor_norm

    Te = cp1d.electrons.temperature
    Ti = cp1d.ion[1].temperature
    pressure_thermal = cp1d.pressure_thermal
    R_pe = cp1d.electrons.pressure ./ pressure_thermal
    Zeff = cp1d.zeff

    psi_cp = cp1d.grid.psi ./ 2pi
    dP_dpsi = gradient(psi_cp, pressure_thermal)
    dTi_dpsi = gradient(psi_cp, Ti)
    dTe_dpsi = gradient(psi_cp, Te)

    fT = interp1d(rho_eq, eqt.profiles_1d.trapped_fraction).(rho)
    I_psi = interp1d(rho_eq, eqt.profiles_1d.f).(rho)

    nue = nuestar(eqt, cp1d)
    nui = nuistar(eqt, cp1d)

    # neo 2021
    f31teff =
        fT ./ (
            1 .+ (0.67 .* (1 .- 0.7 .* fT) .* sqrt.(nue)) ./ (0.56 .+ 0.44 .* Zeff) .+
            (0.52 .+ 0.086 .* sqrt.(nue)) .* (1 .+ 0.87 .* fT) .* nue ./ (1 .+ 1.13 .* (Zeff .- 1) .^ 0.5)
        )
    X = f31teff
    F31 = (
        (1 .+ 0.15 ./ (Zeff .^ 1.2 .- 0.71)) .* X .- 0.22 ./ (Zeff .^ 1.2 .- 0.71) .* X .^ 2 .+ 0.01 ./ (Zeff .^ 1.2 .- 0.71) .* X .^ 3 .+
        0.06 ./ (Zeff .^ 1.2 .- 0.71) .* X .^ 4
    )

    f32eeteff =
        fT ./ (
            1 .+ (0.23 .* (1 .- 0.96 .* fT) .* sqrt.(nue)) ./ Zeff .^ 0.5 .+
            (0.13 .* (1 .- 0.38 .* fT) .* nue ./ Zeff .^ 2) .*
            (sqrt.(1 .+ 2 .* (Zeff .- 1) .^ 0.5) .+ fT .^ 2 .* sqrt.((0.075 .+ 0.25 .* (Zeff .- 1) .^ 2) .* nue))
        )

    X = f32eeteff
    F32ee = (0.1 .+ 0.6 .* Zeff) ./ (Zeff .* (0.77 .+ 0.63 .* (1 .+ (Zeff .- 1) .^ 1.1))) .* (X .- X .^ 4)
    (.+0.7 ./ (1 .+ 0.2 .* Zeff) .* (X .^ 2 .- X .^ 4 .- 1.2 .* (X .^ 3 .- X .^ 4)) .+ 1.3 ./ (1 .+ 0.5 .* Zeff) .* X .^ 4)

    f32eiteff =
        fT ./
        (1 .+ ((0.87 .* (1 .+ 0.39 .* fT) .* sqrt.(nue)) ./ (1 .+ 2.95 .* (Zeff .- 1) .^ 2)) .+ 1.53 .* (1 .- 0.37 .* fT) .* nue .* (2 .+ 0.375 .* (Zeff .- 1)))

    Y = f32eiteff

    F32ei = (
        .-(0.4 .+ 1.93 .* Zeff) ./ (Zeff .* (0.8 .+ 0.6 .* Zeff)) .* (Y .- Y .^ 4) .+
        5.5 ./ (1.5 .+ 2 .* Zeff) .* (Y .^ 2 .- Y .^ 4 .- 0.8 .* (Y .^ 3 .- Y .^ 4)) .- 1.3 ./ (1 .+ 0.5 .* Zeff) .* Y .^ 4
    )

    L_32 = F32ee .+ F32ei

    f34teff = fT ./ ((1 .+ 0.25 .* (1 .- 0.7 .* fT) .* sqrt.(nue) .* (1 .+ 0.45 .* (Zeff .- 1) .^ 0.5)) .+ (0.61 .* (1 .- 0.41 .* fT) .* nue) ./ (Zeff .^ 0.5))

    X = f34teff
    L_34 = (
        (1 .+ 0.15 ./ (Zeff .^ 1.2 .- 0.71)) .* X .- 0.22 ./ (Zeff .^ 1.2 .- 0.71) .* X .^ 2 .+ 0.01 ./ (Zeff .^ 1.2 .- 0.71) .* X .^ 3 .+
        0.06 ./ (Zeff .^ 1.2 .- 0.71) .* X .^ 4
    )
    alpha0 = (.-(0.62 .+ 0.055 .* (Zeff .- 1)) ./ (0.53 .+ 0.17 .* (Zeff .- 1)) .* (1 .- fT) ./ (1 .- (0.31 .- 0.065 .* (Zeff .- 1)) .* fT .- 0.25 .* fT .^ 2))
    alpha =
        ((alpha0 .+ 0.7 .* Zeff .* fT .^ 0.5 .* sqrt.(nui)) ./ (1 .+ 0.18 .* sqrt.(nui)) .- 0.002 .* nui .^ 2 .* fT .^ 6) .*
        (1 ./ (1 .+ 0.004 .* nui .^ 2 .* fT .^ 6))

    bra1 = F31 .* dP_dpsi ./ cp1d.electrons.pressure
    bra2 = L_32 .* dTe_dpsi ./ Te
    bra3 = L_34 .* alpha .* (1 .- R_pe) ./ R_pe .* dTi_dpsi ./ Ti

    equilibrium = top_ids(eqt)
    B0 = get_time_array(equilibrium.vacuum_toroidal_field, :b0, eqt.time)
    j_boot = -I_psi .* cp1d.electrons.pressure .* (bra1 .+ bra2 .+ bra3) ./ B0

    j_boot = abs.(j_boot) .* sign(eqt.global_quantities.ip)
    return j_boot
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

    R = (eqt.profiles_1d.r_outboard + eqt.profiles_1d.r_inboard) / 2.0
    R = interp1d(eqt.profiles_1d.rho_tor_norm, R).(rho)
    a = (eqt.profiles_1d.r_outboard - eqt.profiles_1d.r_inboard) / 2.0
    a = interp1d(eqt.profiles_1d.rho_tor_norm, a).(rho)

    eps = a ./ R

    q = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.q).(rho)

    return 6.921e-18 .* abs.(q) .* R .* ne .* Zeff .* lnLambda_e(ne, Te) ./ (Te .^ 2 .* eps .^ 1.5)
end

function nuistar(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    rho = cp1d.grid.rho_tor_norm
    Zeff = cp1d.zeff

    R = (eqt.profiles_1d.r_outboard + eqt.profiles_1d.r_inboard) / 2.0
    R = interp1d(eqt.profiles_1d.rho_tor_norm, R).(rho)
    a = (eqt.profiles_1d.r_outboard - eqt.profiles_1d.r_inboard) / 2.0
    a = interp1d(eqt.profiles_1d.rho_tor_norm, a).(rho)

    eps = a ./ R

    q = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.q).(rho)
    ne = cp1d.electrons.density
    ni = sum([ion.density for ion in cp1d.ion])
    Ti = cp1d.ion[1].temperature

    Zavg = ne ./ ni

    return 4.90e-18 .* abs.(q) .* R .* ni .* Zeff .^ 4 .* lnLambda_i(ni, Ti, Zavg) ./ (Ti .^ 2 .* eps .^ 1.5)
end

function lnLambda_e(ne, Te)
    return @. 23.5 - log(sqrt(ne / 1e6) * Te ^ (-5.0 / 4.0)) - sqrt(1e-5 + (log(Te) - 2) ^ 2 / 16.0)
end

function lnLambda_i(ni, Ti, Zavg)
    return @. 30.0 - log(Zavg ^ 3 * sqrt(ni) / (Ti ^ 1.5))
end

"""
    nclass_conductivity(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Calculates the neo-classical conductivity in 1/(Ohm*meter) based on the neo 2021 modifcation and stores it in dd
More info see omfit_classes.utils_fusion.py nclass_conductivity function
"""
function nclass_conductivity(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
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


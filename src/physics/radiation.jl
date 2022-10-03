"""
    radiation_losses(sources::IMAS.core_sources)

Evaluate total plasma radiation losses [W] due to both bremsstrahlung and line radiation
Synchlotron radation is not considered since it gets reabsorbed
"""
function radiation_losses(sources::IMAS.core_sources)
    n2i = name_2_index(sources.source)
    radiation_indices = [n2i[name] for name in [:bremsstrahlung, :synchrotron_radiation, :line_radiation]]
    radiation_energy = 0.0
    for source in sources.source
        if source.identifier.index ∈ radiation_indices
            radiation_energy += source.profiles_1d[].electrons.power_inside[end]
        end
    end
    return radiation_energy
end

"""
    bremsstrahlung_source!(dd::IMAS.dd)

Calculates Bremsstrahlung radiation source and modifies dd.core_sources
"""
function bremsstrahlung_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature

    powerDensityBrem = -1.690e-38 .* ne .^ 2 .* cp1d.zeff .* sqrt.(Te)

    index = name_2_index(dd.core_sources.source)[:bremsstrahlung]
    source = resize!(dd.core_sources.source, "identifier.index" => index)
    new_source(source, index, "brem", cp1d.grid.rho_tor_norm, cp1d.grid.volume; electrons_energy=powerDensityBrem)
    return dd
end

"""
    rad_sync(aspect_rat::T, r_min::T, b_ref::T, ne::T, Te::T; reflection_coefficient = 0.8) where {T<:Real}

Synchrotron radiation from Trubnikov, JETP Lett. 16 (1972) 25.0
Transpiled from gacode/tgyro/src/tgyro_rad.f90
"""
function rad_sync(ϵ::T, a::T, B0::T, ne::T, Te::T; reflection_coefficient=0.8) where {T<:Real}
    #---------------------------------------------------
    # MKS to CGS
    aspect_ratio = 1 / ϵ
    r_min = a * 100 # [cm]
    b_ref = B0 * 10000 # [G]
    ne = ne / 1E6 # [1/cm^3]
    Te = Te * 1.0 # [eV]
    e = constants.e * 2.998E9 # [statcoul]
    k = constants.e * 1e7 # [erg / eV]
    m_e = constants.m_e * 1E3 # [g]
    c = constants.c * 1E2  # [cm / s]
    #---------------------------------------------------
    wpe = sqrt(4.0 * pi * ne * e^2 / m_e)
    wce = e * abs(b_ref) / (m_e * c)
    g = k * Te / (m_e * c^2)
    phi = 60.0 * g^1.5 * sqrt((1.0 - reflection_coefficient) * (1.0 + 1.0 / aspect_ratio / sqrt(g)) / (r_min * wpe^2 / c / wce))
    qsync = m_e / (3.0 * pi * c) * g * (wpe * wce)^2 * phi # [erg/cm^3/s]
    return -qsync * 1E-7 * 1E6 #[W/m^3]
end

"""
    synchrotron_source!(dd::IMAS.dd)

Calculates Synchrotron radiation source and modifies dd.core_sources
"""
function synchrotron_source!(dd::IMAS.dd; reflection_coefficient=0.8)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    R = (eq1d.r_outboard + eq1d.r_inboard) / 2.0
    R = interp1d(eq1d.rho_tor_norm, R).(cp1d.grid.rho_tor_norm)
    a = (eq1d.r_outboard - eq1d.r_inboard) / 2.0
    a = interp1d(eq1d.rho_tor_norm, a).(cp1d.grid.rho_tor_norm)
    B0 = abs(@ddtime(eq.vacuum_toroidal_field.b0))
    ϵ = a ./ R

    # Synchrotron radiation
    powerDensitySync = rad_sync.(ϵ, a, B0, ne, Te; reflection_coefficient)

    index = name_2_index(dd.core_sources.source)[:synchrotron_radiation]
    source = resize!(dd.core_sources.source, "identifier.index" => index)
    new_source(source, index, "sync", cp1d.grid.rho_tor_norm, cp1d.grid.volume; electrons_energy=powerDensitySync)
    return dd
end


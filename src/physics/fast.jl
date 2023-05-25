"""
    slowing_down_time(ne::Real, Te::Real, ni::Vector{<:Real}, Ti::Vector{<:Real}, mi::Vector{<:Real}, Zi::Vector{Int}, mf::Real, Zf::Int)

Calculates the slowing down time τ_s [Stix, Plasma Phys. 14 (1972) 367] Eq. 16

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperaturs [eV]

:param mi: list of ion masses [amu]

:param Zi: list of ion charges

:param mf: mass of fast ion [AMU]

:param Zf: fast  ion charge

:return: τ_s: slowing down time
"""
function slowing_down_time(ne::Real, Te::Real, ni::Vector{<:Real}, Ti::Vector{<:Real}, mi::Vector{<:Real}, Zi::Vector{Int}, mf::Real, Zf::Int)
    lnΛ = lnΛ_ei(ne, Te, ni, Ti, mi, Zi)
    ne_cm3 = 1e-6 * ne
    τ_s = 6.27e8 * mf * (Te^1.5) ./ (ne_cm3 * lnΛ * Zf^2)
    return τ_s
end

"""
    slowing_down_time(ne::Real, Te::Real, mf::Real, Zf::Int)

Calculates the slowing down time τ_s for Ti*me/mi < 10Zi^2 eV < Te [Stix, Plasma Phys. 14 (1972) 367] Eq. 16

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param mf: mass of fast ion [AMU]

:param Zf: fast  ion charge

:return: τ_s: slowing down time
"""
function slowing_down_time(ne::Real, Te::Real, mf::Real, Zf::Int)
    lnΛ = lnΛ_ei(ne, Te)
    ne_cm3 = 1e-6 * ne
    τ_s = 6.27e8 * mf * (Te^1.5) / (ne_cm3 * lnΛ * Zf^2)
    return τ_s
end

"""
    _drag_coefficient(n::Real, Z::Int, mf::Real, Zf::Int, lnΛ::Real)

Drag coefficient (Γ) for a fast-ions interacting with a thermal species as defined Eq. 8 in [Gaffey, J. D. (1976). Energetic ion distribution resulting from neutral beam injection in tokamaks. Journal of Plasma Physics, 16(02), 149. doi:10.1017/s0022377800020134]

:param n: density of the thermal species [m^-3]

:param Z: charge of the thermal species

:param mf: fast-ion mass [amu]

:param Zf: fast-ion charge

:param lnΛ: Couloumb logarithm

:return Γ: drag coefficient
"""
function _drag_coefficient(n::Real, Z::Int, mf::Real, Zf::Int, lnΛ::Real)
    n *= 1e-6 #cm^-3
    mf *= constants.m_u # kg
    Γ = (2 * pi * n * (constants.e)^4 * Z^2 * Zf^2 * lnΛ) / (mf^2)
    return Γ
end

"""
    _electron_ion_drag_difference(ne::Real, Te::Real, ni::Vector{<:Real}, Ti::Vector{<:Real}, mi::Vector{<:Real}, Zi::Vector{Int}, Ef::Real, mf::Real, Zf::Int)

Calculates the difference of the electron and ion drag terms in the collision operator defined in Eq. 19 in [Gaffey, J.D (1976). Energetic ion distribution resulting from neutral beam injectioin in tokamaks. Journal of Plasma Physics, 16(02), 149. doi:10.1017/s0022377800020134]

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperaturs [eV]

:param mi: list of ion masses [amu]

:param Zi: list of ion charges

:param Ef: fast ion energy [eV]

:param mf: mass of fast ion [AMU]

:param Zf: fast  ion charge

:return ΔD: drag difference
"""
function _electron_ion_drag_difference(ne::Real, Te::Real, ni::Vector{<:Real}, Ti::Vector{<:Real}, mi::Vector{<:Real}, Zi::Vector{Int}, Ef::Real, mf::Real, Zf::Int)
    m_e = constants.m_e
    m_i = mi .* constants.m_u
    m_f = mf * constants.m_u

    v_f = sqrt(2 * Ef * constants.e / m_f)
    v_e = sqrt(2 * Te * constants.e / m_e)

    lnΛ_fe = lnΛ_ei(ne, Te)
    Γ_fe = _drag_coefficient(ne, -1, mf, Zf, lnΛ_fe)
    electron_drag = ((8 * Γ_fe * m_f) / (3 * sqrt(pi) * m_e * v_e^3)) * v_f^3

    lnΛ_fis = lnΛ_fi(ne, Te, ni, Ti, mi, Zi, v_f / constants.c, mf, Zf; verbose=false)
    Γ_fi = _drag_coefficient.(ni, Zi, mf, Zf, lnΛ_fis)
    ion_drag = 2 * m_f * sum(Γ_fi ./ m_i)

    return electron_drag - ion_drag
end

"""
    critical_energy(ne::Real, Te::Real, ni::Vector{<:Real}, Ti::Vector{<:Real}, mi::Vector{<:Real}, Zi::Vector{Int}, mf::Real, Zf::Int; approximate::Bool=false)

Calculate the critical energy by finding the root of the difference between the electron and ion drag

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperaturs [eV]

:param mi: list of ion masses [amu]

:param Zi: list of ion charges

:param mf: mass of fast ion [AMU]

:param Zf: fast  ion charge

:param approximate: calculate critical energy assuming lnΛ_fe == lnΛ_fi. For DIII-D a correction factor of (lnΛ_fi/lnΛ_fe)^(2/3) ≈ 1.2 can be used.
"""
function critical_energy(ne::Real, Te::Real, ni::Vector{<:Real}, Ti::Vector{<:Real}, mi::Vector{<:Real}, Zi::Vector{Int}, mf::Real, Zf::Int; approximate::Bool=false)
    avg_cmr = sum(ni .* (Zi .^ 2) ./ mi) / ne
    Ec = 14.8 * mf * Te * avg_cmr^(2.0 / 3.0)
    if !(approximate)
        Ec = Roots.find_zero(x -> _electron_ion_drag_difference(ne, Te, ni, Ti, mi, Zi, x, mf, Zf),
            (0.5 * Ec, 2 * Ec))
    end
    return Ec
end

"""
    thermalization_time(v_f, v_c, tau_s)

Calculate thermalization time

:param v_f: fast ion velocity

:param v_c: critical velocity

:param tau_s: slowing down time
"""
function thermalization_time(v_f::Real, v_c::Real, tau_s::Real)
    vf3 = v_f^3
    vc3 = v_c^3
    return tau_s * log((vf3 + vc3) / vc3) / 3
end

"""
    thermalization_time(ne::Real, Te::Real, ni::Vector{<:Real}, Ti::Vector{<:Real}, mi::Vector{<:Real}, Zi::Vector{Int}, Ef::Real, mf::Real, Zf::Int)

Calculate thermalization time of a fast ion with energy Ef and Ti*me/mi < 10Zi^2 eV < Te

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperaturs [eV]

:param mi: list of ion masses [amu]

:param Zi: list of ion charges

:param Ef: fast ion energy [eV]

:param mf: mass of fast ion [AMU]

:param Zf: fast ion charge
"""
function thermalization_time(ne::Real, Te::Real, ni::Vector{<:Real}, Ti::Vector{<:Real}, mi::Vector{<:Real}, Zi::Vector{Int}, Ef::Real, mf::Real, Zf::Int)
    m_f = mf * constants.m_u

    tau_s = slowing_down_time(ne, Te, mf, Zf)
    Ec = critical_energy(ne, Te, ni, Ti, mi, Zi, mf, Zf)

    v_f = sqrt(2 * constants.e * Ef / m_f)
    v_c = sqrt(2 * constants.e * Ec / m_f)

    return thermalization_time(v_f, v_c, tau_s)
end

"""
    fast_density(cs::IMAS.core_sources, cp::IMAS.core_profiles; particle_energy::Real=3.5e6, sourceid::Symbol=:fusion)

Calculates the fast ion density, and adds it to the dd

:param particle_energy: particle energy [eV]

:param sourceid: can be either (:fusion or :nbi)
"""
function fast_density(cs::IMAS.core_sources, cp::IMAS.core_profiles; particle_energy::Real=3.5e6, sourceid::Symbol=:fusion)
    cp1d = cp.profiles_1d[]
    css = cs.source

    ne = cp1d.electrons.density_thermal #electron density profile
    Te = cp1d.electrons.temperature
    Npsi = length(ne)

    if sourceid == :fusion
        particle_specie = "He"
    elseif sourceid == :nbi
        particle_species = ["D", "DT"]
        ion_index = findfirst(ion.label in particle_species for ion in cp1d.ion)
        particle_specie = cp1d.ion[ion_index].label
    else
        error("sourceid can only be `:nbi` or `:fusion`")
    end

    Nions = length(cp1d.ion)
    ion_index = findfirst(ion.label == particle_specie for ion in cp1d.ion)
    if ion_index === nothing
        return
    end

    particle_mass = cp1d.ion[ion_index].element[1].a
    particle_charge = Int(cp1d.ion[ion_index].element[1].z_n)

    ni = zeros(Nions, Npsi)
    Ti = zeros(Nions, Npsi)
    Zi = zeros(Int, Nions)
    mi = zeros(Nions)

    for (idx, ion) in enumerate(cp1d.ion)
        ni[idx, :] = ion.density_thermal
        Ti[idx, :] = ion.temperature
        Zi[idx] = Int(ion.element[1].z_n)
        mi[idx] = ion.element[1].a
    end

    taus = zeros(Npsi)
    taut = zeros(Npsi)
    for i = 1:Npsi
        taus[i] = slowing_down_time(ne[i], Te[i], particle_mass, particle_charge)
        taut[i] = thermalization_time(ne[i], Te[i], ni[:, i], Ti[:, i], mi, Zi, particle_energy, particle_mass, particle_charge)
    end

    cs1ds = findall(sourceid, css)
    cp1d.ion[ion_index].pressure_fast_parallel = zeros(Npsi)
    cp1d.ion[ion_index].pressure_fast_perpendicular = zeros(Npsi)
    cp1d.ion[ion_index].density_fast = zeros(Npsi)

    for cs1d in cs1ds
        qfaste = cs1d.profiles_1d[].electrons.energy
        qfasti = cs1d.profiles_1d[].total_ion_energy
        pressa = taus .* 2.0 ./ 3.0 .* qfaste

        nfast = (qfaste .+ qfasti) ./ (constants.e .* particle_energy) .* taut
        cp1d.ion[ion_index].pressure_fast_parallel += pressa ./ 3.0
        cp1d.ion[ion_index].pressure_fast_perpendicular += pressa ./ 3.0
        cp1d.ion[ion_index].density_fast += nfast
    end
end

"""
    sivukhin_fraction(cp1d::IMAS.core_profiles__profiles_1d, particle_energy::Real, particle_mass::Real)

Compute a low-accuracy but fast approximation to the ion heating fraction (for alpha particles and beam particles)
"""
function sivukhin_fraction(cp1d::IMAS.core_profiles__profiles_1d, particle_energy::Real, particle_mass::Real)
    Te = cp1d.electrons.temperature
    ne = cp1d.electrons.density
    rho = cp1d.grid.rho_tor_norm

    tp = typeof(promote(Te[1], ne[1], rho[1])[1])
    c_a = zeros(tp, length(rho))
    for ion in cp1d.ion
        ni = ion.density_thermal
        Zi = avgZ(ion.element[1].z_n, ion.temperature)
        mi = ion.element[1].a
        c_a .+= (ni ./ ne) .* Zi .^ 2 ./ (mi ./ particle_mass)
    end

    W_crit = Te .* (4.0 .* sqrt.(constants.m_e / (constants.m_p * particle_mass)) ./ (3.0 * sqrt(pi) .* c_a)) .^ (-2.0 / 3.0)
    ion_to_electron_fraction = similar(W_crit)
    x = particle_energy ./ W_crit
    for (idx, x_i) in enumerate(x)
        if x_i > 4.0
            # Large-x asymptotic formula
            f = (2 * pi / 3) / sin(2 * pi / 3) - 2.0 / sqrt(x_i) + 0.5 / x_i^2
            ion_to_electron_fraction[idx] = f / x_i
        elseif x_i < 0.1
            # Small-x asymptotic series
            ion_to_electron_fraction[idx] = 1.0 - 0.4 * x_i^1.5
        else
            y = x_i .* LinRange(0, 1, 12)
            f = integrate(y, 1.0 ./ (1.0 .+ y .^ 1.5))
            ion_to_electron_fraction[idx] = f / x_i
        end
    end

    return ion_to_electron_fraction
end
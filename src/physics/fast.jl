"""
    estrada_I_integrals(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real,  mf::Real, Zf::Int)
 
Returns solution to i2 and i4 integrals from  [Estrada et al.,  Phys of Plasm. 13, 112303 (2006)] Eq. 9 & 10

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperatures [eV]

:param mi: list of ion masses [AMU]

:param Zi: list of ion charges

:param Ef: fast ion energy [eV]

:param mf: mass of fast ion [AMU]

:param Zf: fast ion charge

:return: i2 and i4
"""
                        
function estrada_I_integrals(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real,  mf::Real, Zf::Int)
    Ec = critical_energy(ne, Te, ni, Ti, mi, Zi, mf, Zf)
    i2,i4 = estrada_I_integrals(Ec, Ef)
    return i2,i4
end

"""
    estrada_I_integrals(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real,  mf::Real, Zf::Int)
 
Returns solution to i2 and i4 integrals from  [Estrada et al.,  Phys of Plasm. 13, 112303 (2006)] Eq. 9 & 10

:param Ef: fast ion energy [eV]

:param Ec: critical energy [eV]

:return: i2 and i4
"""
function estrada_I_integrals(Ec::Real, Ef::Real)
    a =  sqrt.(Ec./Ef)
    i2 = (1/3.0) .* log.((1 .+ a.^3) ./ (a.^3))
    i4 = 0.5 .- a.^2 .* ((1/6.0) .* log.((1 .- a .+ a.^2) ./ (1 .+ a).^2) .+
         1 ./ sqrt(3.0) .* (atan.((2 .- a) ./ (a .* sqrt(3.0))) .+ pi/6))
    return i2,i4
end


"""
    slowing_down_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Int)

Calculates the slowing down time τ_s [Stix, Plasma Phys. 14 (1972) 367] Eq. 16

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperatures [eV]

:param mi: list of ion masses [AMU]

:param Zi: list of ion charges

:param mf: mass of fast ion [AMU]

:param Zf: fast ion charge

:return: τ_s: slowing down time
"""
function slowing_down_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Int)
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

:param Zf: fast ion charge

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

:param mf: fast-ion mass [AMU]

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
    _electron_ion_drag_difference(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real, mf::Real, Zf::Int)

Calculates the difference of the electron and ion drag terms in the collision operator defined in Eq. 19 in [Gaffey, J.D (1976). Energetic ion distribution resulting from neutral beam injectioin in tokamaks. Journal of Plasma Physics, 16(02), 149. doi:10.1017/s0022377800020134]

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperatures [eV]

:param mi: list of ion masses [AMU]

:param Zi: list of ion charges

:param Ef: fast ion energy [eV]

:param mf: mass of fast ion [AMU]

:param Zf: fast ion charge

:return ΔD: drag difference
"""
function _electron_ion_drag_difference(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real, mf::Real, Zf::Int)
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
    critical_energy(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Int; approximate::Bool=false)

Calculate the critical energy by finding the root of the difference between the electron and ion drag

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperatures [eV]

:param mi: list of ion masses [AMU]

:param Zi: list of ion charges

:param mf: mass of fast ion [AMU]

:param Zf: fast ion charge

:param approximate: calculate critical energy assuming lnΛ_fe == lnΛ_fi. For DIII-D a correction factor of (lnΛ_fi/lnΛ_fe)^(2/3) ≈ 1.2 can be used.
"""
function critical_energy(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Int; approximate::Bool=false)
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
    thermalization_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real, mf::Real, Zf::Int)

Calculate thermalization time of a fast ion with energy Ef and Ti*me/mi < 10Zi^2 eV < Te

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperatures [eV]

:param mi: list of ion masses [AMU]

:param Zi: list of ion charges

:param Ef: fast ion energy [eV]

:param mf: mass of fast ion [AMU]

:param Zf: fast ion charge
"""
function thermalization_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real, mf::Real, Zf::Int)
    m_f = mf * constants.m_u

    tau_s = slowing_down_time(ne, Te, mf, Zf)
    Ec = critical_energy(ne, Te, ni, Ti, mi, Zi, mf, Zf)

    v_f = sqrt(2 * constants.e * Ef / m_f)
    v_c = sqrt(2 * constants.e * Ec / m_f)

    return thermalization_time(v_f, v_c, tau_s)
end

"""
    fast_particles!(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; verbose::Bool=false)

Calculates the core_profiles fast ion density and pressures resulting from fast ion sources (fusion, nbi)
"""
function fast_particles!(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; verbose::Bool=false)
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature

    Npsi = length(ne)
    Nions = length(cp1d.ion)

    # prepare inputs for slowing_down_time() and thermalization_time() functions
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

    # empty cp1d pressures (expressions)
    empty!(cp1d, :pressure)
    empty!(cp1d, :pressure_parallel)
    empty!(cp1d, :pressure_perpendicular)
    empty!(cp1d, :pressure_ion_total)
    empty!(cp1d, :pressure_thermal)
    # empty all cp1d fast-ion related quantities (expressions)
    for ion in cp1d.ion
        freeze!(ion, :pressure_thermal)
        empty!(ion, :pressure)
        freeze!(ion, :density_thermal)
        empty!(ion, :density)
    end
    # zero out all cp1d fast-ion related quantities
    for ion in cp1d.ion
        ion.pressure_fast_parallel = zeros(Npsi)
        ion.pressure_fast_perpendicular = zeros(Npsi)
        ion.density_fast = zeros(Npsi)
    end

    # go through sources and look for ones that have ion particles source at given energy
    taus = zeros(Npsi)
    taut = zeros(Npsi)
    i4 = zeros(Npsi)
    for source in cs.source
        for sion in source.profiles_1d[].ion
            if !ismissing(sion, :particles) && !ismissing(sion, :fast_particles_energy)

                particle_mass = sion.element[1].a
                particle_charge = Int(sion.element[1].z_n)
                particle_energy = sion.fast_particles_energy

                # find the corresponding thermal ion in core_profiles (we use z_n and a since that's more reliable than labels)
                # NOTE: contribution of non-thermal components in core_profiles comes in through pressure_fast and density_fast
                cindex = findfirst(cion ->
                        (Int(round(cion.element[1].z_n)) == Int(round(particle_charge)) && Int(round(cion.element[1].a)) == Int(round(particle_mass))) ||
                            (Int(round(particle_charge)) == 1 && Int(round(particle_mass)) == 2 && cion.label == "DT") ||
                            (Int(round(particle_charge)) == 1 && Int(round(particle_mass)) == 3 && cion.label == "DT"), cp1d.ion)

                if cindex === nothing
                    if verbose
                        println("$(sion.label) --> ?")
                    end
                else
                    cion = cp1d.ion[cindex]
                    if verbose
                        println("$(sion.label) --> $(cion.label) @ $(sion.fast_particles_energy)")
                    end

                    taus .*= 0.0
                    taut .*= 0.0
                    for i = 1:Npsi
                        taus[i] = slowing_down_time(ne[i], Te[i], particle_mass, particle_charge)
                        taut[i] = @views thermalization_time(ne[i], Te[i], ni[:, i], Ti[:, i], mi, Zi, particle_energy, particle_mass, particle_charge)

                        i2tmp,i4tmp = estrada_I_integrals(ne[i], Te[i], ni[:, i], Ti[:, i], mi, Zi, particle_energy, particle_mass, particle_charge)
                        i4[i] = i4tmp
                    end
                    
                    pressa = i4 .* taus .* 2.0 ./ 3.0 .* (sion.particles .* particle_energy .* constants.e)
                    cion.pressure_fast_parallel += pressa ./ 3.0
                    cion.pressure_fast_perpendicular += pressa ./ 3.0
                    cion.density_fast += sion.particles .* taut
                end
            end
        end
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
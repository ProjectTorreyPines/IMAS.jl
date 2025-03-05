"""
    ENERGY_FLUX_e::T
    ENERGY_FLUX_i::T
    PARTICLE_FLUX_e::T
    PARTICLE_FLUX_i::Vector{T}
    STRESS_TOR_i::T

Structure used to store fluxes information from GA code in Gyrobohm units
"""
struct flux_solution{T<:Real}
    ENERGY_FLUX_e::T
    ENERGY_FLUX_i::T
    PARTICLE_FLUX_e::T
    PARTICLE_FLUX_i::Vector{T}
    STRESS_TOR_i::T
end

function Base.show(io::IO, ::MIME"text/plain", sol::flux_solution)
    txt = """
    Qe =  $(sol.ENERGY_FLUX_e)
    Qi =  $(sol.ENERGY_FLUX_i)
    Γe =  $(sol.PARTICLE_FLUX_e)
    Γi = [$(join( map(string, sol.PARTICLE_FLUX_i),", "))]
    Πi =  $(sol.STRESS_TOR_i)
    """
    return print(io, txt)
end

"""
    c_s(cp1d::IMAS.core_profiles__profiles_1d)

Sounds speed in [cm/s]
"""
function c_s(cp1d::IMAS.core_profiles__profiles_1d)
    return sqrt.(cgs.k .* cp1d.electrons.temperature ./ cgs.md)
end

"""
    rho_s(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

sound gyro radius in [cm]
"""
function rho_s(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    eqt1d = eqt.profiles_1d
    bu = interp1d(eqt1d.rho_tor_norm, abs.(bunit(eqt1d)) .* cgs.T_to_Gauss).(cp1d.grid.rho_tor_norm)
    return c_s(cp1d) ./ (cgs.e .* bu) .* (cgs.md .* cgs.c)
end

"""
    r_min_core_profiles(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Geometric minor radius in [cm] evaluated
"""
function r_min_core_profiles(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    eq1d = eqt.profiles_1d
    return interp1d(eq1d.rho_tor_norm, cgs.m_to_cm * 0.5 * (eq1d.r_outboard - eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
end

##### Gyrobohm normalizations from gacode
"""
    gyrobohm_energy_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Gyrobohm energy flux
"""
function gyrobohm_energy_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    return cp1d.electrons.density_thermal ./ cgs.m³_to_cm³ .* cgs.k .* cp1d.electrons.temperature .*
           c_s(cp1d) .* (rho_s(cp1d, eqt) ./ (eqt.boundary.minor_radius .* cgs.m_to_cm)) .^ 2 .* cgs.Erg_to_J .*
           cgs.m²_to_cm²
end

"""
    gyrobohm_particle_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Gyrobohm particle flux
"""
function gyrobohm_particle_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    norm = mks.e .* cp1d.electrons.temperature
    return gyrobohm_energy_flux(cp1d, eqt) ./ norm
end

"""
    gyrobohm_momentum_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Gyrobohm momentum flux
"""
function gyrobohm_momentum_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    return cp1d.electrons.density_thermal ./ cgs.m³_to_cm³ .* cgs.k .* cp1d.electrons.temperature .*
           eqt.boundary.minor_radius .* cgs.m_to_cm .* (rho_s(cp1d, eqt) ./ (eqt.boundary.minor_radius .* cgs.m_to_cm)) .^ 2 .* cgs.Erg_to_J .*
           cgs.m²_to_cm²
end

"""
    volume_prime_miller_correction(eqt::IMAS.equilibrium__time_slice)

Correction to account for transformation from Miller r grid in GA code equilibrium to Psi grid in FUSE equilibrium
"""
function volume_prime_miller_correction(eqt::IMAS.equilibrium__time_slice)
    a_minor = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    return gradient(a_minor, eqt.profiles_1d.volume) ./ eqt.profiles_1d.surface
end

"""
    flux_gacode_to_fuse(
        flux_types::Tuple{Vararg{Symbol}},
        flux_solutions::Vector{<:IMAS.flux_solution},
        m1d::IMAS.core_transport__model___profiles_1d,
        eqt::IMAS.equilibrium__time_slice,
        cp1d::core_profiles__profiles_1d
    )

Normalizes specified transport fluxes output by GA code via gyrobohm normalization and Miller volume correction
"""
function flux_gacode_to_fuse(
    flux_types::Tuple{Vararg{Symbol}},
    flux_solutions::Vector{<:IMAS.flux_solution},
    m1d::IMAS.core_transport__model___profiles_1d{T},
    eqt::IMAS.equilibrium__time_slice{T},
    cp1d::core_profiles__profiles_1d{T}) where {T<:Real}

    rho_eq_idxs = [argmin(abs.(eqt.profiles_1d.rho_tor_norm .- rho)) for rho in m1d.grid_flux.rho_tor_norm]
    rho_cp_idxs = [argmin(abs.(cp1d.grid.rho_tor_norm .- rho)) for rho in m1d.grid_flux.rho_tor_norm]

    vprime_miller = volume_prime_miller_correction(eqt)

    if :ion_energy_flux in flux_types
        m1d.total_ion_energy.flux = gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idxs] .* [f.ENERGY_FLUX_i for f in flux_solutions] .* vprime_miller[rho_eq_idxs]
    end
    if :electron_energy_flux in flux_types
        m1d.electrons.energy.flux = gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idxs] .* [f.ENERGY_FLUX_e for f in flux_solutions] .* vprime_miller[rho_eq_idxs]
    end
    if :electron_particle_flux in flux_types
        m1d.electrons.particles.flux = gyrobohm_particle_flux(cp1d, eqt)[rho_cp_idxs] .* [f.PARTICLE_FLUX_e for f in flux_solutions] .* vprime_miller[rho_eq_idxs]
    end
    if :momentum_flux in flux_types
        m1d.momentum_tor.flux = gyrobohm_momentum_flux(cp1d, eqt)[rho_cp_idxs] .* [f.STRESS_TOR_i for f in flux_solutions] .* vprime_miller[rho_eq_idxs]
    end

    if :ion_particle_flux in flux_types
        for (kk, ion) in enumerate(cp1d.ion)
            ion = resize!(m1d.ion, "element[1].a" => ion.element[1].z_n, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label)
            ion.particles.flux = gyrobohm_particle_flux(cp1d, eqt)[rho_cp_idxs] .* [pick_ion_flux(f.PARTICLE_FLUX_i, kk) for f in flux_solutions] .* vprime_miller[rho_eq_idxs]
        end
    end
end

"""
    pick_ion_flux(ion_fluxes::AbstractVector{T}, kk::Int) where {T<:Real}

Select which ion flux to take
"""
function pick_ion_flux(ion_fluxes::AbstractVector{T}, kk::Int) where {T<:Real}
    if isempty(ion_fluxes)
        return T(0.0)
    elseif kk <= length(ion_fluxes)
        return ion_fluxes[kk]
    else
        return ion_fluxes[end]
    end
end

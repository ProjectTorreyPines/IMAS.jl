const gacode_units = (
    e=4.8032e-10, # stacoul
    k=1.6022e-12, # erg/eV
    me=9.1094e-28, # g
    c=2.9979e10, # cm/s
    mp=1.6726e-24, # g
    m_to_cm=1e2,
    T_to_Gauss=1e4,
    Erg_to_J=1e-7,
    m³_to_cm³=1e-6
)

struct flux_solution{T<:Real}
    PARTICLE_FLUX_e::T
    STRESS_TOR_i::T
    ENERGY_FLUX_e::T
    ENERGY_FLUX_i::T
end

function Base.show(io::IO, sol::flux_solution)
    txt = """
    Γe = $(sol.PARTICLE_FLUX_e)
    Πi = $(sol.STRESS_TOR_i)
    Qe = $(sol.ENERGY_FLUX_e)
    Qi = $(sol.ENERGY_FLUX_i)
    """
    print(io, txt)
end

function c_s(cp1d::IMAS.core_profiles__profiles_1d)
    return sqrt.(gacode_units.k .* cp1d.electrons.temperature ./ (cp1d.ion[1].element[1].a .* gacode_units.mp))
end

function rho_s(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    bunit = interp1d(eqt.profiles_1d.rho_tor_norm, abs.(IMAS.bunit(eqt)) .* gacode_units.T_to_Gauss).(cp1d.grid.rho_tor_norm)
    return c_s(cp1d) ./ (gacode_units.e .* bunit) .* (cp1d.ion[1].element[1].a .* gacode_units.mp .* gacode_units.c)
end

function r_min_core_profiles(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    eq1d = eqt.profiles_1d
    return interp1d(eq1d.rho_tor_norm, gacode_units.m_to_cm * 0.5 * (eq1d.r_outboard - eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
end

##### Gyrobohm normalizations from gacode
function gyrobohm_energy_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    return cp1d.electrons.density_thermal .* gacode_units.m³_to_cm³ .* gacode_units.k .* cp1d.electrons.temperature .*
           c_s(cp1d) .* (rho_s(cp1d, eqt) ./ (eqt.boundary.minor_radius .* gacode_units.m_to_cm)) .^ 2 .* gacode_units.Erg_to_J .* gacode_units.m_to_cm^2
end

function gyrobohm_particle_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    norm =  constants.e .* cp1d.electrons.temperature
    return gyrobohm_energy_flux(cp1d, eqt) ./ norm
end

function gyrobohm_momentum_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    return cp1d.electrons.density_thermal .* gacode_units.m³_to_cm³ .* gacode_units.k .* cp1d.electrons.temperature .*
           eqt.boundary.minor_radius .* gacode_units.m_to_cm .* (rho_s(cp1d, eqt) ./ (eqt.boundary.minor_radius .* gacode_units.m_to_cm)) .^ 2 .* gacode_units.Erg_to_J .* gacode_units.m_to_cm^2
end

"""
    volume_prime_miller_correction(eqt::IMAS.equilibrium__time_slice)

Correction to account for transformation from Miller r grid in GA code equilibrium to Psi grid in FUSE equilibrium
"""
function volume_prime_miller_correction(eqt::IMAS.equilibrium__time_slice)
    a_minor = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    return IMAS.gradient(a_minor, eqt.profiles_1d.volume) ./ eqt.profiles_1d.surface
end

"""
    flux_gacode_to_fuse(flux_types::Vector{Symbol}, flux_solutions::Vector{<:IMAS.flux_solution}, m1d::IMAS.core_transport__model___profiles_1d, eqt::IMAS.equilibrium__time_slice, cp1d::core_profiles__profiles_1d)

Normalizes specified transport fluxes output by GA code via gyrobohm normalization and Miller volume correction
"""
function flux_gacode_to_fuse(flux_types::Vector{Symbol}, flux_solutions::Vector{<:IMAS.flux_solution}, m1d::IMAS.core_transport__model___profiles_1d, eqt::IMAS.equilibrium__time_slice, cp1d::core_profiles__profiles_1d)

    rho_eq_idxs = [argmin(abs.(eqt.profiles_1d.rho_tor_norm .- rho)) for rho in m1d.grid_flux.rho_tor_norm]
    rho_cp_idxs = [argmin(abs.(cp1d.grid.rho_tor_norm .- rho)) for rho in m1d.grid_flux.rho_tor_norm]

    ga_tr = Dict(
        :ion_energy_flux=>[m1d.total_ion_energy, gyrobohm_energy_flux, :ENERGY_FLUX_i],
        :electron_energy_flux=>[m1d.electrons.energy, gyrobohm_energy_flux, :ENERGY_FLUX_e],
        :electron_particle_flux=>[m1d.electrons.particles, gyrobohm_particle_flux, :PARTICLE_FLUX_e],
        :momentum_flux=>[m1d.momentum_tor, gyrobohm_momentum_flux, :STRESS_TOR_i])

    vprime_miller = volume_prime_miller_correction(eqt)    

    for flux_type in flux_types
        result = ga_tr[flux_type][2](cp1d, eqt)[rho_cp_idxs] .* 
        [getproperty(f, ga_tr[flux_type][3]) for f in flux_solutions] .* vprime_miller[rho_eq_idxs]

        setproperty!(ga_tr[flux_type][1], :flux, result)
    end

end
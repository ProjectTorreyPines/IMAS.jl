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

##### GA code to FUSE flux normalizations: gyrobohm normalization + Miller volume correction 
# Volume correction accounts for transformation from Miller r grid in GA code equilibrium to Psi grid in FUSE equilibrium
function energy_flux_gacode_to_fuse(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice, rho::Real)
    rho_cp_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho))
    rho_eq_idx = argmin(abs.(eqt.profiles_1d.rho_tor_norm .- rho))

    a_minor = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    volume_prime_miller_correction = IMAS.gradient(a_minor, eqt.profiles_1d.volume) ./ eqt.profiles_1d.surface
    return volume_prime_miller_correction[rho_eq_idx] * gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx]
end

function particle_flux_gacode_to_fuse(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice, rho::Real)
    rho_cp_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho))
    rho_eq_idx = argmin(abs.(eqt.profiles_1d.rho_tor_norm .- rho))

    a_minor = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    volume_prime_miller_correction = IMAS.gradient(a_minor, eqt.profiles_1d.volume) ./ eqt.profiles_1d.surface
    return volume_prime_miller_correction[rho_eq_idx] * gyrobohm_particle_flux(cp1d, eqt)[rho_cp_idx]
end

function momentum_flux_gacode_to_fuse(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice, rho::Real)
    rho_cp_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho))
    rho_eq_idx = argmin(abs.(eqt.profiles_1d.rho_tor_norm .- rho))

    a_minor = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    volume_prime_miller_correction = IMAS.gradient(a_minor, eqt.profiles_1d.volume) ./ eqt.profiles_1d.surface
    return volume_prime_miller_correction[rho_eq_idx] * gyrobohm_momentum_flux(cp1d, eqt)[rho_cp_idx]
end
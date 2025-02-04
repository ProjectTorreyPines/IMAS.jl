document[Symbol("Physics profiles")] = Symbol[]

"""
    pressure_thermal(cp1d::IMAS.core_profiles__profiles_1d)

core_profiles thermal pressure
"""
function pressure_thermal(cp1d::IMAS.core_profiles__profiles_1d)
    p = pressure_thermal(cp1d.electrons)
    p .+= pressure_thermal(cp1d.ion)
    return p
end

"""
    pressure_thermal(cp1de::IMAS.core_profiles__profiles_1d___electrons)

electrons thermal pressure
"""
function pressure_thermal(cp1de::IMAS.core_profiles__profiles_1d___electrons)
    return cp1de.temperature .* cp1de.density_thermal .* mks.e
end

"""
    pressure_thermal(ion::IMAS.core_profiles__profiles_1d___ion)

ion thermal pressure
"""
function pressure_thermal(ion::IMAS.core_profiles__profiles_1d___ion)
    return ion.temperature .* ion.density_thermal .* mks.e
end

"""
    pressure_thermal(cp1di::IMAS.IDSvector{IMAS.core_profiles__profiles_1d___ion{T}}) where {T<:Real}

thermal pressure for all ions
"""
function pressure_thermal(cp1di::IMAS.IDSvector{IMAS.core_profiles__profiles_1d___ion{T}}) where {T<:Real}
    p = cp1di[1].temperature .* 0.0
    for ion in cp1di
        p .+= ion.temperature .* ion.density_thermal
    end
    return p .* mks.e
end

@compat public pressure_thermal
push!(document[Symbol("Physics profiles")], :pressure_thermal)

"""
    beta_tor_thermal_norm(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d)

Normalised toroidal beta from thermal pressure only, defined as 100 * beta_tor_thermal * a[m] * B0 [T] / ip [MA]
"""
function beta_tor_thermal_norm(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d)
    return beta_tor(eq, cp1d; norm=true, thermal=true)
end

@compat public beta_tor_thermal_norm
push!(document[Symbol("Physics profiles")], :beta_tor_thermal_norm)

"""
    beta_tor_norm(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d)

Normalised toroidal beta from total pressure, defined as 100 * beta_tor * a[m] * B0 [T] / ip [MA]
"""
function beta_tor_norm(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d)
    return beta_tor(eq, cp1d; norm=true, thermal=false)
end

@compat public beta_tor_norm
push!(document[Symbol("Physics profiles")], :beta_tor_norm)

"""
    beta_tor(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d; norm::Bool, thermal::Bool)

Toroidal beta, defined as the volume-averaged total perpendicular pressure divided by (B0^2/(2*mu0)), i.e. beta_toroidal = 2 mu0 int(p dV) / V / B0^2
"""
function beta_tor(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d; norm::Bool, thermal::Bool)
    dd = top_dd(eq)
    B0 = get_time_array(eq.vacuum_toroidal_field, :b0, dd.global_time, :constant)
    ip = Ip(cp1d)

    if thermal
        pressure = cp1d.pressure_thermal
    else
        pressure = cp1d.pressure
    end
    @assert !any(isnan.(pressure))

    volume = cp1d.grid.volume
    pressure_avg = trapz(volume, pressure) / volume[end]

    beta_tor = 2.0 * mks.μ_0 * pressure_avg / B0^2

    if norm
        eqt = eq.time_slice[dd.global_time]
        out = beta_tor * eqt.boundary.minor_radius * abs(B0) / abs(ip / 1e6) * 1.0e2
    else
        out = beta_tor
    end

    return out
end

@compat public beta_tor
push!(document[Symbol("Physics profiles")], :beta_tor)

function list_ions!(ct::IMAS.core_transport, ions::Vector{Symbol}; time0::Float64)
    for model in ct.model
        if isempty(model.profiles_1d)
            continue
        end
        ct1d = model.profiles_1d[time0]
        for ion in ct1d.ion
            push!(ions, Symbol(ion.label))
        end
    end
    return ions
end

function list_ions!(cp::IMAS.core_profiles, ions::Vector{Symbol}; time0::Float64)
    if isempty(cp.profiles_1d)
        return ions
    end
    cp1d = cp.profiles_1d[time0]
    for ion in cp1d.ion
        push!(ions, Symbol(ion.label))
    end
    return ions
end

function list_ions!(cs::IMAS.core_sources, ions::Vector{Symbol}; time0::Float64)
    for source in cs.source
        if isempty(source.profiles_1d)
            continue
        end
        sc1d = source.profiles_1d[time0]
        for ion in sc1d.ion
            push!(ions, Symbol(ion.label))
        end
    end
    return ions
end

"""
    list_ions(ids::IDS, idss::Vararg{<:IDS}; time0::Float64)

List of ions mentioned in multiple IDSs at a given time
"""
function list_ions(ids1::IDS, idss::Vararg{<:IDS}; time0::Float64)
    ions = Symbol[]
    for ids in [ids1; idss...]
        append!(ions, list_ions!(ids, ions; time0))
        sort!(unique!(ions))
    end
    return ions
end

@compat public list_ions
push!(document[Symbol("Physics profiles")], :list_ions)

"""
    ion_element!(
        ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
        ion_symbol::Symbol)

Fills the `ion.element` structure with the a and z_n information, also updates the `ion.label`
"""
function ion_element!(
    ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
    ion_symbol::Symbol;
    fast::Bool=false)

    z_n, a, label = ion_properties(ion_symbol; fast)

    element = resize!(ion.element, 1)[1]

    element.z_n = z_n
    element.a = a
    ion.label = label

    return ion.element
end

function ion_element!(
    ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
    ion_string::AbstractString;
    fast::Bool=false)
    return ion_element!(ion, Symbol(ion_string); fast)
end

function ion_element!(
    ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
    ion_z::Int, ion_a::Float64; fast::Bool=false)
    if ion_z == 1 && ion_a == 1.0
        ion_symbol = :H
    elseif ion_z == 1 && ion_a == 2.0
        ion_symbol = :D
    elseif ion_z == 1 && ion_a == 2.5
        ion_symbol = :DT
    elseif ion_z == 1 && ion_a == 3.0
        ion_symbol = :T
    elseif ion_z == 2 && ion_a == 4.0
        ion_symbol = :α
    else
        ion_symbol = elements[Int(ion_z)].symbol
    end
    return ion_element!(ion, ion_symbol; fast)
end

function ion_element!(
    ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
    ion_z::Int;
    fast::Bool=false)
    return ion_element!(ion, elements[ion_z].symbol; fast)
end

@compat public ion_element!
push!(document[Symbol("Physics profiles")], :ion_element!)

"""
    ion_properties( ion_symbol::Symbol; fast::Bool=false)

Returns named tuple with z_n, a, and label information of a given ion
"""
function ion_properties(ion_symbol::Symbol; fast::Bool=false)
    # H isotopes
    if ion_symbol ∈ (:H, :H1)
        z_n = 1.0
        a = 1.00797
        label = "H"

    elseif ion_symbol ∈ (:D, :H2)
        z_n = 1.0
        a = 2.014
        label = "D"

    elseif ion_symbol ∈ (:T, :H3)
        z_n = 1.0
        a = 3.016
        label = "T"

    elseif ion_symbol == :α
        z_n = 2.0
        a = 4.007
        label = "α"

    elseif ion_symbol == :DT
        z_n = 1.0
        a = (2.014 + 3.016) / 2.0
        label = "DT"

    else
        # all other ions
        ion_name, ion_a = match(r"(.*?)(\d*)$", string(ion_symbol))
        element_ion = elements[Symbol(ion_name)]
        z_n = float(element_ion.number)
        if isempty(ion_a)
            a = element_ion.atomic_mass.val
        else
            z = element_ion.number
            n = parse(Int, ion_a) - z
            a = atomic_mass(z, n)
        end
        label = "$(ion_name)$(Int(floor(a)))"
    end

    if fast
        label = "$(label)_fast"
    end

    return (z_n=z_n, a=a, label=label)
end

@compat public ion_properties
push!(document[Symbol("Physics profiles")], :ion_properties)

"""
    binding_energy(Z::Int, N::Int)

Return the estimate binding energy in MeV
"""
function binding_energy(Z::Int, N::Int)
    # Constants (in MeV)
    a_v = 15.8
    a_s = 18.3
    a_c = 0.714
    a_a = 23.2

    # Total number of nucleons
    A = Z + N

    # Calculate the pairing term
    if A % 2 != 0
        δ = 0.0
    elseif Z % 2 == 0
        δ = 12.0 / sqrt(A)
    else
        δ = -12.0 / sqrt(A)
    end

    # Calculate the binding energy (in MeV)
    E = a_v * A - a_s * A^(2 / 3) - a_c * Z * (Z - 1) / A^(1 / 3) - a_a * (N - Z)^2 / A + δ

    return E
end

@compat public binding_energy
push!(document[Symbol("Physics profiles")], :binding_energy)

"""
    atomic_mass(Z::Int, N::Int)

Returns the estimated nucleus mass including the estimated effect of binding energy
"""
function atomic_mass(Z::Int, N::Int)
    # Mass of proton and neutron (in amu)
    mass_proton = mks.m_p / mks.m_u
    mass_neutron = mks.m_n / mks.m_u

    # Conversion factor from MeV/c^2 to atomic mass units (amu)
    mev_to_amu = mks.m_u * mks.c^2 / (mks.e * 1E6)

    # Estimate binding energy
    E = binding_energy(Z, N)

    # Calculate mass of nucleus (in amu)
    mass = Z * mass_proton + N * mass_neutron - E / mev_to_amu

    return mass
end

@compat public atomic_mass
push!(document[Symbol("Physics profiles")], :atomic_mass)

"""
    energy_thermal(cp1d::IMAS.core_profiles__profiles_1d)

Calculates the thermal stored energy
"""
function energy_thermal(cp1d::IMAS.core_profiles__profiles_1d)
    return 3.0 / 2.0 * trapz(cp1d.grid.volume, cp1d.pressure_thermal)
end

@compat public energy_thermal
push!(document[Symbol("Physics profiles")], :energy_thermal)

"""
    energy_thermal_ped(cp1d::IMAS.core_profiles__profiles_1d, su::IMAS.summary)

Calculates the pedestal contribution to the thermal stored energy by integrating over entire domain but with pedestal pressure in the core
"""
function energy_thermal_ped(cp1d::IMAS.core_profiles__profiles_1d, su::IMAS.summary)
    rho_index = argmin(abs.(cp1d.grid.rho_tor_norm .- su.local.pedestal.position.rho_tor_norm))
    pressure = cp1d.pressure_thermal
    pressure[1:rho_index] .= cp1d.pressure_thermal[rho_index]
    return 3.0 / 2.0 * trapz(cp1d.grid.volume, pressure)
end

@compat public energy_thermal_ped
push!(document[Symbol("Physics profiles")], :energy_thermal_ped)

"""
    tau_e_thermal(cp1d::IMAS.core_profiles__profiles_1d, sources::IMAS.core_sources; subtract_radiation_losses::Bool=true)

Evaluate thermal energy confinement time
"""
function tau_e_thermal(cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.core_sources; subtract_radiation_losses::Bool=true)
    dd = top_dd(cp1d)
    total_source = total_sources(cs, cp1d; time0=dd.global_time, fields=[:power_inside, :total_ion_power_inside])
    total_power_inside = total_source.electrons.power_inside[end] + total_source.total_ion_power_inside[end]
    if subtract_radiation_losses
        total_power_inside -= radiation_losses(cs)
    end
    total_power_inside = max(0.0, total_power_inside)
    return energy_thermal(cp1d) / total_power_inside
end

"""
    tau_e_thermal(dd::IMAS.dd; time0::Float64=dd.global_time, subtract_radiation_losses::Bool=true)
"""
function tau_e_thermal(dd::IMAS.dd; time0::Float64=dd.global_time, subtract_radiation_losses::Bool=true)
    return tau_e_thermal(dd.core_profiles.profiles_1d[time0], dd.core_sources; subtract_radiation_losses)
end

@compat public tau_e_thermal
push!(document[Symbol("Physics profiles")], :tau_e_thermal)

"""
    tau_e_h98(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.plot_core_sources; subtract_radiation_losses::Bool=true)

H98y2 ITER elmy H-mode confinement time scaling

NOTE: H98y2 uses aereal elongation

See Table 5 in https://iopscience.iop.org/article/10.1088/0029-5515/39/12/302/pdf and https://iopscience.iop.org/article/10.1088/0029-5515/48/9/099801/pdf for additional correction with plasma_volume
"""
function tau_e_h98(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.core_sources; subtract_radiation_losses::Bool=true)
    dd = top_dd(cp1d)
    total_source = total_sources(cs, cp1d; time0=dd.global_time, fields=[:power_inside, :total_ion_power_inside])
    total_power_inside = total_source.electrons.power_inside[end] + total_source.total_ion_power_inside[end]
    if subtract_radiation_losses
        total_power_inside -= radiation_losses(cs)
    end
    total_power_inside = max(0.0, total_power_inside)

    isotope_factor =
        trapz(cp1d.grid.volume, sum(ion.density .* ion.element[1].a for ion in cp1d.ion if ion.element[1].z_n == 1.0)) /
        trapz(cp1d.grid.volume, sum(ion.density for ion in cp1d.ion if ion.element[1].z_n == 1.0))

    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

    κ_areal = areal_elongation(eqt)

    ne_line = geometric_midplane_line_averaged_density(eqt, cp1d)

    tau98 = (
        0.0562 *
        abs(eqt.global_quantities.ip / 1e6)^0.93 *
        abs(B0)^0.15 *
        (total_power_inside / 1e6)^-0.69 *
        (ne_line / 1e19)^0.41 *
        isotope_factor^0.19 *
        R0^1.97 *
        (R0 / eqt.boundary.minor_radius)^-0.58 *
        κ_areal^0.78
    )
    return tau98
end

"""
    tau_e_h98(dd::IMAS.dd; time0::Float64=dd.global_time, subtract_radiation_losses::Bool=true)
"""
function tau_e_h98(dd::IMAS.dd; time0::Float64=dd.global_time, subtract_radiation_losses::Bool=true)
    eqt = dd.equilibrium.time_slice[time0]
    cp1d = dd.core_profiles.profiles_1d[time0]
    cs = dd.core_sources
    return tau_e_h98(eqt, cp1d, cs; subtract_radiation_losses)
end

@compat public tau_e_h98
push!(document[Symbol("Physics profiles")], :tau_e_h98)

"""
    tau_e_ds03(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.core_sources; subtract_radiation_losses::Bool=true)

Petty's 2003 confinement time scaling

NOTE: Petty uses elongation at the separatrix and makes no distinction between volume and line-average density
"""
function tau_e_ds03(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.core_sources; subtract_radiation_losses::Bool=true)
    dd = top_dd(cp1d)
    total_source = total_sources(cs, cp1d; time0=dd.global_time, fields=Symbol[:power_inside, :total_ion_power_inside])
    total_power_inside = total_source.electrons.power_inside[end] + total_source.total_ion_power_inside[end]
    if subtract_radiation_losses
        total_power_inside -= radiation_losses(cs)
    end
    total_power_inside = max(0.0, total_power_inside)

    isotope_factor =
        trapz(cp1d.grid.volume, sum(ion.density .* ion.element[1].a for ion in cp1d.ion if ion.element[1].z_n == 1.0)) /
        trapz(cp1d.grid.volume, sum(ion.density for ion in cp1d.ion if ion.element[1].z_n == 1.0))

    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

    ne_line = geometric_midplane_line_averaged_density(eqt, cp1d)
    ne_vol = ne_vol_avg(cp1d)

    tauds03 = (
        0.028 *
        abs(eqt.global_quantities.ip / 1e6)^0.83 *
        abs(B0)^0.07 *
        (total_power_inside / 1e6)^-0.55 *
        (0.5 * (ne_line + ne_vol) / 1e19)^0.49 *
        isotope_factor^0.14 *
        R0^2.11 *
        (R0 / eqt.boundary.minor_radius)^-0.30 *
        eqt.boundary.elongation^0.75
    )

    return tauds03
end

"""
    tau_e_ds03(dd::IMAS.dd; time0::Float64=dd.global_time, subtract_radiation_losses::Bool=true)
"""
function tau_e_ds03(dd::IMAS.dd; time0::Float64=dd.global_time, subtract_radiation_losses::Bool=true)
    eqt = dd.equilibrium.time_slice[time0]
    cp1d = dd.core_profiles.profiles_1d[time0]
    cs = dd.core_sources
    return tau_e_ds03(eqt, cp1d, cs; subtract_radiation_losses)
end

@compat public tau_e_ds03
push!(document[Symbol("Physics profiles")], :tau_e_ds03)

"""
    bunit(eqt1d::IMAS.equilibrium__time_slice___profiles_1d)

Calculate bunit from equilibrium
"""
function bunit(eqt1d::IMAS.equilibrium__time_slice___profiles_1d)
    rmin = 0.5 .* (eqt1d.r_outboard .- eqt1d.r_inboard)
    phi = eqt1d.phi
    return gradient(2π * rmin, phi) ./ rmin
end

function bunit(eqt::IMAS.equilibrium__time_slice)
    return bunit(eqt.profiles_1d)
end

@compat public bunit
push!(document[Symbol("Physics profiles")], :bunit)

"""
    greenwald_density(eqt::IMAS.equilibrium__time_slice)

Simple greenwald line-averaged density limit
"""
function greenwald_density(eqt::IMAS.equilibrium__time_slice)
    return greenwald_density(eqt.global_quantities.ip, eqt.boundary.minor_radius)
end

"""
    greenwald_density(ip::T, minor_radius::T) where {T<:Real}
"""
function greenwald_density(ip::T, minor_radius::T) where {T<:Real}
    return abs(ip / 1e6) / (pi * minor_radius^2) * 1e20
end

"""
    greenwald_density(dd::IMAS.dd)
"""
function greenwald_density(dd::IMAS.dd)
    return greenwald_density(dd.equilibrium.time_slice[])
end

"""
    greenwald_density(ps::IMAS.pulse_schedule; time0=global_time(ps))
"""
function greenwald_density(ps::IMAS.pulse_schedule; time0=global_time(ps))
    ip = get_time_array(ps.flux_control.i_plasma, :reference, time0, :linear)
    minor_radius = get_time_array(ps.position_control.minor_radius, :reference, time0, :linear)
    return greenwald_density(ip, minor_radius)
end

@compat public greenwald_density
push!(document[Symbol("Physics profiles")], :greenwald_density)

"""
    greenwald_fraction(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Greewald fraction
"""
function greenwald_fraction(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    nel = IMAS.geometric_midplane_line_averaged_density(eqt, cp1d)
    ngw = greenwald_density(eqt)
    return nel / ngw
end

"""
    greenwald_fraction(dd::IMAS.dd)
"""
function greenwald_fraction(dd::IMAS.dd)
    return greenwald_fraction(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])
end

@compat public greenwald_fraction
push!(document[Symbol("Physics profiles")], :greenwald_fraction)

"""
    geometric_midplane_line_averaged_density(::Nothing, cp1d::IMAS.core_profiles__profiles_1d)

Calculates the averaged density along rho_tor_norm (to be used when equilibrium information is not available)
"""
function geometric_midplane_line_averaged_density(::Nothing, cp1d::IMAS.core_profiles__profiles_1d)
    return trapz(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)
end

"""
    geometric_midplane_line_averaged_density(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Calculates the line averaged density from a midplane horizantal line
"""
function geometric_midplane_line_averaged_density(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    return geometric_midplane_line_averaged_density(eqt, cp1d.electrons.density, cp1d.grid.rho_tor_norm)
end

"""
    geometric_midplane_line_averaged_density(eqt::IMAS.equilibrium__time_slice, ne_profile::AbstractVector{<:Real}, rho_ne::AbstractVector{<:Real})
"""
function geometric_midplane_line_averaged_density(eqt::IMAS.equilibrium__time_slice, ne_profile::AbstractVector{<:Real}, rho_ne::AbstractVector{<:Real})
    a_cp = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard).(rho_ne)
    return trapz(a_cp, ne_profile) / a_cp[end]
end

"""
    geometric_midplane_line_averaged_density(dd::IMAS.dd)
"""
function geometric_midplane_line_averaged_density(dd::IMAS.dd)
    return geometric_midplane_line_averaged_density(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])
end

@compat public greenwald_fraction
push!(document[Symbol("Physics profiles")], :greenwald_fraction)

"""
    ne_line(ps::IMAS.pulse_schedule; time0=global_time(ps))

returns n_e_line from pulse_schedule looking first in `pulse_schedule.density_control.ne_line.reference` and then `pulse_schedule.density_control.greenwald_fraction.reference`
"""
function ne_line(ps::IMAS.pulse_schedule; time0=global_time(ps))
    if !ismissing(ps.density_control.n_e_line, :reference)
        return get_time_array(ps.density_control.n_e_line, :reference, time0, :linear)
    elseif !ismissing(ps.density_control.n_e_greenwald_fraction, :reference)
        return get_time_array(ps.density_control.n_e_greenwald_fraction, :reference, time0, :linear) * greenwald_density(ps; time0)
    else
        error("neither `pulse_schedule.density_control.ne_line.reference` or `pulse_schedule.density_control.greenwald_fraction.reference` have data")
    end
end

@compat public ne_line
push!(document[Symbol("Physics profiles")], :ne_line)

"""
    ne_vol_avg(cp1d::IMAS.core_profiles__profiles_1d)

Volume averaged electron density
"""
function ne_vol_avg(cp1d::IMAS.core_profiles__profiles_1d)
    return trapz(cp1d.grid.volume, cp1d.electrons.density) / cp1d.grid.volume[end]
end

@compat public ne_vol_avg
push!(document[Symbol("Physics profiles")], :ne_vol_avg)

"""
    beta_tor(pressure_average::Real, Bt::Real)

Calculates Beta_tor from pressure and Bt
"""
function beta_tor(pressure_average::Real, Bt::Real)
    return pi * 8.0e-7 * pressure_average / Bt^2
end

@compat public beta_tor
push!(document[Symbol("Physics profiles")], :beta_tor)

"""
    beta_n(beta_tor::Real, minor_radius::Real, Bt::Real, Ip::Real)

Calculates BetaN from beta_tor
"""
function beta_n(beta_tor::Real, minor_radius::Real, Bt::Real, Ip::Real)
    return beta_tor * minor_radius * abs(Bt) / abs(Ip / 1e6) * 1.0e2 # [%]
end

@compat public beta_n
push!(document[Symbol("Physics profiles")], :beta_n)

"""
    pressure_avg_from_beta_n(beta_n::Real, minor_radius::Real, Bt::Real, Ip::Real)

Calculates average pressure from BetaN
"""
function pressure_avg_from_beta_n(beta_n::Real, minor_radius::Real, Bt::Real, Ip::Real)
    return beta_n * abs(Bt) * abs(Ip / 1e6) / (minor_radius * pi * 8.0e-7 * 1.0e2)
end

@compat public pressure_avg_from_beta_n
push!(document[Symbol("Physics profiles")], :pressure_avg_from_beta_n)

"""
    Hmode_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, width::Real; offset::Real=0.0)

Generate H-mode density and temperature profiles evenly spaced in the radial coordinate

  - `edge`: separatrix value
  - `ped`: pedestal value
  - `core`: on-axis value
  - `ngrid`: number of radial grid points
  - `expin`: inner core exponent for H-mode pedestal profile
  - `expout`: outer core exponent for H-mode pedestal profile
  - `width`: full width of pedestal (from the separatrix to the pedestal itself)
  - `offset`: offset of the pedestal center
"""
function Hmode_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, width::Real; offset::Real=0.0)
    @assert edge >= 0.0 "invalid edge = $edge"
    @assert ped >= 0.0 "invalid ped = $ped"
    @assert core >= 0.0 "invalid core = $core"
    @assert expin >= 0.0 "invalid expin = $expin"
    @assert expout >= 0.0 "invalid expout = $expout"
    @assert 0.0 < width < 1.0 "invalid width = $width"

    xpsi = range(0.0, 1.0, ngrid)

    widthp = 0.5 * width  # width as defined in eped
    xphalf = 1.0 - widthp - offset
    pconst = 1.0 - tanh((1.0 - xphalf) / widthp)
    a_t = 2.0 * (ped - edge) / (1.0 + tanh(1.0) - pconst)
    coretanh = 0.5 * a_t * (1.0 - tanh(-xphalf / widthp) - pconst) + edge

    # edge tanh part
    val = @. 0.5 * a_t * (1.0 - tanh((xpsi - xphalf) / widthp) - pconst) + edge + core * 0.0

    # core tanh+polynomial part
    xped = xphalf - widthp
    xtoped = xpsi ./ xped
    for i in 1:ngrid
        if xtoped[i] < 1.0
            @inbounds val[i] += (core - coretanh) * (1.0 - xtoped[i]^expin)^expout
        end
    end

    return val
end

"""
    Hmode_profiles(edge::Real, ped::Real, ngrid::Int, expin::Real, expout::Real, width::Real)

NOTE: The core value is allowed to float
"""
function Hmode_profiles(edge::Real, ped::Real, ngrid::Int, expin::Real, expout::Real, width::Real)
    @assert edge >= 0.0
    @assert ped >= 0.0
    @assert expin >= 0.0
    @assert expout >= 0.0
    @assert 0.0 < width < 1.0 "pedestal width cannot be $width"

    xpsi = range(0.0, 1.0, ngrid)

    widthp = 0.5 * width  # width as defined in eped
    xphalf = 1.0 - widthp
    pconst = 1.0 - tanh((1.0 - xphalf) / widthp)
    a_t = 2.0 * (ped - edge) / (1.0 + tanh(1.0) - pconst)

    # edge tanh part
    val = @. 0.5 * a_t * (1.0 - tanh((xpsi - xphalf) / widthp) - pconst) + edge

    # core tanh+polynomial part
    xped = xphalf - widthp
    xtoped = xpsi ./ xped
    factor = 0.5 * a_t * (1.0 - tanh((xped - xphalf) / widthp) - pconst) + edge
    index = xtoped .< 1.0
    val[index] += factor .* (1.0 .- xtoped[index] .^ expin) .^ expout

    return val
end

@compat public Hmode_profiles
push!(document[Symbol("Physics profiles")], :Hmode_profiles)

"""
    Lmode_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, width::Real)

Generate L-mode density and temperature profiles evenly spaced in the radial coordinate

  - `edge`: separatrix value
  - `ped`: pedestal value
  - `ngrid`: number of radial grid points
  - `expin`: inner core exponent for H-mode pedestal profile
  - `expout`: outer core exponent for H-mode pedestal profile
  - `width`: width of pedestal
"""
function Lmode_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, width::Real)
    rho = range(0.0, 1.0, ngrid)
    rho_ped = 1.0 - width
    rho_ped_idx = argmin(abs.(rho .- rho_ped))

    f(x) = abs.(1.0 - x .^ expin) .^ expout
    profile = f.(rho ./ rho_ped) .* (core - ped) .+ ped
    profile[rho_ped_idx:end] .= range(ped, edge, ngrid - rho_ped_idx + 1)

    res = Optim.optimize(α -> cost_WPED_α!(rho, profile, α, ped, rho_ped), -500, 500, Optim.GoldenSection(); rel_tol=1E-3)
    cost_WPED_α!(rho, profile, res.minimizer, ped, rho_ped)

    return profile
end

@compat public Lmode_profiles
push!(document[Symbol("Physics profiles")], :Lmode_profiles)

"""
    A_effective(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}

A_effective towards L to H scaling see G. Birkenmeier et al 2022 Nucl. Fusion 62 086005
"""
function A_effective(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}
    numerator = T[]
    denominator = T[]
    for ion in cp1d.ion
        if ion.element[1].z_n == 1
            n_int = trapz(cp1d.grid.volume, ion.density)
            push!(numerator, n_int * ion.element[1].a)
            push!(denominator, n_int)
        end
    end
    return reduce(+, numerator) / reduce(+, denominator)
end

@compat public A_effective
push!(document[Symbol("Physics profiles")], :A_effective)

"""
    scaling_L_to_H_power(A_effective::Real, ne_volume::Real, B0::Real, surface_area::Real)

L to H transition power scaling for metal walls and isotope effect according to : G. Birkenmeier et al 2022 Nucl. Fusion 62 086005

inputs in SI and returns power in W
"""
function scaling_L_to_H_power(A_effective::Real, ne_volume::Real, B0::Real, surface_area::Real)
    return 1e6 * 0.8 * 2.0 / A_effective * 0.049 * (ne_volume / 1e20)^0.72 * abs(B0)^0.8 * surface_area^0.94
end

"""
    scaling_L_to_H_power(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
"""
function scaling_L_to_H_power(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    B0 = B0_geo(eqt)
    return scaling_L_to_H_power(
        A_effective(cp1d),
        trapz(cp1d.grid.volume, cp1d.electrons.density) / cp1d.grid.volume[end],
        B0,
        eqt.profiles_1d.surface[end]
    )
end

"""
    scaling_L_to_H_power(dd::IMAS.dd)
"""
function scaling_L_to_H_power(dd::IMAS.dd)
    return scaling_L_to_H_power(dd.core_profiles.profiles_1d[], dd.equilibrium.time_slice[])
end

@compat public scaling_L_to_H_power
push!(document[Symbol("Physics profiles")], :scaling_L_to_H_power)

"""
    L_H_threshold(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Returns ratio of Psol to Plh
"""
function L_H_threshold(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    Psol = power_sol(cs, cp1d)
    Plh = scaling_L_to_H_power(cp1d, eqt)
    return Psol / Plh
end

"""
    L_H_threshold(dd::IMAS.dd)
"""
function L_H_threshold(dd::IMAS.dd)
    return L_H_threshold(dd.core_sources, dd.core_profiles.profiles_1d[], dd.equilibrium.time_slice[])
end

@compat public L_H_threshold
push!(document[Symbol("Physics profiles")], :L_H_threshold)

"""
    satisfies_h_mode_conditions(dd::IMAS.dd)

Returns `true` if the plasma is diverted, has positive triangularity, and `Psol>Plh`
"""
function satisfies_h_mode_conditions(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    diverted = length(eqt.boundary.x_point) > 0
    Psol_gt_Plh = IMAS.L_H_threshold(dd) > getproperty(dd.requirements, :lh_power_threshold_fraction, 1.0)
    positive_triangularity = eqt.boundary.triangularity > 0.0
    if Psol_gt_Plh && diverted && positive_triangularity
        return true
    else
        return false
    end
end

@compat public satisfies_h_mode_conditions
push!(document[Symbol("Physics profiles")], :satisfies_h_mode_conditions)

"""
    ITB(rho0::T, width::T, height::T, rho::AbstractVector{T}) where {T<:Real}

tanh profile to be added to existing profiles to model Internal Transport Barrier (ITB).
The ITB is centered at `rho0`, with a full width of `width` and given `height`
"""
function ITB(rho0::T, width::T, height::T, rho::AbstractVector{T}) where {T<:Real}
    return @. (tanh((-rho + rho0) / width * pi) + 1.0) * 0.5 * height
end

@compat public ITB
push!(document[Symbol("Physics profiles")], :ITB)

"""
    ITB_profile(rho::AbstractVector{T}, input_profile::AbstractVector{T}, rho0::T, width::T, height_ratio::T) where {T<:Real}

Add ITB to existing profile. The ITB is centered at `rho0`, with a full width of `width`, and a height expressed as a ratio of the input_profile evaluated on axis
"""
function ITB_profile(rho::AbstractVector{T}, input_profile::AbstractVector{T}, rho0::T, width::T, height_ratio::T) where {T<:Real}
    height = interp1d(rho, input_profile).(0.0) * height_ratio
    itb = ITB(rho0, width, height, rho)
    return input_profile .+ itb
end

@compat public ITB_profile
push!(document[Symbol("Physics profiles")], :ITB_profile)

"""
    species(cp1d::IMAS.core_profiles__profiles_1d; only_electrons_ions::Symbol=:all, only_thermal_fast::Symbol=:all)

Returns species index and names (followed by "_fast" if density_fast is present), for example:

    (0, :electrons)
    (0, :electrons_fast)
    (1, :DT)
    (2, :Kr83)
    (3, :He4)
    (1, :DT_fast)
    (3, :He4_fast)
"""
function species(cp1d::IMAS.core_profiles__profiles_1d; only_electrons_ions::Symbol=:all, only_thermal_fast::Symbol=:all)
    @assert only_electrons_ions ∈ (:all, :electrons, :ions) "only_electrons_ions can be one of (:all, :electrons, :ions)"
    @assert only_thermal_fast ∈ (:all, :thermal, :fast) "only_thermal_fast can be one of (:all, :thermal, :fast)"
    out = []
    if only_electrons_ions ∈ (:all, :electrons)
        if only_thermal_fast ∈ (:all, :thermal) && sum(cp1d.electrons.density_thermal) > 0.0
            push!(out, (0, :electrons))
        end
        if only_thermal_fast ∈ (:all, :fast) && sum(cp1d.electrons.density_fast) > 0.0
            push!(out, (0, :electrons_fast))
        end
    end
    if only_electrons_ions ∈ (:all, :ions)
        if only_thermal_fast ∈ (:all, :thermal)
            dd_thermal = ((k, Symbol(ion.label)) for (k, ion) in enumerate(cp1d.ion) if sum(ion.density_thermal) > 0.0)
            for item in dd_thermal
                push!(out, item)
            end
        end
        if only_thermal_fast ∈ (:all, :fast)
            dd_fast = ((k, Symbol("$(ion.label)_fast")) for (k, ion) in enumerate(cp1d.ion) if sum(ion.density_fast) > 0.0)
            for item in dd_fast
                push!(out, item)
            end
        end
    end
    return out
end

@compat public species
push!(document[Symbol("Physics profiles")], :species)

"""
    is_quasi_neutral(cp1d::IMAS.core_profiles__profiles_1d; rtol::Float64=0.001)

Returns true if quasi neutrality is satisfied within a relative tolerance
"""
function is_quasi_neutral(cp1d::IMAS.core_profiles__profiles_1d; rtol::Float64=0.001)
    Nis = sum(sum(ion.density .* ion.z_ion for ion in cp1d.ion))
    Ne = sum(cp1d.electrons.density)
    if (1.0 + rtol) >= (Ne / Nis) >= (1.0 - rtol)
        return true
    else
        return false
    end
end

"""
    is_quasi_neutral(dd::IMAS.dd)
"""
function is_quasi_neutral(dd::IMAS.dd; rtol::Float64=0.001)
    return is_quasi_neutral(dd.core_profiles.profiles_1d[]; rtol)
end

@compat public is_quasi_neutral
push!(document[Symbol("Physics profiles")], :is_quasi_neutral)

"""
    enforce_quasi_neutrality!(cp1d::IMAS.core_profiles__profiles_1d, species::Symbol)

If `species` is `:electrons` then updates `electrons.density_thermal` to meet quasi neutrality condtion.

If `species` is a ion species, it evaluates the difference in number of charges needed to reach quasineutrality, and assigns positive difference to target ion density_thermal species and negative difference to electrons density_thermal

Also, sets `density` to the original expression
"""
function enforce_quasi_neutrality!(cp1d::IMAS.core_profiles__profiles_1d, species::Symbol)
    # Make sure expressions are used for total densities
    empty!(cp1d.electrons, :density)
    for ion in cp1d.ion
        empty!(ion, :density)
    end

    if species == :electrons
        cp1d.electrons.density_thermal = sum(ion.density .* ion.z_ion for ion in cp1d.ion) .- cp1d.electrons.density_fast

    else
        # identify ion species
        species_indx = findfirst(Symbol(ion.label) == species for ion in cp1d.ion)
        @assert species_indx !== nothing
        ion0 = cp1d.ion[species_indx]

        # evaluate the difference in number of thermal charges needed to reach quasineutrality
        ne = cp1d.electrons.density .+ cp1d.electrons.density_fast
        q_density_difference = ne .- sum(ion.density .* ion.z_ion for ion in cp1d.ion if ion !== ion0) .- ion0.density_fast .* ion0.z_ion

        # positive difference is assigned to target ion density_thermal
        index = q_density_difference .> 0.0
        ion0.density_thermal = zero(cp1d.electrons.density_thermal)
        ion0.density_thermal[index] .= q_density_difference[index] ./ ion0.z_ion

        # negative difference is assigned to electrons density_thermal
        index = q_density_difference .< 0.0
        cp1d.electrons.density_thermal[index] .+= abs.(q_density_difference[index])
    end

    return nothing
end

"""
    enforce_quasi_neutrality!(dd::IMAS.dd, species::Symbol)
"""
function enforce_quasi_neutrality!(dd::IMAS.dd, species::Symbol)
    return enforce_quasi_neutrality!(dd.core_profiles.profiles_1d[], species)
end

@compat public enforce_quasi_neutrality!
push!(document[Symbol("Physics profiles")], :enforce_quasi_neutrality!)

"""
    lump_ions_as_bulk_and_impurity(ions::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion}, rho_tor_norm::Vector{<:Real})

Changes core_profiles.ion to 2 species, one bulk species (H, D, T) and one combined impurity species
"""
function lump_ions_as_bulk_and_impurity(cp1d::IMAS.core_profiles__profiles_1d{T})::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion{T}} where {T<:Real}

    ions = cp1d.ion

    if length(ions) < 2
        error("lump_ions_as_bulk_and_impurity() requires at least two ion species")
    end

    zs = [ion.element[1].z_n for ion in ions]
    as = [ion.element[1].a for ion in ions]
    bulk_index = findall(zs .== 1)
    impu_index = findall(zs .!= 1)

    rho_tor_norm = cp1d.grid.rho_tor_norm
    ratios = zeros(length(rho_tor_norm), length(ions))
    ntot = zeros(length(rho_tor_norm))
    for index in (bulk_index, impu_index)
        ntot .*= 0.0
        for ix in index
            tmp = ions[ix].density_thermal
            ratios[:, ix] = tmp
            ntot .+= tmp
        end
        for ix in index
            ratios[:, ix] ./= ntot
        end
        @assert all(ntot .> 0.0) "Species $([ions[ix].label for ix in index]) have zero density: $(ntot)"
    end

    ne = cp1d.electrons.density_thermal
    n1 = zero(ne)
    for ion in cp1d.ion
        if ion.element[1].z_n == 1
            n1 .+= ion.density_thermal
        end
    end
    Zi = (n1 .- ne .* cp1d.zeff) ./ (n1 .- ne)
    Zi = sum(Zi) / length(Zi) # make Zi constant as function of radius
    ni = (ne .- n1) ./ Zi

    ions2 = IMAS.IDSvector{IMAS.core_profiles__profiles_1d___ion{T}}()

    # bulk ions
    push!(ions2, IMAS.core_profiles__profiles_1d___ion{T}())
    bulk = ions2[end]
    resize!(bulk.element, 1)
    bulk.label = "bulk"
    bulk.element[1].z_n = 1.0
    IMAS.setproperty!(bulk, :density_thermal, n1; error_on_missing_coordinates=false)

    # impurity ions
    push!(ions2, IMAS.core_profiles__profiles_1d___ion{T}())
    impu = ions2[end]
    resize!(impu.element, 1)
    impu.label = "impurity"
    impu.element[1].z_n = Zi
    IMAS.setproperty!(impu, :density_thermal, ni; error_on_missing_coordinates=false)

    # weight different ion quantities based on their density
    for (index, ion2) in ((bulk_index, bulk), (impu_index, impu))
        ion2.element[1].a = 0.0
        for ix in index
            ion2.element[1].a += as[ix] * sum(ratios[:, ix]) / length(ratios[:, ix])
        end
        for item in (:temperature, :rotation_frequency_tor)
            value = rho_tor_norm .* 0.0
            IMAS.setproperty!(ion2, item, value; error_on_missing_coordinates=false)
            for ix in index
                tmp = getproperty(ions[ix], item, T[])
                if !isempty(tmp)
                    value .+= tmp .* ratios[:, ix]
                end
            end
        end
    end

    return ions2
end

@compat public lump_ions_as_bulk_and_impurity
push!(document[Symbol("Physics profiles")], :lump_ions_as_bulk_and_impurity)

"""
    zeff(cp1d::IMAS.core_profiles__profiles_1d; temperature_dependent_ionization_state::Bool=true)

Returns plasma effective charge

`temperature_dependent_ionization_state` evaluates Zeff with average ionization state of an ion at a given temperature
"""
function zeff(cp1d::IMAS.core_profiles__profiles_1d; temperature_dependent_ionization_state::Bool=true)
    z = zero(cp1d.grid.rho_tor_norm)
    for ion in cp1d.ion
        if temperature_dependent_ionization_state
            Zi = avgZ(ion.element[1].z_n, ion.temperature)
        else
            Zi = ion.element[1].z_n
        end
        z .+= ion.density .* Zi .^ 2
    end
    ne = cp1d.electrons.density
    for k in eachindex(z)
        z[k] = max(1.0, z[k] / ne[k])
    end
    return z
end

@compat public zeff
push!(document[Symbol("Physics profiles")], :zeff)

Memoize.@memoize function avgZinterpolator(filename::String)
    txt = open(filename, "r") do io
        return strip(read(io, String))
    end

    iion = Vector{Int}()
    Ti = Vector{Float64}()
    for (k, line) in enumerate(split(txt, "\n"))
        if k == 1
            continue
        elseif k == 2
            continue
        elseif k == 3
            append!(iion, collect(map(x -> parse(Int, x), split(line))))
        elseif k == 4
            append!(Ti, collect(map(x -> parse(Float64, x), split(line))))
        end
    end

    data = zeros(length(Ti), length(iion))
    for (k, line) in enumerate(split(txt, "\n"))
        if k > 4
            data[k-4, :] = map(x -> parse(Float64, x), split(line))
        end
    end

    return Interpolations.extrapolate(Interpolations.interpolate((log10.(Ti), iion), log10.(data .+ 1.0), Interpolations.Gridded(Interpolations.Linear())), Interpolations.Flat())
end

"""
    avgZ(Z::Float64,Ti::T)::T

Returns average ionization state of an ion at a given temperature
"""
function avgZ(Z::Float64, Ti::T)::T where {T}
    func = avgZinterpolator(joinpath(@__DIR__, "..", "..", "data", "Zavg_z_t.dat"))
    return 10.0 .^ (func.(log10.(Ti ./ 1E3), Z)) .- 1.0
end

@compat public avgZ
push!(document[Symbol("Physics profiles")], :avgZ)

"""
    t_i_average(cp1d::IMAS.core_profiles__profiles_1d)::Vector{<:Real}

Returns the average ion temperature weighted by each species density over the total number of ion particles
"""
function t_i_average(cp1d::IMAS.core_profiles__profiles_1d)::Vector{<:Real}
    n = length(cp1d.grid.rho_tor_norm)
    t_i_a = zeros(n)
    ntot = zeros(n)
    for ion in cp1d.ion
        t_i_a += ion.density_thermal .* ion.temperature
        ntot += ion.density_thermal
    end
    return t_i_a ./= ntot
end

@compat public t_i_average
push!(document[Symbol("Physics profiles")], :t_i_average)

"""
    exponential_profile(x::AbstractArray{<:Real}, x0::Real, T0::Real, T1::Real, alpha::Real)

Function for edge blending using exponential function
"""
function exponential_profile(x::AbstractArray{<:Real}, x0::Real, T0::Real, T1::Real, alpha::Real)
    @assert x[1] == 0.0
    @assert x[end] == 1.0
    @assert 0.0 < x0 < 1.0
    if alpha == 0.0
        sigma = 1E10
    else
        sigma = 1 / alpha
    end
    y = exponential_profile(x, x0, sigma)
    y0 = exponential_profile(x0, x0, sigma)
    y1 = exponential_profile(1.0, x0, sigma)
    return @. (y - y1) / (y0 - y1) * (T0 - T1) + T1
end

function exponential_profile(x::Union{Real,AbstractArray}, x0::Real, sigma::Real)
    return exp.(.-(x .- x0) / sigma)
end

@compat public edge_profile
push!(document[Symbol("Physics profiles")], :edge_profile)

"""
    core_edge_energy(cp1d::IMAS.core_profiles__profiles_1d, rho_ped::Real; thermal::Bool=true)

Evaluate stored energy in the core (< rho_ped) or the pedestal (> rho_ped)
"""
function core_edge_energy(cp1d::IMAS.core_profiles__profiles_1d, rho_ped::Real; thermal::Bool=true)
    if thermal
        p = pressure_thermal(cp1d)
    else
        p = pressure(cp1d)
    end
    rho_tor_norm = cp1d.grid.rho_tor_norm
    rho_bound_idx = argmin(abs(rho - rho_ped) for rho in rho_tor_norm)
    pedge = p[rho_bound_idx]
    fedge = (k, x) -> (k <= rho_bound_idx) ? pedge : p[k]
    fcore = (k, x) -> (k <= rho_bound_idx) ? (p[k] - pedge) : 0.0
    volume = cp1d.grid.volume
    core_value = 1.5 * trapz(volume, fcore)
    edge_value = 1.5 * trapz(volume, fedge)
    return (core_value=core_value, edge_value=edge_value)
end

@compat public core_edge_energy
push!(document[Symbol("Physics profiles")], :core_edge_energy)

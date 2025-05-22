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
        if hasdata(ion, :density_thermal)
            p .+= ion.temperature .* ion.density_thermal
        end
    end
    return p .* mks.e
end

@compat public pressure_thermal
push!(document[Symbol("Physics profiles")], :pressure_thermal)

"""
    beta_tor_thermal_norm(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Normalized toroidal beta from thermal pressure only, defined as 100 * beta_tor_thermal * a[m] * B0 [T] / ip [MA]
"""
function beta_tor_thermal_norm(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    return beta_tor(eqt, cp1d; norm=true, thermal=true)
end

@compat public beta_tor_thermal_norm
push!(document[Symbol("Physics profiles")], :beta_tor_thermal_norm)

"""
    beta_tor_norm(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Normalized toroidal beta from total pressure, defined as 100 * beta_tor * a[m] * B0 [T] / ip [MA]
"""
function beta_tor_norm(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    return beta_tor(eqt, cp1d; norm=true, thermal=false)
end

@compat public beta_tor_norm
push!(document[Symbol("Physics profiles")], :beta_tor_norm)

"""
    beta_tor(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; norm::Bool, thermal::Bool)

Toroidal beta, defined as the volume-averaged total perpendicular pressure divided by (B0^2/(2*mu0)), i.e. beta_toroidal = 2 mu0 int(p dV) / V / B0^2
"""
function beta_tor(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; norm::Bool, thermal::Bool)
    if thermal
        pressure = cp1d.pressure_thermal
    else
        pressure = cp1d.pressure
    end
    @assert !any(isnan.(pressure))
    volume = cp1d.grid.volume
    pressure_avg = trapz(volume, pressure) / volume[end]

    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    ip = eqt.global_quantities.ip

    beta_tor = 2.0 * mks.μ_0 * pressure_avg / B0^2

    if norm
        out = beta_tor * eqt.boundary.minor_radius * abs(B0) / abs(ip / 1e6) * 1.0e2
    else
        out = beta_tor
    end

    return out
end

@compat public beta_tor
push!(document[Symbol("Physics profiles")], :beta_tor)

"""
    ions_sort_map::Dict

Dictionary that returns an integer number, used to sort ions by their Z
with fast-ion species listed at the end
"""
const ions_sort_map = Dict(
    sym => z for (z, sym) in enumerate([
        "H", "D", "DT", "T", "He", "α ",
        "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ])
)
for element in collect(keys(ions_sort_map))
    ions_sort_map["$(element)_fast"] = 1000 + ions_sort_map[element]
end

function list_ions!(ct::IMAS.core_transport, ions::Set{Symbol}; time0::Float64)
    for model in ct.model
        if isempty(model.profiles_1d)
            continue
        end
        time = [ids.time for ids in model.profiles_1d]
        index = nearest_causal_time(time, time0; bounds_error=false).index
        ct1d = model.profiles_1d[index]
        for ion in ct1d.ion
            push!(ions, Symbol(ion.label))
        end
    end
    return ions
end

function list_ions!(cp::IMAS.core_profiles, ions::Set{Symbol}; time0::Float64)
    if isempty(cp.profiles_1d)
        return ions
    end
    time = [ids.time for ids in cp.profiles_1d]
    index = nearest_causal_time(time, time0; bounds_error=false).index
    cp1d = cp.profiles_1d[index]
    for ion in cp1d.ion
        push!(ions, Symbol(ion.label))
    end
    return ions
end

function list_ions!(cs::IMAS.core_sources, ions::Set{Symbol}; time0::Float64)
    for source in cs.source
        if isempty(source.profiles_1d)
            continue
        end
        time = [ids.time for ids in source.profiles_1d]
        index = nearest_causal_time(time, time0; bounds_error=false).index
        sc1d = source.profiles_1d[index]
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
    ion_set = Set{Symbol}()
    for ids in (ids1, idss...)
        list_ions!(ids, ion_set; time0)
    end
    ions_sorted = sort!(collect(ion_set); by=x -> ions_sort_map[replace(string(x), r"[0-9]+" => "")])
    return ions_sorted
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
    rho_index = argmin_abs(cp1d.grid.rho_tor_norm, su.local.pedestal.position.rho_tor_norm)
    pressure = cp1d.pressure_thermal
    pressure[1:rho_index] .= cp1d.pressure_thermal[rho_index]
    return 3.0 / 2.0 * trapz(cp1d.grid.volume, pressure)
end

@compat public energy_thermal_ped
push!(document[Symbol("Physics profiles")], :energy_thermal_ped)

"""
    tau_e_thermal(cp1d::IMAS.core_profiles__profiles_1d, sources::IMAS.core_sources; , include_radiation::Bool=true, include_time_derivative::Bool=true

Evaluate thermal energy confinement time

NOTE: This can go to infinity if there's more power coming out of the plasma than there is going in
"""
function tau_e_thermal(dd::IMAS.dd; time0::Float64=dd.global_time, include_radiation::Bool=true, include_time_derivative::Bool=true)
    cp1d = dd.core_profiles.profiles_1d[time0]
    tot_pow_in = total_power_inside(dd.core_sources, cp1d; time0, include_radiation, include_time_derivative)
    tot_pow_in = max(0.0, tot_pow_in)
    return energy_thermal(cp1d) / tot_pow_in
end


@compat public tau_e_thermal
push!(document[Symbol("Physics profiles")], :tau_e_thermal)

"""
    tau_e_h98(dd::IMAS.dd; time0::Float64=dd.global_time, include_radiation::Bool=true, include_time_derivative::Bool=true)

H98y2 ITER elmy H-mode confinement time scaling

NOTE: H98y2 uses aereal elongation

NOTE: (IPB_Chap_2.pdf, pg. 72) ...for practical reasons, the power lost by radiation inside
the separatrix of the existing devices has been neglected when deriving the scalings.
However, for ITER, such radiation is subtracted from the loss power when calculating the
projected energy confinement time.

See Table 5 in https://iopscience.iop.org/article/10.1088/0029-5515/39/12/302/pdf and https://iopscience.iop.org/article/10.1088/0029-5515/48/9/099801/pdf for additional correction with plasma_volume
"""
function tau_e_h98(dd::IMAS.dd; time0::Float64=dd.global_time, include_radiation::Bool=true, include_time_derivative::Bool=true)
    eqt = dd.equilibrium.time_slice[time0]
    cp1d = dd.core_profiles.profiles_1d[time0]
    cs = dd.core_sources

    tot_pow_in = total_power_inside(cs, cp1d; time0, include_radiation, include_time_derivative)
    tot_pow_in = max(0.0, tot_pow_in)

    isotope_factor =
        trapz(cp1d.grid.volume, sum(ion.density_thermal .* ion.element[1].a for ion in cp1d.ion if ion.element[1].z_n == 1.0)) /
        trapz(cp1d.grid.volume, sum(ion.density_thermal for ion in cp1d.ion if ion.element[1].z_n == 1.0))

    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

    κ_areal = areal_elongation(eqt)

    nel = ne_line(eqt, cp1d)

    tau98 = (
        0.0562 *
        abs(eqt.global_quantities.ip / 1e6)^0.93 *
        abs(B0)^0.15 *
        (tot_pow_in / 1e6)^-0.69 *
        (nel / 1e19)^0.41 *
        isotope_factor^0.19 *
        R0^1.97 *
        (R0 / eqt.boundary.minor_radius)^-0.58 *
        κ_areal^0.78
    )
    return tau98
end

@compat public tau_e_h98
push!(document[Symbol("Physics profiles")], :tau_e_h98)

"""
    tau_e_ds03(dd::IMAS.dd; time0::Float64=dd.global_time, include_radiation::Bool=true, include_time_derivative::Bool=true)

Petty's 2003 confinement time scaling

NOTE: Petty uses elongation at the separatrix and makes no distinction between volume and line-average density
"""
function tau_e_ds03(dd::IMAS.dd; time0::Float64=dd.global_time, include_radiation::Bool=true, include_time_derivative::Bool=true)
    eqt = dd.equilibrium.time_slice[time0]
    cp1d = dd.core_profiles.profiles_1d[time0]
    cs = dd.core_sources

    tot_pow_in = total_power_inside(cs, cp1d; time0, include_radiation, include_time_derivative)
    tot_pow_in = max(0.0, tot_pow_in)

    isotope_factor =
        trapz(cp1d.grid.volume, sum(ion.density_thermal .* ion.element[1].a for ion in cp1d.ion if ion.element[1].z_n == 1.0)) /
        trapz(cp1d.grid.volume, sum(ion.density_thermal for ion in cp1d.ion if ion.element[1].z_n == 1.0))

    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

    nel = ne_line(eqt, cp1d)
    ne_vol = ne_vol_avg(cp1d)

    tauds03 = (
        0.028 *
        abs(eqt.global_quantities.ip / 1e6)^0.83 *
        abs(B0)^0.07 *
        (tot_pow_in / 1e6)^-0.55 *
        (0.5 * (nel + ne_vol) / 1e19)^0.49 *
        isotope_factor^0.14 *
        R0^2.11 *
        (R0 / eqt.boundary.minor_radius)^-0.30 *
        eqt.boundary.elongation^0.75
    )

    return tauds03
end

@compat public tau_e_ds03
push!(document[Symbol("Physics profiles")], :tau_e_ds03)

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
    nel = ne_line(eqt, cp1d)
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
    ne_line(::Nothing, cp1d::IMAS.core_profiles__profiles_1d)

Calculates the averaged density along rho_tor_norm (to be used when equilibrium information is not available)
"""
function ne_line(::Nothing, cp1d::IMAS.core_profiles__profiles_1d)
    return trapz(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)
end

"""
    ne_line(eqt::IMAS.equilibrium__time_slice, ne_profile::AbstractVector{<:Real}, rho_ne::AbstractVector{<:Real})
"""
function ne_line(eqt::IMAS.equilibrium__time_slice, ne_profile::AbstractVector{<:Real}, rho_ne::AbstractVector{<:Real})
    a_cp = interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard).(rho_ne)
    return trapz(a_cp, ne_profile) / a_cp[end]
end

"""
    ne_line(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Calculates the line averaged density from the equilibrium midplane horizantal line
"""
function ne_line(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    return ne_line(eqt, cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm)
end

"""
    ne_line(dd::IMAS.dd; time0::Float64=dd.global_time)
"""
function ne_line(dd::IMAS.dd; time0::Float64=dd.global_time)
    return ne_line(dd.equilibrium.time_slice[time0], dd.core_profiles.profiles_1d[time0])
end

"""
    ne_line(ps::IMAS.pulse_schedule; time0=global_time(ps))

returns n_e_line from pulse_schedule looking first in `pulse_schedule.density_control.ne_line.reference` and then `pulse_schedule.density_control.greenwald_fraction.reference`
"""
function ne_line(ps::IMAS.pulse_schedule; time0::Float64=global_time(ps))
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
    return trapz(cp1d.grid.volume, cp1d.electrons.density_thermal) / cp1d.grid.volume[end]
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
    @assert edge >= 0.0 "invalid edge = $edge"
    @assert ped >= 0.0 "invalid ped = $ped"
    @assert expin >= 0.0 "invalid expin = $expin"
    @assert expout >= 0.0 "invalid expout = $expout"
    @assert 0.0 < width < 1.0 "invalid width = $width"

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
    for k in eachindex(xtoped)
        if xtoped[k] < 1.0
            @inbounds val[k] += factor * (1.0 - xtoped[k] ^ expin) ^ expout
        end
    end

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
    rho_ped_idx = argmin_abs(rho, rho_ped)

    f(x) = abs.(1.0 - x .^ expin) .^ expout
    profile = f.(rho ./ rho_ped) .* (core - ped) .+ ped
    profile[rho_ped_idx:end] .= range(ped, edge, ngrid - rho_ped_idx + 1)

    res = Optim.optimize(α -> cost_WPED_α!(rho, profile, α, ped, rho_ped), -500, 500, Optim.Brent(); rel_tol=1E-3)
    cost_WPED_α!(rho, profile, res.minimizer, ped, rho_ped)

    return profile
end

@compat public Lmode_profiles
push!(document[Symbol("Physics profiles")], :Lmode_profiles)

"""
    A_effective(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}

A_effective towards L to H scaling see: G. Birkenmeier et al 2022 Nucl. Fusion 62 086005
"""
function A_effective(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}
    numerator = zero(T)
    denominator = zero(T)
    for ion in cp1d.ion
        if ion.element[1].z_n == 1
            n_int = trapz(cp1d.grid.volume, ion.density_thermal)
            numerator += n_int * ion.element[1].a
            denominator += n_int
        end
    end

    return numerator / denominator
end

@compat public A_effective
push!(document[Symbol("Physics profiles")], :A_effective)

"""
    scaling_L_to_H_power(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice; metallic_wall:Bool=true)

L to H transition power threshold for metal walls and isotope effect according to: G. Birkenmeier et al 2022 Nucl. Fusion 62 086005

See also: `Reducing the L-H Transition Power Threshold in DIII-D ITER-Similar-Shape Hydrogen Plasmas, L. Schmitz IAEA FEC 2020`

NOTE: This scaling does not reflect a well documented nonmonotonic density dependence.

L_H_threshold doubles when ion ∇B drift is away from X-point

returns Infinity if the plasma is diverted or in negative triangularity
"""
function scaling_L_to_H_power(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice; metallic_wall::Bool=true)
    if isempty(eqt.boundary.x_point)
        return Inf
    end
    if eqt.boundary.triangularity < 0.0
        return Inf
    end

    Bgeo = B0_geo(eqt)

    if Bgeo < 0
        # COCOS 11: Negative Bt is clockwise --> B×∇B is downward
        # This checks out with ITER. Negative B in COCOS 11 and lower single null.
        ∇B_drift_direction = -1
    else
        # COCOS 11: Positive Bt is counter clockwise --> B×∇B is upward
        # This checks out with DIII-D shot 200204, where at 1.2 s
        # they go from lower to upper single null to go into H-mode
        ∇B_drift_direction = +1
    end
    # L_H_threshold doubles when ion ∇B drift is away from primary X-point
    if sign(eqt.boundary.x_point[1].z) == ∇B_drift_direction
        ∇B_drift_multiplier = 1.0
    else
        ∇B_drift_multiplier = 2.0
    end

    # The Martin scaling is only valid for plasma densities above the power threshold minimum
    ne_volume = trapz(cp1d.grid.volume, cp1d.electrons.density_thermal) / cp1d.grid.volume[end] / 1E20
    Rgeo = eqt.boundary.geometric_axis.r
    ageo = eqt.boundary.minor_radius
    ne_min = 0.7 * abs(eqt.global_quantities.ip / 1e6)^0.34 * abs(Bgeo)^0.62 * ageo^-0.95 * (Rgeo / ageo)^0.4 # in 1e19 m^-3
    ne_min *= 0.1 # [10^20 m⁻³]
    ne_volume = max(ne_min, ne_volume)

    surface_area = eqt.profiles_1d.surface[end]

    if metallic_wall
        # PLH is at least 20% lower in a metallic wall compared to the original ITPA scaling, which was derived from carbon wall data
        wall_factor = 0.8
    else
        wall_factor = 1.0
    end

    isotope_effect = 2.0 / A_effective(cp1d)

    power_threshold = 1e6 * wall_factor * isotope_effect * 0.049 * ne_volume^0.72 * abs(Bgeo)^0.8 * surface_area^0.94 * ∇B_drift_multiplier

    return power_threshold
end

"""
    scaling_L_to_H_power(dd::IMAS.dd; time0::Float64=dd.global_time)
"""
function scaling_L_to_H_power(dd::IMAS.dd; time0::Float64=dd.global_time)
    return scaling_L_to_H_power(dd.core_profiles.profiles_1d[time0], dd.equilibrium.time_slice[time0])
end

@compat public scaling_L_to_H_power
push!(document[Symbol("Physics profiles")], :scaling_L_to_H_power)

"""
    L_H_threshold(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice; time0::Float64=dd.global_time)

Returns ratio of Psol to Plh
"""
function L_H_threshold(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice; time0::Float64=dd.global_time)
    Psol = power_sol(cs, cp1d; time0)
    Plh = scaling_L_to_H_power(cp1d, eqt)
    return Psol / Plh
end

"""
    L_H_threshold(dd::IMAS.dd)
"""
function L_H_threshold(dd::IMAS.dd; time0::Float64=dd.global_time)
    return L_H_threshold(dd.core_sources, dd.core_profiles.profiles_1d[time0], dd.equilibrium.time_slice[time0]; time0)
end

@compat public L_H_threshold
push!(document[Symbol("Physics profiles")], :L_H_threshold)

"""
    satisfies_h_mode_conditions(dd::IMAS.dd; threshold_multiplier::Float64=1.0)

Returns `true` if the plasma is diverted, has positive triangularity, and `Psol > Plh * threshold_multiplier`
"""
function satisfies_h_mode_conditions(dd::IMAS.dd; threshold_multiplier::Float64=1.0)
    Psol_gt_Plh = L_H_threshold(dd) > threshold_multiplier
    if Psol_gt_Plh
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
    out = @NamedTuple{index::Int64, name::Symbol}[]
    if only_electrons_ions ∈ (:all, :electrons)
        if only_thermal_fast ∈ (:all, :thermal) && hasdata(cp1d.electrons, :density_thermal) && sum(cp1d.electrons.density_thermal) > 0.0
            push!(out, (index=0, name=:electrons))
        end
        if only_thermal_fast ∈ (:all, :fast) && hasdata(cp1d.electrons, :density_fast) && sum(cp1d.electrons.density_fast) > 0.0
            push!(out, (index=0, name=:electrons_fast))
        end
    end
    if only_electrons_ions ∈ (:all, :ions)
        if only_thermal_fast ∈ (:all, :thermal)
            dd_thermal = ((index=k, name=Symbol(ion.label)) for (k, ion) in enumerate(cp1d.ion) if hasdata(ion, :density_thermal) && sum(ion.density_thermal) > 0.0)
            for item in dd_thermal
                push!(out, item)
            end
        end
        if only_thermal_fast ∈ (:all, :fast)
            dd_fast = ((index=k, name=Symbol("$(ion.label)_fast")) for (k, ion) in enumerate(cp1d.ion) if hasdata(ion, :density_fast) && sum(ion.density_fast) > 0.0)
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
    IMAS.unfreeze!(cp1d.electrons, :density)
    for ion in cp1d.ion
        IMAS.unfreeze!(ion, :density)
    end

    if species == :electrons
        cp1d.electrons.density_thermal = sum(ion.density .* ion.z_ion for ion in cp1d.ion) .- cp1d.electrons.density_fast

    else
        # identify ion species
        species_indx = findfirst(Symbol(ion.label) == species for ion in cp1d.ion)
        @assert species_indx !== nothing
        ion0 = cp1d.ion[species_indx]

        # evaluate the difference in number of thermal charges needed to reach quasineutrality
        if hasdata(cp1d.electrons, :density_fast)
            ne = cp1d.electrons.density_thermal .+ cp1d.electrons.density_fast
        else
            ne = cp1d.electrons.density_thermal
        end
        if hasdata(ion0, :density_fast)
            q_density_difference = ne .- sum(ion.density .* ion.z_ion for ion in cp1d.ion if ion !== ion0) .- ion0.density_fast .* ion0.z_ion
        else
            q_density_difference = ne .- sum(ion.density .* ion.z_ion for ion in cp1d.ion if ion !== ion0)
        end

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
    setproperty!(bulk, :density_thermal, n1; error_on_missing_coordinates=false)

    # impurity ions
    push!(ions2, IMAS.core_profiles__profiles_1d___ion{T}())
    impu = ions2[end]
    resize!(impu.element, 1)
    impu.label = "impurity"
    impu.element[1].z_n = Zi
    setproperty!(impu, :density_thermal, ni; error_on_missing_coordinates=false)

    # weight different ion quantities based on their density
    for (index, ion2) in ((bulk_index, bulk), (impu_index, impu))
        ion2.element[1].a = 0.0
        for ix in index
            ion2.element[1].a += as[ix] * sum(ratios[:, ix]) / length(ratios[:, ix])
        end
        for item in (:temperature, :rotation_frequency_tor)
            value = rho_tor_norm .* 0.0
            setproperty!(ion2, item, value; error_on_missing_coordinates=false)
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
    new_impurity_fraction!(cp1d::IMAS.core_profiles__profiles_1d, impurity_name::Symbol, impurity_ne_fraction::Real)

Add thermal impurity to core_profiles given impurity fraction.

The new impurity will have the same density profile shape as the electron density profile and the same ion temperature as the main ion specie.

The main ion specie will be adjusted to have quasineutrality.
"""
function new_impurity_fraction!(cp1d::IMAS.core_profiles__profiles_1d, impurity_name::Symbol, impurity_ne_fraction::Real)
    impurity_label = ion_properties(impurity_name).label
    species_indx = findfirst(ion.label == impurity_label for ion in cp1d.ion)
    if species_indx === nothing
        resize!(cp1d.ion, length(cp1d.ion) + 1)
        new_ion = cp1d.ion[end]
    else
        new_ion = cp1d.ion[species_indx]
        empty!(new_ion)
    end
    IMAS.ion_element!(new_ion, impurity_name)

    new_ion.density_thermal = cp1d.electrons.density_thermal .* impurity_ne_fraction
    new_ion.density_fast = zeros(length(cp1d.electrons.density_thermal))
    new_ion.temperature = cp1d.ion[1].temperature

    main_ion = cp1d.ion[1]
    IMAS.enforce_quasi_neutrality!(cp1d, Symbol(main_ion.label))

    return cp1d
end

@compat public new_impurity_fraction!
push!(document[Symbol("Physics profiles")], :new_impurity_fraction!)

"""
    new_impurity_radiation!(dd::IMAS.dd, impurity_name::Symbol, total_radiated_power::Real)

Add thermal impurity to core_profiles to match a total radiated power.

The new impurity will have the same density profile shape as the electron density profile and the same ion temperature as the main ion specie.

The main ion specie will be adjusted to have quasineutrality.
"""
function new_impurity_radiation!(dd::IMAS.dd, impurity_name::Symbol, total_radiated_power::Real)
    @assert total_radiated_power <= 0.0 "Radiated power must be <= 0.0"
    cp1d = dd.core_profiles.profiles_1d[]

    function cost(impurity_ne_fraction)
        new_impurity_fraction!(cp1d, impurity_name, impurity_ne_fraction)
        radiation_source!(dd)
        c = ((total_radiation_sources(dd).electrons.power_inside[end] - total_radiated_power) / total_radiated_power) .^ 2
        return c
    end

    impurity_ne_fraction = Optim.optimize(cost, 0.0, 1 / ion_properties(impurity_name).z_n, Optim.Brent()).minimizer
    cost(impurity_ne_fraction)

    return dd
end

"""
    new_impurity_radiation!(dd::IMAS.dd, impurity_name::Symbol, time_total_radiated_power::AbstractVector{Float64}, value_total_radiated_power::AbstractVector{<:Real})
"""
function new_impurity_radiation!(dd::IMAS.dd, impurity_name::Symbol, time_total_radiated_power::AbstractVector{Float64}, value_total_radiated_power::AbstractVector{<:Real})
    total_radiated_power_itp = IMAS.interp1d(time_total_radiated_power, value_total_radiated_power)

    global_time_bkp = dd.global_time
    for time0 in dd.core_profiles.time
        dd.global_time = time0
        try
            IMAS.new_impurity_radiation!(dd, impurity_name, total_radiated_power_itp(time0))
        catch e
            @warn("$(typeof(e)): Could not set new_impurity_radiation!(dd, $(impurity_name), $(total_radiated_power_itp(time0))) at time $(dd.global_time)")
        end
    end

    dd.global_time = global_time_bkp
    return dd
end

@compat public new_impurity_radiation!
push!(document[Symbol("Physics profiles")], :new_impurity_radiation!)


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
        z .+= ion.density_thermal .* Zi .^ 2
    end
    ne = cp1d.electrons.density_thermal
    @. z = max(one(eltype(z)), z / ne)
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
    avgZ(Z::Real,Ti::T)::T

Returns average ionization state of an ion at a given temperature
"""
function avgZ(Z::Real, Ti::T) where {T}
    func = avgZinterpolator(joinpath(@__DIR__, "..", "..", "data", "Zavg_z_t.dat"))
    return @. 10.0 ^ (func(log10(Ti / 1E3), Z)) - 1.0
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
    return core_edge_energy(p, cp1d.grid.rho_tor_norm, cp1d.grid.volume, rho_ped)
end

"""
    core_edge_energy(eqt::IMAS.equilibrium__time_slice, rho_ped::Real)

Evaluate stored energy in the core (< rho_ped) or the pedestal (> rho_ped)
"""
function core_edge_energy(eqt::IMAS.equilibrium__time_slice, rho_ped::Real)
    eqt1d = eqt.profiles_1d
    return core_edge_energy(eqt1d.pressure, eqt1d.rho_tor_norm, eqt1d.volume, rho_ped)
end

function core_edge_energy(pressure::AbstractVector{T}, rho::AbstractVector{T}, volume::AbstractVector{T}, rho_ped::Real) where {T<:Real}
    rho_bound_idx = argmin_abs(rho, rho_ped)
    pedge = pressure[rho_bound_idx]
    fedge = (k, x) -> (k <= rho_bound_idx) ? pedge : pressure[k]
    fcore = (k, x) -> (k <= rho_bound_idx) ? (pressure[k] - pedge) : 0.0
    core_value = 1.5 * trapz(volume, fcore)
    edge_value = 1.5 * trapz(volume, fedge)
    return (core_value=core_value, edge_value=edge_value)
end

@compat public core_edge_energy
push!(document[Symbol("Physics profiles")], :core_edge_energy)

"""
    scale_ion_densities_to_target_zeff(cp1d::IMAS.core_profiles__profiles_1d{T}, rho_scale::Real, target_zeff::Real) where {T<:Real}

Returns scale coefficients for main ions and impurity density profiles needed to achieve a target zeff at a specific radial location
"""
function scale_ion_densities_to_target_zeff(cp1d::IMAS.core_profiles__profiles_1d{T}, rho_scale::Real, target_zeff::Real) where {T<:Real}
    rho_index = argmin_abs(cp1d.grid.rho_tor_norm, rho_scale)
    ne = cp1d.electrons.density_thermal[rho_index]
    Ti = cp1d.t_i_average[rho_index]
    nh = zero(T)
    nim_Z = zero(T)
    nim_Z2 = zero(T)
    for ion in cp1d.ion
        if !ismissing(ion, :density_thermal)
            if ion.element[1].z_n == 1.0
                nh += ion.density_thermal[rho_index]
            else
                nimp = ion.density_thermal[rho_index]
                Zimp = IMAS.avgZ(ion.element[1].z_n, Ti)
                nim_Z += nimp .* Zimp
                nim_Z2 += nimp .* Zimp^2
            end
        end
    end

    impurity_scale = (target_zeff .- 1.0) .* ne ./ (nim_Z2 .- nim_Z)
    manion_scale = (ne .- impurity_scale .* nim_Z) ./ nh
    @assert all(impurity_scale .>= 0.0) "all(impurity_scale .>= 0.0) ", string(impurity_scale)
    @assert all(manion_scale .>= 0.0) "all(manion_scale .>= 0.0) ", string(manion_scale)
    original_zeff = (nh .+ nim_Z2) ./ (nh .+ nim_Z)
    new_zeff = (manion_scale .* nh .+ impurity_scale .* nim_Z2) ./ (manion_scale .* nh .+ impurity_scale .* nim_Z)
    original_quasineutrality = (ne .- nh .- nim_Z) ./ ne
    new_quasineutrality = (ne .- manion_scale .* nh .- impurity_scale .* nim_Z) ./ ne

    return (
        impurity_scale=impurity_scale,
        manion_scale=manion_scale,
        original_zeff=original_zeff,
        new_zeff=new_zeff,
        original_quasineutrality=original_quasineutrality,
        new_quasineutrality=new_quasineutrality
    )
end

"""
    scale_ion_densities_to_target_zeff(cp1d::IMAS.core_profiles__profiles_1d{T}, target_zeff::Vector{T}) where {T<:Real}

Returns scale coefficients for main ions and impurity density profiles needed to achieve a target zeff profile
"""
function scale_ion_densities_to_target_zeff(cp1d::IMAS.core_profiles__profiles_1d{T}, target_zeff::Vector{T}) where {T<:Real}
    ne = cp1d.electrons.density_thermal
    Ti = cp1d.t_i_average
    nh = zero(ne)
    nim_Z = zero(ne)
    nim_Z2 = zero(ne)
    for ion in cp1d.ion
        if !ismissing(ion, :density_thermal)
            if ion.element[1].z_n == 1.0
                nh .+= ion.density_thermal
            else
                nimp = ion.density_thermal
                Zimp = IMAS.avgZ(ion.element[1].z_n, Ti)
                nim_Z .+= nimp .* Zimp
                nim_Z2 .+= nimp .* Zimp .^ 2
            end
        end
    end

    impurity_scale = (target_zeff .- 1.0) .* ne ./ (nim_Z2 .- nim_Z)
    manion_scale = (ne .- impurity_scale .* nim_Z) ./ nh
    @assert all(impurity_scale .>= 0.0) "all(impurity_scale .>= 0.0) ", string(impurity_scale)
    @assert all(manion_scale .>= 0.0) "all(manion_scale .>= 0.0) ", string(manion_scale)
    original_zeff = (nh .+ nim_Z2) ./ (nh .+ nim_Z)
    new_zeff = (manion_scale .* nh .+ impurity_scale .* nim_Z2) ./ (manion_scale .* nh .+ impurity_scale .* nim_Z)
    original_quasineutrality = (ne .- nh .- nim_Z) ./ ne
    new_quasineutrality = (ne .- manion_scale .* nh .- impurity_scale .* nim_Z) ./ ne

    return (
        impurity_scale=impurity_scale,
        manion_scale=manion_scale,
        original_zeff=original_zeff,
        new_zeff=new_zeff,
        original_quasineutrality=original_quasineutrality,
        new_quasineutrality=new_quasineutrality
    )
end

@compat public scale_ion_densities_to_target_zeff
push!(document[Symbol("Physics profiles")], :scale_ion_densities_to_target_zeff)

"""
    scale_ion_densities_to_target_zeff!(cp1d::IMAS.core_profiles__profiles_1d{T}, rho_scale::Real, target_zeff::Real) where {T<:Real}

Scale main ions and impurity density profiles in place to achieve a target zeff at a specific radial location
"""
function scale_ion_densities_to_target_zeff!(cp1d::IMAS.core_profiles__profiles_1d{T}, rho_scale::Real, target_zeff::Real) where {T<:Real}
    @assert 0.0 <= rho_scale <= 1.0
    scales = scale_ion_densities_to_target_zeff(cp1d, rho_scale, target_zeff)
    for ion in cp1d.ion
        if !ismissing(ion, :density_thermal)
            if ion.element[1].z_n == 1.0
                ion.density_thermal .= ion.density_thermal .* scales.manion_scale
            else
                ion.density_thermal .= ion.density_thermal .* scales.impurity_scale
            end
        end
    end
    return cp1d
end

"""
    scale_ion_densities_to_target_zeff!(cp1d::IMAS.core_profiles__profiles_1d{T}, target_zeff::Vector{T}) where {T<:Real}

Scale main ions and impurity density profiles in place to achieve a target zeff profile
"""
function scale_ion_densities_to_target_zeff!(cp1d::IMAS.core_profiles__profiles_1d{T}, target_zeff::Vector{T}) where {T<:Real}
    scales = scale_ion_densities_to_target_zeff(cp1d, target_zeff)
    for ion in cp1d.ion
        if !ismissing(ion, :density_thermal)
            if ion.element[1].z_n == 1.0
                ion.density_thermal .= ion.density_thermal .* scales.manion_scale
            else
                ion.density_thermal .= ion.density_thermal .* scales.impurity_scale
            end
        end
    end
    return cp1d
end

@compat public scale_ion_densities_to_target_zeff!
push!(document[Symbol("Physics profiles")], :scale_ion_densities_to_target_zeff!)

"""
    flatten_profile!(prof::AbstractVector{T}, rho::AbstractVector{T}, volume_or_area::AbstractVector{T}, rho0::T, width::T) where {T<:Real}

Flatten input profile inside of rho0 with given width, while keeping integral (over volume or area) unchanged
"""
function flatten_profile!(prof::AbstractVector{T}, rho::AbstractVector{T}, volume_or_area::AbstractVector{T}, rho0::T, width::T) where {T<:Real}
    i_inversion = argmin_abs(rho, rho0)
    if i_inversion > 2
        prof_tanh = 0.5 .* tanh.(-(rho .- rho[i_inversion]) ./ width) .+ 0.5

        en = trapz(volume_or_area, prof)
        en_tanh = trapz(volume_or_area, prof_tanh)
        prof[1:i_inversion] .= prof[i_inversion]

        en_cut = trapz(volume_or_area, prof)
        prof .+= (en .- en_cut) ./ en_tanh .* prof_tanh
    end

    return prof
end

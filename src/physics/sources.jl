

"""
    function slowing_down_time(ne::AbstractVector{<:Real}, te::AbstractVector{<:Real}, mfast::Real, zfast::Real)

Calculates the slowing down time taus [Stix, Plasma Phys. 14 (1972) 367] Eq. 16

:param ne: electron density [m^-3]

:param te: electron temperature [eV]

:param efast: energy of birth ion [eV]

:param mfast: mass of fast ion [AMU]

:param zfast: fast  ion charge

:return: taus: slowing down time
""" 
function slowing_down_time(ne::AbstractVector{<:Real}, te::AbstractVector{<:Real}, mfast::Real, zfast::Real)

    ne_cm3 = 1e-6 .* ne
    
    loglam = 24.0 .- log.(sqrt.(ne_cm3) ./ te)

    taus = 6.27e8 .* mfast .* (te.^1.5) ./ (ne_cm3 .* loglam .* zfast^2)
    
    return taus
end


"""
    function critical_energy(ni::AbstractMatrix{<:Real}, zi::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, ne::AbstractVector{<:Real}, te::AbstractVector{<:Real}, mfast::Real)

Alpha heating coefficients [Stix, Plasma Phys. 14 (1972) 367]
See in particular Eqs. 15 and 17.
This function is a direct translation to Python of what in is in TGYRO

:param ni: list with thermal ions densities [m^-3]

:param zi: list with thermal ions charges

:param mi: list with thermal ions masses [AMU]

:param ne: electron density [m^-3]

:param te: electron temperature [eV]

:param efast: energy of birth ion [eV]

:param mfast: mass of fast ion [AMU]

:return: ecrit, the critical energy
"""

function critical_energy(ni::AbstractMatrix{<:Real}, zi::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, ne::AbstractVector{<:Real}, te::AbstractVector{<:Real}, mfast::Real)
    me = constants.m_e

    c_a =  zeros(size(ne))
    for k in 1:size(ni, 1)
        c_a += (ni[k,:] ./ ne) .* (zi[k] .^ 2) ./ (mi[k] .* constants.m_p ./ mfast)
    end

    ecrit = te .* (4.0 .* sqrt.(me ./ mfast) ./ (3.0 .* sqrt.(pi) .* c_a)) .^ (-2.0 / 3.0)

    return ecrit
end


"""
Calculates the fast ion density, and adds it to the dd

:param particle_energy: particle energy [eV]

:param particle_specie: particle specie 
"""
function fast_density(dd::IMAS.dd; particle_energy::Real=3.5e6, sourceid::Symbol=:fusion)

    cp1d = dd.core_profiles.profiles_1d[]

    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature

    if sourceid == :fusion
        particle_specie = "He"
    elseif sourceid == :nbi
        particle_species = ["D", "DT"]
        ion_index = findfirst(ion.label in particle_species for ion in cp1d.ion)
        particle_specie = cp1d.ion[ion_index].label
    end
    ion_index = findfirst(ion.label == particle_specie for ion in cp1d.ion)

    particle_mass = cp1d.ion[ion_index].element[1].a
    particle_charge = cp1d.ion[ion_index].element[1].z_n

    ni = zeros((length(cp1d.ion),length(cp1d.electrons.density)))
    Zi = zeros(length(cp1d.ion))
    mi = zeros(length(cp1d.ion))

    for (idx,ion) in enumerate(cp1d.ion)
        ni[idx,:] = ion.density_thermal
        Zi[idx] = ion.element[1].z_n
        mi[idx] = ion.element[1].a
    end

    taus = slowing_down_time(ne, Te, particle_charge, particle_mass)
    Ecrit = critical_energy(ni, Zi, mi, ne, Te, particle_mass*constants.m_p)

    vfrac = sqrt.(Ecrit ./ particle_energy)

    encapf = log.(1.0 .+ 1.0 ./ vfrac.^3) ./ 3.0  # assume no neutrals

    cs1ds = findall(sourceid, dd.core_sources.source)
    cp1d.ion[ion_index].pressure_fast_parallel  = zeros(length(cp1d.electrons.density))
    cp1d.ion[ion_index].pressure_fast_perpendicular  = zeros(length(cp1d.electrons.density))
    cp1d.ion[ion_index].density_fast  = zeros(length(cp1d.electrons.density))

    for cs1d in cs1ds
        qfaste = cs1d.profiles_1d[].electrons.energy
        qfasti = cs1d.profiles_1d[].total_ion_energy
        pressa = taus .* 2.0 ./ 3.0 .* qfaste

        nfast = (qfaste .+ qfasti) ./ (constants.e .* particle_energy) .* encapf .* taus
        cp1d.ion[ion_index].pressure_fast_parallel += pressa ./ 3.0
        cp1d.ion[ion_index].pressure_fast_perpendicular += pressa ./ 3.0
        cp1d.ion[ion_index].density_fast += nfast
    end
end

"""
    sivukhin_fraction(cp1d::IMAS.core_profiles__profiles_1d, particle_energy::Real, particle_mass::Real)

Compute a low-accuracy but fast approximation to the ion heating fraction (for alpha particles and beam particles).
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

"""
    reactivity(Ti::AbstractVector{<:Real}, model::String="D-T"; polarized_fuel_fraction::Real=0.0)

Fusion reactivity coming from H.-S. Bosch and G.M. Hale, Nucl. Fusion 32 (1992) 611.
"""

function reactivity(Ti::AbstractVector{<:Real}, model::String="D-T"; polarized_fuel_fraction::Real=0.0)
    spf = 1 #default value for non spin polarized fuel 
    if model == "D-T"
        # Table VII
        c1 = 1.17302e-9
        c2 = 1.51361e-2
        c3 = 7.51886e-2
        c4 = 4.60643e-3
        c5 = 1.3500e-2
        c6 = -1.06750e-4
        c7 = 1.36600e-5
        bg = 34.3827
        er = 1.124656e6
        if polarized_fuel_fraction > 0.0
            spf = 1.5 #spin polarization factor - 1.5 chosen according to GACP 20010393
        end
    elseif model == "D-He3"
        bg = 68.7508
        mc2 = 1124572.0
        c1 = 5.51036e-10
        c2 = 6.41918e-3
        c3 = -2.02896e-3
        c4 = -1.91080e-5
        c5 = 1.35776e-4
        c6 = 0.0
        c7 = 0.0
        er = 18.3e6
        if polarized_fuel_fraction > 0.0
            spf = 1.5 #also 1.5 for D-He3 according to Kulsrud (1982), PRL 49(17), 1248-1251
        end
    elseif model == "D-DtoT"
        bg = 31.3970
        mc2 = 937814.0
        c1 = 5.65718e-12
        c2 = 3.41267e-3
        c3 = 1.99167e-3
        c4 = 0.0
        c5 = 1.05060e-5
        c6 = 0.0
        c7 = 0.0
        er = 4.03e6
        if polarized_fuel_fraction > 0.0
            error("Sorry, spin polarized fuel option is not available for $(model)")
        end
    elseif model == "D-DtoHe3"
        bg = 31.3970
        mc2 = 937814.0
        c1 = 5.43360e-12
        c2 = 5.85778e-3
        c3 = 7.68222e-3
        c4 = 0.0
        c5 = -2.96400e-6
        c6 = 0.0
        c7 = 0.0
        er = 0.82e6
        if polarized_fuel_fraction > 0.0
            error("Sorry, spin polarized fuel option is not available for $(model)")
        end
    else
        error("Reactivity model can be either [\"D-T\",\"D-He3\",\"D-DtoT\", \"D-DtoHe3\"]")
    end

    Ti = Ti ./ 1e3  # from eV to keV

    r0 = Ti .* (c2 .+ Ti .* (c4 .+ Ti .* c6)) ./ (1.0 .+ Ti .* (c3 .+ Ti .* (c5 .+ Ti .* c7)))
    theta = Ti ./ (1.0 .- r0)
    xi = (bg .^ 2 ./ (4.0 .* theta)) .^ (1.0 ./ 3.0)
    sigv = c1 .* theta .* sqrt.(xi ./ (er .* Ti .^ 3)) .* exp.(-3.0 .* xi)

    @assert 0.0 <= polarized_fuel_fraction <= 1.0 "Polarized fuel fraction should be between 0.0 and 1.0"
    return ((1.0 .- polarized_fuel_fraction) .+ spf * polarized_fuel_fraction) .* sigv / 1e6  # m^3/s
end

"""
    alpha_power(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)

Volumetric heating source of α particles coming from DT reaction [W m⁻³]
"""
function alpha_heating(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)

    fast_helium_energy = 3.5e6 * 1.6022e-19  # Joules
    # Find the right D-T density
    ion_list = [ion.label for ion in cp1d.ion]
    result = zero(cp1d.electrons.density)
    if "D" in ion_list && "T" in ion_list
        D_index = findfirst(ion -> isequal(ion, "D"), ion_list)
        n_deuterium = cp1d.ion[D_index].density
        T_index = findfirst(ion -> isequal(ion, "T"), ion_list)
        n_tritium = cp1d.ion[T_index].density
        Ti = (cp1d.ion[D_index].temperature .+ cp1d.ion[T_index].temperature) ./ 2.0
        sigv = reactivity(Ti, "D-T"; polarized_fuel_fraction)
        result .= n_deuterium .* n_tritium .* sigv .* fast_helium_energy  # J/m^3/s = W/m^3

    elseif "DT" in ion_list
        DT_index = findfirst(ion -> isequal(ion, "DT"), ion_list)
        n_deuterium = n_tritium = cp1d.ion[DT_index].density ./ 2
        Ti = cp1d.ion[DT_index].temperature
        sigv = reactivity(Ti, "D-T"; polarized_fuel_fraction)
        result .= n_deuterium .* n_tritium .* sigv .* fast_helium_energy  # J/m^3/s = W/m^3
    end

    return result
end

"""
    alpha_power(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)

Total power in α particles [W]
"""
function alpha_power(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)
    return integrate(cp1d.grid.volume, alpha_heating(cp1d; polarized_fuel_fraction))
end

function fusion_power(dd::IMAS.dd)
    return fusion_power(dd.core_profiles)
end

"""
    fusion_power(cp::IMAS.core_profiles)

Calculates the fusion power in [W]
"""
function fusion_power(cp::IMAS.core_profiles)
    cp1d = cp.profiles_1d[]
    polarized_fuel_fraction = getproperty(cp.global_quantities, :polarized_fuel_fraction, 0.0)
    return alpha_power(cp1d; polarized_fuel_fraction) * 5.0
end

function DT_fusion_source!(dd::IMAS.dd)
    return DT_fusion_source!(dd.core_sources, dd.core_profiles)
end

"""
    DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)

Calculates DT fusion heating with an estimation of the alpha slowing down to the ions and electrons, modifies dd.core_sources
"""
function DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)
    cp1d = cp.profiles_1d[]

    polarized_fuel_fraction = getproperty(cp.global_quantities, :polarized_fuel_fraction, 0.0)
    α = alpha_heating(cp1d; polarized_fuel_fraction)
    if sum(α) == 0
        deleteat!(cs.source, "identifier.index" => 6)
        return cs
    end
    ion_to_electron_fraction = sivukhin_fraction(cp1d, 3.5e6, 4.0)

    index = name_2_index(cs.source)[:fusion]
    source = resize!(cs.source, "identifier.index" => index; allow_multiple_matches=true)
    new_source(
        source,
        index,
        "α",
        cp1d.grid.rho_tor_norm,
        cp1d.grid.volume;
        electrons_energy=α .* (1.0 .- ion_to_electron_fraction),
        total_ion_energy=α .* ion_to_electron_fraction
    )

    fast_density(dd)
    
    return source
end

"""
    D_D_to_He3_reactions(dd::IMAS.dd)

Calculates the number of D-D thermal fusion reactions to He3 in [reactions/m³/s]
"""
function D_D_to_He3_reactions(cp1d::IMAS.core_profiles__profiles_1d)
    index = findfirst(ion.label == "D" for ion in cp1d.ion)
    @assert index !== nothing "There is no Deuterium only species in dd.core_profiles"
    dd_deut = cp1d.ion[index]
    return dd_deut.density_thermal .^ 2 .* reactivity(dd_deut.temperature, "D-DtoHe3"; polarized_fuel_fraction=0.0) #  reactions/m³/s
end

"""
    D_D_to_He3_source!(dd::IMAS.dd)

Calculates the He-3 heating source from D-D fusion reactions, estimates energy transfer to ions and electrons, modifies dd.core_sources
"""
function D_D_to_He3_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    he3_energy = 0.82e6 * constants.e  # eV to Joules
    reactivity = D_D_to_He3_reactions(cp1d)
    energy =  reactivity .* he3_energy 
    ion_to_electron_fraction = sivukhin_fraction(cp1d, 0.82e6, 3.0)
    index = name_2_index(dd.core_sources.source)[:fusion]
    source = resize!(dd.core_sources.source, "identifier.index" => index; allow_multiple_matches=true)
    new_source(
        source,
        index,
        "He3",
        cp1d.grid.rho_tor_norm,
        cp1d.grid.volume;
        electrons_energy=energy .* (1.0 .- ion_to_electron_fraction),
        total_ion_energy=energy .* ion_to_electron_fraction
    )
    return source
end

"""
    collisional_exchange_source!(dd::IMAS.dd)

Calculates collisional exchange source and modifies dd.core_sources
"""
function collisional_exchange_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature
    Ti = cp1d.ion[1].temperature

    index = name_2_index(dd.core_sources.source)[:collisional_equipartition]
    if all(Te .≈ Ti)
        deleteat!(dd.core_sources.source, "identifier.index" => index)
    else
        nu_exch = collision_frequencies(dd)[3]
        delta = 1.5 .* nu_exch .* ne .* constants.e .* (Te .- Ti)
        source = resize!(dd.core_sources.source, "identifier.index" => index; allow_multiple_matches=true)
        new_source(source, index, "exchange", cp1d.grid.rho_tor_norm, cp1d.grid.volume; electrons_energy=-delta, total_ion_energy=delta)
        return source
    end
end

"""
    ohmic_source!(dd::IMAS.dd)

Calculates the ohmic source and modifies dd.core_sources
"""
function ohmic_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    j_ohmic = getproperty(cp1d, :j_ohmic, missing)
    if j_ohmic !== missing
        powerDensityOhm = j_ohmic .^ 2 ./ cp1d.conductivity_parallel
        index = name_2_index(dd.core_sources.source)[:ohmic]
        source = resize!(dd.core_sources.source, "identifier.index" => index)
        new_source(source, index, "ohmic", cp1d.grid.rho_tor_norm, cp1d.grid.volume; electrons_energy=powerDensityOhm, j_parallel=j_ohmic)
        return source
    end
end

"""
    bootstrap_source!(dd::IMAS.dd)

Calculates the bootsrap current source and modifies dd.core_sources
"""
function bootstrap_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    j_bootstrap = getproperty(cp1d, :j_bootstrap, missing)
    if j_bootstrap !== missing
        index = name_2_index(dd.core_sources.source)[:bootstrap_current]
        source = resize!(dd.core_sources.source, "identifier.index" => index)
        new_source(source, index, "bootstrap", cp1d.grid.rho_tor_norm, cp1d.grid.volume; j_parallel=j_bootstrap)
        return source
    end
end

"""
    sources!(dd::IMAS.dd)

Calculates intrisic sources and sinks and adds them to dd.core_sources
"""
function sources!(dd::IMAS.dd)
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)
    IMAS.collisional_exchange_source!(dd)
    IMAS.bremsstrahlung_source!(dd)
    IMAS.line_radiation_source!(dd)
    IMAS.synchrotron_source!(dd)
    IMAS.DT_fusion_source!(dd)
end

"""
    total_power_source(source::IMAS.core_sources__source___profiles_1d)

Returns the total power (electron + ion) for a single source
"""
function total_power_source(source::IMAS.core_sources__source___profiles_1d)
    return getproperty(source.electrons, :power_inside, [0.0])[end] + getproperty(source, :total_ion_power_inside, [0.0])[end]
end

"""
    total_power_time(core_sources::IMAS.core_sources, include_indexes::Vector{<:Integer})

Returns tuple of vectors with the total thermal power and time_array for given set of sources selected by identifier.index
"""
function total_power_time(core_sources::IMAS.core_sources, include_indexes::Vector{<:Integer})
    sources = []
    for index in include_indexes
        append!(sources, findall(core_sources.source, "identifier.index" => index))
    end
    time_array = core_sources.time
    total_power = zeros(length(time_array))
    for source in sources
        total_power .+= [total_power_source(source.profiles_1d[t]) for t in time_array]
    end
    return total_power, time_array
end

function total_sources(dd::IMAS.dd)
    total_sources(dd.core_sources, dd.core_profiles.profiles_1d[])
end

"""
    total_sources(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; include_indexes=missing, exclude_indexes=missing)

Returns core_sources__source___profiles_1d with sources totals and possiblity to explicitly include/exclude certain sources based on their unique index identifier.
"""
function total_sources(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; include_indexes=missing, exclude_indexes=missing)
    total_source1d = IMAS.core_sources__source___profiles_1d()
    total_source1d.grid.rho_tor_norm = rho = cp1d.grid.rho_tor_norm
    total_source1d.time = cp1d.time

    for prop in [:volume, :area, :surface]
        value = getproperty(cp1d.grid, prop, missing)
        if value === missing
            for source in core_sources.source
                value = getproperty(source.profiles_1d[Float64(cp1d.time)].grid, prop, missing)
                if value !== missing
                    break
                end
            end
        end
        if value !== missing
            setproperty!(total_source1d.grid, prop, value)
        end
    end

    all_indexes = [source.identifier.index for source in core_sources.source]

    for source in core_sources.source
        if include_indexes !== missing && source.identifier.index ∈ include_indexes
            # pass
        elseif include_indexes !== missing
            continue
        end

        if source.identifier.index in [0]
            @debug "total_sources() skipping unspecified source with index $(source.identifier.index)"
            continue
        elseif 107 >= source.identifier.index >= 100
            @debug "total_sources() skipping combination source with index $(source.identifier.index)"
            continue
        elseif (source.identifier.index) in [1] && any(all_indexes .> 1)
            @debug "total_sources() skipping total source with index $(source.identifier.index)"
            continue
        elseif (source.identifier.index) in [200] && any(300 > all_indexes > 200)
            @debug "total_sources() skipping total radiation source with index $(source.identifier.index)"
            continue
        elseif exclude_indexes !== missing && source.identifier.index ∈ exclude_indexes
            continue
        end
        source_name = ismissing(source.identifier, :name) ? "?" : source.identifier.name

        if isempty(source.profiles_1d)
            continue
        end
        @debug "total_sources() including $source_name source with index $(source.identifier.index)"
        source1d = source.profiles_1d[Float64(cp1d.time)]
        for sub in [nothing, :electrons]
            ids1 = total_source1d
            ids2 = source1d
            if sub !== nothing
                ids1 = getproperty(ids1, sub)
                ids2 = getproperty(ids2, sub)
            end
            for field in keys(ids1)
                y = getproperty(ids2, field, missing)
                if typeof(y) <: AbstractVector{<:Real}
                    if typeof(getraw(ids1, field)) <: Union{Missing,Function}
                        setproperty!(ids1, field, zeros(length(total_source1d.grid.rho_tor_norm)))
                    end
                    old_value = getproperty(ids1, field)
                    x = source1d.grid.rho_tor_norm
                    setproperty!(ids1, field, old_value .+ interp1d(x, y).(rho))
                end
            end
        end
    end

    # assign zeros to missing fields of total_sources
    for sub in [nothing, :electrons]
        ids1 = total_source1d
        if sub !== nothing
            ids1 = getproperty(ids1, sub)
        end
        for field in keys(ids1)
            if ismissing(ids1, field)
                setproperty!(ids1, field, zeros(size(rho)))
            end
        end
    end

    return total_source1d
end

"""
    new_source(
        source::IMAS.core_sources__source,
        index::Int,
        name::String,
        rho::Union{AbstractVector,AbstractRange},
        volume::Union{AbstractVector,AbstractRange};
        electrons_energy::Union{AbstractVector,Missing}=missing,
        electrons_power_inside::Union{AbstractVector,Missing}=missing,
        total_ion_energy::Union{AbstractVector,Missing}=missing,
        total_ion_power_inside::Union{AbstractVector,Missing}=missing,
        electrons_particles::Union{AbstractVector,Missing}=missing,
        electrons_particles_inside::Union{AbstractVector,Missing}=missing,
        j_parallel::Union{AbstractVector,Missing}=missing,
        current_parallel_inside::Union{AbstractVector,Missing}=missing,
        momentum_tor::Union{AbstractVector,Missing}=missing,
        torque_tor_inside::Union{AbstractVector,Missing}=missing
    )

Populates the IMAS.core_sources__source with given heating, particle, current, momentun profiles
"""
function new_source(
    source::IMAS.core_sources__source,
    index::Int,
    name::String,
    rho::Union{AbstractVector,AbstractRange},
    volume::Union{AbstractVector,AbstractRange};
    electrons_energy::Union{AbstractVector,Missing}=missing,
    electrons_power_inside::Union{AbstractVector,Missing}=missing,
    total_ion_energy::Union{AbstractVector,Missing}=missing,
    total_ion_power_inside::Union{AbstractVector,Missing}=missing,
    electrons_particles::Union{AbstractVector,Missing}=missing,
    electrons_particles_inside::Union{AbstractVector,Missing}=missing,
    j_parallel::Union{AbstractVector,Missing}=missing,
    current_parallel_inside::Union{AbstractVector,Missing}=missing,
    momentum_tor::Union{AbstractVector,Missing}=missing,
    torque_tor_inside::Union{AbstractVector,Missing}=missing
)

    source.identifier.name = name
    source.identifier.index = index
    cs1d = resize!(source.profiles_1d)
    cs1d.grid.rho_tor_norm = rho
    cs1d.grid.volume = volume

    if electrons_energy !== missing
        cs1d.electrons.energy = interp1d(LinRange(0, 1, length(electrons_energy)), electrons_energy).(cs1d.grid.rho_tor_norm)
    end
    if electrons_power_inside !== missing
        cs1d.electrons.power_inside = interp1d(LinRange(0, 1, length(electrons_power_inside)), electrons_power_inside).(cs1d.grid.rho_tor_norm)
    end

    if total_ion_energy !== missing
        cs1d.total_ion_energy = interp1d(LinRange(0, 1, length(total_ion_energy)), total_ion_energy).(cs1d.grid.rho_tor_norm)
    end
    if total_ion_power_inside !== missing
        cs1d.total_ion_power_inside = interp1d(LinRange(0, 1, length(total_ion_power_inside)), total_ion_power_inside).(cs1d.grid.rho_tor_norm)
    end

    if electrons_particles !== missing
        cs1d.electrons.particles = interp1d(LinRange(0, 1, length(electrons_particles)), electrons_particles).(cs1d.grid.rho_tor_norm)
    end
    if electrons_particles_inside !== missing
        cs1d.electrons.particles_inside = interp1d(LinRange(0, 1, length(electrons_particles_inside)), electrons_particles_inside).(cs1d.grid.rho_tor_norm)
    end

    if j_parallel !== missing
        cs1d.j_parallel = interp1d(LinRange(0, 1, length(j_parallel)), j_parallel).(cs1d.grid.rho_tor_norm)
    end
    if current_parallel_inside !== missing
        cs1d.current_parallel_inside = interp1d(LinRange(0, 1, length(current_parallel_inside)), current_parallel_inside).(cs1d.grid.rho_tor_norm)
    end

    if momentum_tor !== missing
        cs1d.momentum_tor = interp1d(LinRange(0, 1, length(momentum_tor)), momentum_tor).(cs1d.grid.rho_tor_norm)
    end
    if torque_tor_inside !== missing
        cs1d.torque_tor_inside = interp1d(LinRange(0, 1, length(torque_tor_inside)), torque_tor_inside).(cs1d.grid.rho_tor_norm)
    end

    return source
end

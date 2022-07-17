"""
    sivukhin_fraction(cp1d::IMAS.core_profiles__profiles_1d, particle_energy::Real, particle_mass::Real)

Compute a low-accuracy but fast approximation to the ion heating fraction (for alpha particles and beam particles).
"""
function sivukhin_fraction(cp1d::IMAS.core_profiles__profiles_1d, particle_energy::Real, particle_mass::Real)
    Te = cp1d.electrons.temperature
    ne = cp1d.electrons.density
    rho = cp1d.grid.rho_tor_norm

    particle_mass = particle_mass * constants.m_p

    tp = typeof(promote(Te[1], ne[1], rho[1])[1])
    c_a = zeros(tp, length(rho))
    W_crit = similar(c_a)
    ion_elec_fraction = similar(W_crit)
    for ion in cp1d.ion
        ni = ion.density
        Zi = ion.element[1].z_n
        mi = ion.element[1].a * constants.m_p
        c_a .+= (ni ./ ne) .* Zi .^ 2 .* (mi ./ particle_mass)
    end

    W_crit = Te .* (4.0 .* sqrt.(constants.m_e / particle_mass) ./ (3.0 * sqrt(pi) .* c_a)) .^ (-2.0 / 3.0)

    x = particle_energy ./ W_crit
    for (idx, x_i) in enumerate(x)
        y = x_i .* rho
        f = integrate(y, 1.0 ./ (1.0 .+ y .^ 1.5))
        ion_elec_fraction[idx] = f / x_i
    end

    return ion_elec_fraction
end

"""
    alpha_power(cp1d::IMAS.core_profiles__profiles_1d)

Volumetric heating source of α particles coming from DT reaction [W m⁻³]
Based on Table VII of H.-S. Bosch and G.M. Hale, Nucl. Fusion 32 (1992) 611.
"""
function alpha_heating(cp1d::IMAS.core_profiles__profiles_1d)
    fast_helium_energy = 3.5e6 * 1.6022e-12 * 1e-7  # Joules

    c1 = 1.17302e-9
    c2 = 1.51361e-2
    c3 = 7.51886e-2
    c4 = 4.60643e-3
    c5 = 1.3500e-2
    c6 = -1.06750e-4
    c7 = 1.36600e-5
    bg = 34.3827
    er = 1.124656e6

    # Find the right D-T density
    ion_list = [ion.label for ion in cp1d.ion]
    if "D" in ion_list && "T" in ion_list && length(findall(ion -> isequal(ion, "T"), ion_list)) < 2
        D_index = findfirst(ion -> isequal(ion, "D"), ion_list)
        n_deuterium = cp1d.ion[D_index].density
        T_index = findfirst(ion -> isequal(ion, "T"), ion_list)
        n_tritium = cp1d.ion[T_index].density
        Ti = (cp1d.ion[D_index].temperature + cp1d.ion[T_index].temperature) ./ 2.0 .* 1e-3 # keV
    elseif "DT" in ion_list
        DT_index = findfirst(ion -> isequal(ion, "DT"), ion_list)
        n_deuterium = n_tritium = cp1d.ion[DT_index].density ./ 2
        Ti = cp1d.ion[DT_index].temperature .* 1e-3 # keV
    else
        return cp1d.electrons.density .* 0.0
    end

    r0 = Ti .* (c2 .+ Ti .* (c4 .+ Ti .* c6)) ./ (1.0 .+ Ti .* (c3 .+ Ti .* (c5 .+ Ti .* c7)))
    theta = Ti ./ (1.0 .- r0)
    xi = (bg .^ 2 ./ (4.0 .* theta)) .^ (1.0 ./ 3.0)
    sigv = c1 .* theta .* sqrt.(xi ./ (er .* Ti .^ 3)) .* exp.(-3.0 .* xi)

    reactivity = sigv / 1e6  # m^3/s

    return n_deuterium .* n_tritium .* reactivity .* fast_helium_energy  # J/m^3/s = W/m^3
end

"""
    alpha_power(cp1d::IMAS.core_profiles__profiles_1d)

Total power in α particles [W]
"""
function alpha_power(cp1d::IMAS.core_profiles__profiles_1d)
    return integrate(cp1d.grid.volume, alpha_heating(cp1d))
end

"""
    fusion_power(cp1d::IMAS.core_profiles__profiles_1d)

Total fusion power [W]
"""
function fusion_power(cp1d::IMAS.core_profiles__profiles_1d)
    return alpha_power(cp1d) * 5
end

"""
    DT_fusion_source!(dd::IMAS.dd)

Calculates DT fusion heating with an estimation of the alpha slowing down to the ions and electrons, modifies dd.core_sources
"""
function DT_fusion_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]

    α = alpha_heating(cp1d)
    if sum(α) == 0
        deleteat!(dd.core_sources.source, "identifier.index" => 6)
        return dd
    end
    ion_electron_fraction = sivukhin_fraction(cp1d, 3.5e6, 4.0)

    source = resize!(dd.core_sources.source, "identifier.index" => 6; allow_multiple_matches=true)
    new_source(
        source,
        6,
        "α",
        cp1d.grid.rho_tor_norm,
        cp1d.grid.volume;
        electrons_energy=α .* (1 .- ion_electron_fraction),
        total_ion_energy=α .* ion_electron_fraction
    )
    @ddtime(dd.summary.fusion.power.value = source.profiles_1d[].total_ion_power_inside[end] + source.profiles_1d[].electrons.power_inside[end])

    return dd
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

    if all(Te .≈ Ti)
        deleteat!(dd.core_sources.source, "identifier.index" => 11)
    else
        nu_exch = collision_frequencies(dd)[3]
        delta = 1.5 .* nu_exch .* ne .* constants.e .* (Te .- Ti)
        source = resize!(dd.core_sources.source, "identifier.index" => 11; allow_multiple_matches=true)
        new_source(source, 11, "exchange", cp1d.grid.rho_tor_norm, cp1d.grid.volume; electrons_energy=-delta, total_ion_energy=delta)
    end

    return dd
end

"""
    bremsstrahlung_source!(dd::IMAS.dd)

Calculates Bremsstrahlung radiation source and modifies dd.core_sources
"""
function bremsstrahlung_source!(dd::IMAS.dd)
    # Plasma estimated at ellipsoid torus for volume contribution (triangularity is small correction)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature

    # Bremsstrahlung radiation
    powerDensityBrem = -1.690e-38 .* ne .^ 2 .* cp1d.zeff .* sqrt.(Te)
    source = resize!(dd.core_sources.source, "identifier.index" => 8)
    new_source(source, 8, "brem", cp1d.grid.rho_tor_norm, cp1d.grid.volume; electrons_energy=powerDensityBrem)
    return dd
end

"""
    radiation_losses(sources::IMAS.core_sources)

Evaluate total plasma radiation losses [W] due to both bremsstrahlung and line radiation
Synchlotron radation is not considered since it gets reabsorbed
"""
function radiation_losses(sources::IMAS.core_sources)
    radiation_indices = [8, 10] # [brehm, line]
    radiation_energy = 0.0
    for source in sources.source
        if source.identifier.index ∈ radiation_indices
            radiation_energy += source.profiles_1d[].electrons.power_inside[end]
        end
    end
    return radiation_energy
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
        source = resize!(dd.core_sources.source, "identifier.index" => 7)
        new_source(source, 7, "ohmic", cp1d.grid.rho_tor_norm, cp1d.grid.volume; electrons_energy=powerDensityOhm, j_parallel=j_ohmic)
    end
    return dd
end

"""
    bootstrap_source!(dd::IMAS.dd)

Calculates the bootsrap current source and modifies dd.core_sources
"""
function bootstrap_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    j_bootstrap = getproperty(cp1d, :j_bootstrap, missing)
    if j_bootstrap !== missing
        source = resize!(dd.core_sources.source, "identifier.index" => 13)
        new_source(source, 13, "bootstrap", cp1d.grid.rho_tor_norm, cp1d.grid.volume; j_parallel=j_bootstrap)
    end
    return dd
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

function total_sources(dd)
    total_sources(dd.core_sources, dd.core_profiles.profiles_1d[])
end

"""
    total_sources(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; include_indexes=missing, exclude_indexes=missing)

Returns core_sources__source___profiles_1d with sources totals and possiblity to explicitly include/exclude certain sources based on their unique index identifier.
"""
function total_sources(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; include_indexes=missing, exclude_indexes=missing)
    total_source1d = IMAS.core_sources__source___profiles_1d()
    total_source1d.grid.rho_tor_norm = rho = cp1d.grid.rho_tor_norm
    if !ismissing(cp1d.grid, :volume)
        total_source1d.grid.volume = cp1d.grid.volume
    end
    if !ismissing(cp1d.grid, :area)
        total_source1d.grid.area = cp1d.grid.area
    end
    total_source1d.time = cp1d.time

    all_indexes = [source.identifier.index for source in core_sources.source]

    for source in core_sources.source
        if include_indexes !== missing && source.identifier.index ∈ include_indexes
            # pass
        elseif include_indexes !== missing
            continue
        end

        if source.identifier.index in [0]
            @warn "total_sources() skipping unspecified source with index $(source.identifier.index)"
            continue
        elseif 107 >= source.identifier.index >= 100
            @warn "total_sources() skipping combination source with index $(source.identifier.index)"
            continue
        elseif (source.identifier.index) in [1] && any(all_indexes .> 1)
            @warn "total_sources() skipping total source with index $(source.identifier.index)"
            continue
        elseif (source.identifier.index) in [200] && any(300 > all_indexes > 200)
            @warn "total_sources() skipping total radiation source with index $(source.identifier.index)"
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
                initialized = false
                if !ismissing(ids2, field)
                    y = getproperty(ids2, field)
                    if typeof(y) <: AbstractVector{<:Real}
                        if typeof(getfield(ids1, field)) <: Union{Missing,Function}
                            setproperty!(ids1, field, zeros(length(total_source1d.grid.rho_tor_norm)))
                            initialized = true
                        end
                        old_value = getproperty(ids1, field)
                        x = source1d.grid.rho_tor_norm
                        setproperty!(ids1, field, old_value .+ interp1d(x, y).(rho))
                    end
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
    resize!(source.profiles_1d)
    cs1d = source.profiles_1d[]
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
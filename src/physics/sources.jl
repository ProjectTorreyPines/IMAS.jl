"""
    fusion_source!(dd::IMAS.dd)

Calculates fusion source from D-T and D-D reactions and modifies dd.core_sources
"""
function fusion_source!(dd::IMAS.dd)
    deleteat!(dd.core_sources.source, :fusion)
    D_T_to_He4_source!(dd.core_sources, dd.core_profiles)
    D_D_to_He3_source!(dd.core_sources, dd.core_profiles)
    D_D_to_T_source!(dd.core_sources, dd.core_profiles)
    fast_density(cs, cp)
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
        deleteat!(dd.core_sources.source, :collisional_equipartition)
    else
        nu_exch = collision_frequencies(dd)[3]
        delta = 1.5 .* nu_exch .* ne .* constants.e .* (Te .- Ti)
        source = resize!(dd.core_sources.source, :collisional_equipartition; allow_multiple_matches=true)
        new_source(source, source.identifier.index, "exchange", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_energy=-delta, total_ion_energy=delta)
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
        source = resize!(dd.core_sources.source, :ohmic)
        new_source(source, source.identifier.index, "ohmic", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_energy=powerDensityOhm, j_parallel=j_ohmic)
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
        source = resize!(dd.core_sources.source, :bootstrap_current)
        new_source(source, source.identifier.index, "bootstrap", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; j_parallel=j_bootstrap)
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
    IMAS.fusion_source!(dd)
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

    for prop in (:volume, :area, :surface)
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

        if source.identifier.index == 0
            @debug "total_sources() skipping unspecified source with index $(source.identifier.index)"
            continue
        elseif 107 >= source.identifier.index >= 100
            @debug "total_sources() skipping combination source with index $(source.identifier.index)"
            continue
        elseif (source.identifier.index) == 1 && any(all_indexes .> 1)
            @debug "total_sources() skipping total source with index $(source.identifier.index)"
            continue
        elseif (source.identifier.index) == 200 && any(300 > all_indexes > 200)
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
        for sub in (nothing, :electrons)
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
    for sub in (nothing, :electrons)
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
        volume::Union{AbstractVector,AbstractRange},
        area::Union{AbstractVector,AbstractRange};
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
    volume::Union{AbstractVector,AbstractRange},
    area::Union{AbstractVector,AbstractRange};
    electrons_energy::Union{AbstractVector,Missing}=missing,
    electrons_power_inside::Union{AbstractVector,Missing}=missing,
    total_ion_energy::Union{AbstractVector,Missing}=missing,
    total_ion_power_inside::Union{AbstractVector,Missing}=missing,
    electrons_particles::Union{AbstractVector,Missing}=missing,
    electrons_particles_inside::Union{AbstractVector,Missing}=missing,
    j_parallel::Union{AbstractVector,Missing}=missing,
    current_parallel_inside::Union{AbstractVector,Missing}=missing,
    momentum_tor::Union{AbstractVector,Missing}=missing,
    torque_tor_inside::Union{AbstractVector,Missing}=missing)

    source.identifier.name = name
    source.identifier.index = index
    cs1d = resize!(source.profiles_1d)
    cs1d.grid.rho_tor_norm = rho
    cs1d.grid.volume = volume
    cs1d.grid.area = area

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

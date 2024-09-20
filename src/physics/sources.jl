"""
    fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)

Calculates fusion source from D-T and D-D reactions and adds them to `dd.core_sources`

If D+T plasma, then D+D is neglected
"""
function fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)
    cp1d = cp.profiles_1d[]
    ion_list = (ion.label for ion in cp1d.ion)
    if "T" in ion_list || "DT" in ion_list
        D_T_to_He4_source!(cs, cp)
    else
        D_D_to_He3_source!(cs, cp)
        D_D_to_T_source!(cs, cp)
    end
    return fast_particles!(cs, cp.profiles_1d[])
end

function fusion_source!(dd::IMAS.dd)
    return fusion_source!(dd.core_sources, dd.core_profiles)
end

"""
    collisional_exchange_source!(dd::IMAS.dd)

Calculates collisional exchange source and adds it to `dd.core_sources`
"""
function collisional_exchange_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature
    Ti = cp1d.t_i_average

    nu_exch = collision_frequencies(dd).nu_exch
    delta = 1.5 .* nu_exch .* ne .* constants.e .* (Te .- Ti)

    source = resize!(dd.core_sources.source, :collisional_equipartition; wipe=false)
    new_source(source, source.identifier.index, "exchange", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area;
        electrons_energy=-delta, total_ion_energy=delta)
    return source
end

"""
    ohmic_source!(dd::IMAS.dd)

Calculates the ohmic source from data in `dd.core_profiles` and adds it to `dd.core_sources`
"""
function ohmic_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    if !ismissing(cp1d, :j_ohmic)
        eqt = dd.equilibrium.time_slice[]
        eqt1d = eqt.profiles_1d
        rho_tor_norm = cp1d.grid.rho_tor_norm
        rho_eq = eqt1d.rho_tor_norm
        gm1 = interp1d(rho_eq, eqt1d.gm1, :cubic).(rho_tor_norm)
        gm9 = interp1d(rho_eq, eqt1d.gm9, :cubic).(rho_tor_norm)
        f = interp1d(rho_eq, eqt1d.f, :cubic).(rho_tor_norm)
        powerDensityOhm = (cp1d.j_tor .* gm9) .* (cp1d.j_ohmic .* eqt.global_quantities.vacuum_toroidal_field.b0) ./ (f .* gm1 .* cp1d.conductivity_parallel)
        source = resize!(dd.core_sources.source, :ohmic; wipe=false)
        new_source(source, source.identifier.index, "ohmic", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area;
            electrons_energy=powerDensityOhm, j_parallel=cp1d.j_ohmic)
        return source
    end
end

"""
    bootstrap_source!(dd::IMAS.dd)

Calculates the bootsrap current source from data in `dd.core_profiles` and adds it to `dd.core_sources`
"""
function bootstrap_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    if !ismissing(cp1d, :j_bootstrap)
        source = resize!(dd.core_sources.source, :bootstrap_current; wipe=false)
        new_source(source, source.identifier.index, "bootstrap", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area;
            j_parallel=cp1d.j_bootstrap)
        return source
    end
end

"""
    sources!(dd::IMAS.dd; bootstrap::Bool=true, ohmic::Bool=true)

Calculates intrisic sources and sinks and adds them to `dd.core_sources`
"""
function sources!(dd::IMAS.dd; bootstrap::Bool=true, ohmic::Bool=true)
    if bootstrap
        IMAS.bootstrap_source!(dd)
    end
    if ohmic
        IMAS.ohmic_source!(dd)
    end
    IMAS.collisional_exchange_source!(dd)
    IMAS.bremsstrahlung_source!(dd)
    IMAS.line_radiation_source!(dd)
    IMAS.synchrotron_source!(dd)
    IMAS.fusion_source!(dd)
    return nothing
end

"""
    time_derivative_source!(dd::IMAS.dd, cp1d_old::IMAS.core_profiles__profiles_1d, Δt::Float64)

Calculates the time dependent sources and sinks and adds them to `dd.core_sources`
"""
function time_derivative_source!(dd::IMAS.dd, cp1d_old::IMAS.core_profiles__profiles_1d, Δt::Float64)
    cp1d = dd.core_profiles.profiles_1d[]
    eqt1d = dd.equilibrium.time_slice[].profiles_1d

    R_flux_avg = IMAS.interp1d(eqt1d.rho_tor_norm, eqt1d.gm8).(cp1d.grid.rho_tor_norm)
    ddt_sources = time_derivative_source!(cp1d, cp1d_old, Δt, R_flux_avg)

    source = resize!(dd.core_sources.source, :time_derivative; wipe=false)
    new_source(source, source.identifier.index, "∂/∂t term", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area;
        electrons_energy=ddt_sources.Qe,
        total_ion_energy=ddt_sources.Qi,
        electrons_particles=ddt_sources.Sne,
        momentum_tor=ddt_sources.PI)

    return source
end

"""
    time_derivative_source!(cp1d_new::IMAS.core_profiles__profiles_1d, cp1d_old::IMAS.core_profiles__profiles_1d, Δt::Float64, R_flux_avg::Vector)

These are the ∂/∂t term in the transport equations
"""
function time_derivative_source!(cp1d_new::IMAS.core_profiles__profiles_1d, cp1d_old::IMAS.core_profiles__profiles_1d, Δt::Float64, R_flux_avg::Vector)
    Sne = -(cp1d_new.electrons.density .- cp1d_old.electrons.density) / Δt

    Qe = -1.5 * (cp1d_new.electrons.pressure .- cp1d_old.electrons.pressure) / Δt

    Qi = -1.5 * (cp1d_new.pressure_ion_total .- cp1d_old.pressure_ion_total) / Δt

    d1 = cp1d_new.rotation_frequency_tor_sonic .* total_mass_density(cp1d_new) .* R_flux_avg / Δt
    d2 = cp1d_old.rotation_frequency_tor_sonic .* total_mass_density(cp1d_old) .* R_flux_avg / Δt
    PI = d2 - d1

    return (Sne=Sne, Qi=Qi, Qe=Qe, PI=PI)
end

"""
    total_mass_density(cp1d::IMAS.core_profiles__profiles_1d)
"""
function total_mass_density(cp1d::IMAS.core_profiles__profiles_1d)
    mass_density = constants.m_e * cp1d.electrons.density
    for ion in cp1d.ion
        mass_density .+= ion.density * ion.element[1].a * constants.m_p
    end
    return mass_density
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
    sources = IMAS.core_sources__source[]
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

function total_sources(dd::IMAS.dd; time0::Float64=dd.global_time, kw...)
    return total_sources(dd.core_sources, dd.core_profiles.profiles_1d[time0]; kw...)
end

"""
    total_sources(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; include_indexes=missing, exclude_indexes=missing)

Returns core_sources__source___profiles_1d with sources totals and possiblity to

  - include/exclude certain sources based on their unique index identifier
  - include only certain fields among these: [:particles_inside, :energy, :power_inside, :momentum_tor, :total_ion_power_inside, :total_ion_energy, :j_parallel, :torque_tor_inside, :current_parallel_inside, :particles]
"""
function total_sources(
    core_sources::IMAS.core_sources{T},
    cp1d::IMAS.core_profiles__profiles_1d{T};
    include_indexes::Vector{Int}=Int[],
    exclude_indexes::Vector{Int}=Int[],
    fields::Vector{Symbol}=Symbol[],
    only_positive_negative::Int=0
) where {T<:Real}

    total_source1d = IMAS.core_sources__source___profiles_1d{T}()
    total_source1d.grid.rho_tor_norm = rho = cp1d.grid.rho_tor_norm
    total_source1d.time = cp1d.time

    matching = Dict{Symbol,Symbol}()
    matching[:power_inside] = :energy
    matching[:energy] = :power_inside
    matching[:total_ion_power_inside] = :total_ion_energy
    matching[:total_ion_energy] = :total_ion_power_inside
    matching[:particles_inside] = :particles
    matching[:particles] = :particles_inside
    matching[:current_parallel_inside] = :j_parallel
    matching[:j_parallel] = :current_parallel_inside
    matching[:torque_tor_inside] = :momentum_tor
    matching[:momentum_tor] = :torque_tor_inside

    @assert isempty(fields) || all(field in keys(matching) for field in fields) "Supported fields are $(collect(keys(matching)))"

    for prop in (:volume, :area, :surface)
        value = getproperty(cp1d.grid, prop, missing)
        if value === missing
            for source in core_sources.source
                value = getproperty(source.profiles_1d[Float64(cp1d.time)].grid, prop, missing)
                if value !== missing
                    break
                end
            end
        else
            setproperty!(total_source1d.grid, prop, value)
        end
    end

    # initialize ions
    total_source1d_ions = IMAS.core_sources__source___profiles_1d___ion[]
    for source in core_sources.source
        if source.profiles_1d[1].time <= Float64(cp1d.time)
            source1d = source.profiles_1d[Float64(cp1d.time)]
            for ion in source1d.ion
                l = length(total_source1d.ion)
                tmp = resize!(total_source1d.ion, "element[1].a" => ion.element[1].z_n, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label)
                if l != length(total_source1d.ion)
                    push!(total_source1d_ions, tmp)
                end
            end
        end
    end

    # zero out total_sources
    for ids1 in [[total_source1d, total_source1d.electrons]; total_source1d_ions]
        for field in keys(ids1)
            if field in keys(matching)
                setproperty!(ids1, field, zeros(T, size(rho)))
            end
        end
    end

    # start accumulating
    all_indexes = [source.identifier.index for source in core_sources.source]
    for source in core_sources.source
        if isempty(include_indexes) || source.identifier.index ∈ include_indexes
            # pass
        else
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
        elseif (source.identifier.index) == 200 && any(300 .> all_indexes .> 200)
            @debug "total_sources() skipping total radiation source with index $(source.identifier.index)"
            continue
        elseif source.identifier.index ∈ exclude_indexes
            continue
        end
        if isempty(source.profiles_1d)
            continue # skip sources that have no profiles_1d
        end
        if source.profiles_1d[1].time > Float64(cp1d.time)
            continue # skip sources that start after time of interest
        end

        source_name = ismissing(source.identifier, :name) ? "?" : source.identifier.name
        @debug "total_sources() including $source_name source with index $(source.identifier.index)"

        source1d = source.profiles_1d[Float64(cp1d.time)]

        # ions that this source contributes to
        ion_ids1_ids2 = []
        for total_source1d_ion in total_source1d_ions
            for source1d_ion in source1d.ion
                if total_source1d_ion.label == source1d_ion.label
                    push!(ion_ids1_ids2, (total_source1d_ion, source1d_ion))
                end
            end
        end

        # add to the tallies for this source
        x = source1d.grid.rho_tor_norm
        if rho == x
            rho = x
        end
        for (ids1, ids2) in [[(total_source1d, source1d), (total_source1d.electrons, source1d.electrons)]; ion_ids1_ids2]
            for field in keys(ids1)
                if (isempty(fields) || field ∈ fields || (field ∈ keys(matching) && matching[field] ∈ fields)) && field ∈ keys(matching)
                    if hasdata(ids2, field) || hasdata(ids2, matching[field])
                        y = getproperty(ids2, field)
                        if only_positive_negative != 0 && any((sign(yy) ∉ (0, sign(only_positive_negative)) for yy in y))
                            continue
                        end
                        if rho === x
                            rho_data = y
                        else
                            rho_data = DataInterpolations.LinearInterpolation(y, x).(rho)
                        end
                        setproperty!(ids1, field, getproperty(ids1, field) .+ rho_data)
                    end
                end
            end
        end
    end

    # ion.energy source is always zero
    for total_source1d_ion in total_source1d_ions
        empty!(total_source1d_ion, :energy)
    end

    return total_source1d
end

function total_radiation_sources(dd::IMAS.dd; time0::Float64=dd.global_time, kw...)
    return total_radiation_sources(dd.core_sources, dd.core_profiles.profiles_1d[time0]; kw...)
end

function total_radiation_sources(
    core_sources::IMAS.core_sources{T},
    cp1d::IMAS.core_profiles__profiles_1d{T};
    include_indexes::Vector{Int}=Int[],
    exclude_indexes::Vector{Int}=Int[]
) where {T<:Real}

    # we need to exclude the collisional_equipartition term
    index = IMAS.name_2_index(core_sources.source)[:collisional_equipartition]
    push!(exclude_indexes, index)

    fields = [:power_inside, :energy]
    only_positive_negative = -1
    return total_sources(core_sources, cp1d; include_indexes, exclude_indexes, fields, only_positive_negative)
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
        cs1d.electrons.energy = value = interp1d(range(0, 1, length(electrons_energy)), electrons_energy).(cs1d.grid.rho_tor_norm)
        cs1d.electrons.power_inside = cumtrapz(volume, value)
    elseif electrons_power_inside !== missing
        cs1d.electrons.power_inside = value = interp1d(range(0, 1, length(electrons_power_inside)), electrons_power_inside).(cs1d.grid.rho_tor_norm)
        cs1d.electrons.energy = gradient(volume, value)
    else
        cs1d.electrons.energy = zero(volume)
        cs1d.electrons.power_inside = zero(volume)
    end

    if total_ion_energy !== missing
        cs1d.total_ion_energy = value = interp1d(range(0, 1, length(total_ion_energy)), total_ion_energy).(cs1d.grid.rho_tor_norm)
        cs1d.total_ion_power_inside = cumtrapz(volume, value)
    elseif total_ion_power_inside !== missing
        cs1d.total_ion_power_inside = value = interp1d(range(0, 1, length(total_ion_power_inside)), total_ion_power_inside).(cs1d.grid.rho_tor_norm)
        cs1d.total_ion_energy = gradient(volume, value)
    else
        cs1d.total_ion_energy = zero(volume)
        cs1d.total_ion_power_inside = zero(volume)
    end

    if electrons_particles !== missing
        cs1d.electrons.particles = value = interp1d(range(0, 1, length(electrons_particles)), electrons_particles).(cs1d.grid.rho_tor_norm)
        cs1d.electrons.particles_inside = cumtrapz(volume, value)
    elseif electrons_particles_inside !== missing
        cs1d.electrons.particles_inside = value = interp1d(range(0, 1, length(electrons_particles_inside)), electrons_particles_inside).(cs1d.grid.rho_tor_norm)
        cs1d.electrons.particles = gradient(volume, value)
    else
        cs1d.electrons.particles = zero(volume)
        cs1d.electrons.particles_inside = zero(volume)
    end

    if j_parallel !== missing
        cs1d.j_parallel = value = interp1d(range(0, 1, length(j_parallel)), j_parallel).(cs1d.grid.rho_tor_norm)
        cs1d.current_parallel_inside = cumtrapz(volume, value)
    elseif current_parallel_inside !== missing
        cs1d.current_parallel_inside = value = interp1d(range(0, 1, length(current_parallel_inside)), current_parallel_inside).(cs1d.grid.rho_tor_norm)
        cs1d.j_parallel = gradient(volume, value)
    else
        cs1d.j_parallel = zero(volume)
        cs1d.current_parallel_inside = zero(volume)
    end

    if momentum_tor !== missing
        cs1d.momentum_tor = value = interp1d(range(0, 1, length(momentum_tor)), momentum_tor).(cs1d.grid.rho_tor_norm)
        cs1d.torque_tor_inside = cumtrapz(volume, value)
    elseif torque_tor_inside !== missing
        cs1d.torque_tor_inside = value = interp1d(range(0, 1, length(torque_tor_inside)), torque_tor_inside).(cs1d.grid.rho_tor_norm)
        cs1d.momentum_tor = gradient(volume, value)
    else
        cs1d.momentum_tor = zero(volume)
        cs1d.torque_tor_inside = zero(volume)
    end

    return source
end

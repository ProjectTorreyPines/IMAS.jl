document[Symbol("Physics sources")] = Symbol[]

"""
    fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles; DD_fusion::Bool=false)

Calculates fusion source from D-T and D-D reactions and adds them to `dd.core_sources`

If D+T plasma, then D+D is neglected

If D+D plasma fusion is included depending on `DD_fusion` switch
"""
function fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles; DD_fusion::Bool=false)
    cp1d = cp.profiles_1d[]
    ion_list = (ion.label for ion in cp1d.ion)
    if "T" in ion_list || "DT" in ion_list
        D_T_to_He4_source!(cs, cp; combine_DT=("DT" in ion_list))
    elseif DD_fusion
        D_D_to_He3_source!(cs, cp)
        D_D_to_T_source!(cs, cp)
    end
end

"""
    fusion_source!(dd::IMAS.dd; DD_fusion::Bool=false)
"""
function fusion_source!(dd::IMAS.dd; DD_fusion::Bool=false)
    return fusion_source!(dd.core_sources, dd.core_profiles; DD_fusion)
end

@compat public fusion_source!
push!(document[Symbol("Physics sources")], :fusion_source!)

"""
    collisional_exchange_source!(dd::IMAS.dd)

Calculates collisional exchange source and adds it to `dd.core_sources`
"""
function collisional_exchange_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature
    Ti = cp1d.t_i_average

    nu_exch = collision_frequencies(cp1d).nu_exch
    delta = 1.5 .* nu_exch .* ne .* mks.e .* (Te .- Ti)

    source = resize!(dd.core_sources.source, :collisional_equipartition; wipe=false)
    new_source(source, source.identifier.index, "exchange", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area;
        electrons_energy=-delta, total_ion_energy=delta)

    return source
end

@compat public collisional_exchange_source!
push!(document[Symbol("Physics sources")], :collisional_exchange_source!)

"""
    ohmic_source!(dd::IMAS.dd)

Calculates the ohmic source from data in `dd.core_profiles` and adds it to `dd.core_sources`
"""
function ohmic_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    j_ohmic = getproperty(cp1d, :j_ohmic, missing)
    if !ismissing(j_ohmic)
        eqt = dd.equilibrium.time_slice[]
        eqt1d = eqt.profiles_1d
        x = cp1d.grid.rho_tor_norm
        rho_eq = eqt1d.rho_tor_norm
        gm1_itp = cubic_interp1d(rho_eq, eqt1d.gm1)
        gm9_itp = cubic_interp1d(rho_eq, eqt1d.gm9)
        f_itp = cubic_interp1d(rho_eq, eqt1d.f)
        b0 = eqt.global_quantities.vacuum_toroidal_field.b0
        powerDensityOhm = @. (cp1d.j_tor * gm9_itp(x)) * (j_ohmic * b0) / (f_itp(x) * gm1_itp(x) * cp1d.conductivity_parallel)
        source = resize!(dd.core_sources.source, :ohmic; wipe=false)
        new_source(source, source.identifier.index, "ohmic", x, cp1d.grid.volume, cp1d.grid.area;
            electrons_energy=powerDensityOhm,
            j_parallel=j_ohmic)
        return source
    end
end

@compat public ohmic_source!
push!(document[Symbol("Physics sources")], :ohmic_source!)

"""
    bootstrap_source!(dd::IMAS.dd)

Calculates the bootsrap current from profiles in `dd.core_profiles`
and adds it both as source in `dd.core_sources` and `cp1d.j_bootstrap`
"""
function bootstrap_source!(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    j_bootstrap = Sauter_neo2021_bootstrap(eqt, cp1d)
    source = resize!(dd.core_sources.source, :bootstrap_current; wipe=false)
    new_source(source, source.identifier.index, "bootstrap", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area;
        j_parallel=j_bootstrap)
    cp1d.j_bootstrap = j_bootstrap
    return source
end

@compat public bootstrap_source!
push!(document[Symbol("Physics sources")], :bootstrap_source!)

"""
    radiation_source!(dd::IMAS.dd)

Calculates the total radiation by calling:

  - bremsstrahlung_source!()
  - line_radiation_source!()
  - synchrotron_source!()
"""
function radiation_source!(dd::IMAS.dd)
    bremsstrahlung_source!(dd) # electron energy

    line_radiation_source!(dd) # electron energy

    synchrotron_source!(dd) # electron energy (ion synchrotron not calculated)

    return dd
end

@compat public radiation_source!
push!(document[Symbol("Physics sources")], :radiation_source!)

"""
    sources!(dd::IMAS.dd; bootstrap::Bool=true, DD_fusion::Bool=false)

Calculates intrisic sources and sinks, and adds them to `dd.core_sources`
"""
function sources!(dd::IMAS.dd; bootstrap::Bool=true, DD_fusion::Bool=false)
    if bootstrap
        bootstrap_source!(dd) # current
    end

    ohmic_source!(dd) # electron energy, current

    collisional_exchange_source!(dd) # electron and ion energy

    radiation_source!(dd) # calls bremsstrahlung_source!() line_radiation_source!() synchrotron_source!()

    fusion_source!(dd; DD_fusion) # electron and ion energy, particles

    return nothing
end

@compat public sources!
push!(document[Symbol("Physics sources")], :sources!)

"""
    time_derivative_source!(dd::IMAS.dd, cp1d_old::IMAS.core_profiles__profiles_1d, Δt::Float64; zero_out::Bool)

Calculates time dependent sources and sinks, and adds them to `dd.core_sources`

These are the ∂/∂t term in the transport equations
"""
function time_derivative_source!(dd::IMAS.dd, cp1d_old::IMAS.core_profiles__profiles_1d, Δt::Float64; name::String="∂/∂t")
    cp1d = dd.core_profiles.profiles_1d[]
    eqt1d = dd.equilibrium.time_slice[].profiles_1d
    R_flux_avg = interp1d(eqt1d.rho_tor_norm, eqt1d.gm8).(cp1d.grid.rho_tor_norm)

    if Δt == 0.0
        electrons_particles = zero(cp1d.grid.rho_tor_norm)
        electrons_energy = zero(cp1d.grid.rho_tor_norm)
        total_ion_energy = zero(cp1d.grid.rho_tor_norm)
        momentum_tor = zero(cp1d.grid.rho_tor_norm)
    else
        electrons_particles = -(cp1d.electrons.density_thermal .- cp1d_old.electrons.density_thermal) / Δt

        electrons_energy = -1.5 * (pressure_thermal(cp1d.electrons) .- pressure_thermal(cp1d_old.electrons)) / Δt

        total_ion_energy = -1.5 * (pressure_thermal(cp1d.ion) .- pressure_thermal(cp1d_old.ion)) / Δt

        d_new = cp1d.rotation_frequency_tor_sonic .* total_mass_density(cp1d) .* R_flux_avg / Δt
        d_old = cp1d_old.rotation_frequency_tor_sonic .* total_mass_density(cp1d_old) .* R_flux_avg / Δt
        momentum_tor = d_new .- d_old
    end

    # fill in time derivative sources for electrons_particles, electrons_energy, total_ion_energy, momentum_tor
    source = resize!(dd.core_sources.source, :time_derivative, "identifier.name" => name; wipe=false)
    new_source(source, source.identifier.index, name, cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area;
        electrons_energy, total_ion_energy, electrons_particles, momentum_tor)

    # fill in time derivative sources for ion_particles
    cs1d = source.profiles_1d[]
    resize!(cs1d.ion, length(cp1d.ion); wipe=false)
    for (k, (ion, ion_old)) in enumerate(zip(cp1d.ion, cp1d_old.ion))
        cs1d.ion[k].label = ion.label
        fill!(cs1d.ion[k].element, ion.element)
        if Δt == 0.0 || ismissing(ion, :density_thermal)
            cs1d.ion[k].particles_inside = cs1d.ion[k].particles = zero(cp1d.grid.volume)
        else
            particles = -(ion.density_thermal .- ion_old.density_thermal) / Δt
            cs1d.ion[k].particles = particles
            cs1d.ion[k].particles_inside = cumtrapz(cp1d.grid.volume, particles)
        end
    end

    return source
end

"""
    time_derivative_source!(dd::IMAS.dd; zero_out::Bool=false)
"""
function time_derivative_source!(dd::IMAS.dd; name::String="∂/∂t", zero_out::Bool=false)
    if zero_out
        time_derivative_source!(dd, dd.core_profiles.profiles_1d[], 0.0; name)
    else
        cp_time = dd.core_profiles.time
        i = nearest_causal_time(cp_time, dd.global_time).index
        if i > 1 && !isinf(cp_time[i-1])
            time_derivative_source!(dd, dd.core_profiles.profiles_1d[i-1], cp_time[i] - cp_time[i-1]; name)
        end
    end
end

@compat public time_derivative_source!
push!(document[Symbol("Physics sources")], :time_derivative_source!)

"""
    total_mass_density(cp1d::IMAS.core_profiles__profiles_1d)

Finds the total mass density [kg/m^-3]
"""
function total_mass_density(cp1d::IMAS.core_profiles__profiles_1d)
    mass_density = mks.m_e * cp1d.electrons.density_thermal
    if hasdata(cp1d.electrons, :density_fast)
        mass_density .+= mks.m_e * cp1d.electrons.density_fast
    end
    for ion in cp1d.ion
        for field in (:density_thermal, :density_fast)
            if hasdata(ion, field)
                mass_density .+= getproperty(ion, field) * ion.element[1].a * mks.m_p
            end
        end
    end
    return mass_density
end

@compat public total_mass_density
push!(document[Symbol("Physics sources")], :total_mass_density)

"""
    total_power_inside(
        core_sources::IMAS.core_sources,
        cp1d::IMAS.core_profiles__profiles_1d;
        time0::Float64=global_time(cp1d),
        include_radiation::Bool=true,
        include_time_derivative::Bool=true
    )

Returns total power inside of the separatrix
"""
function total_power_inside(
    core_sources::IMAS.core_sources,
    cp1d::IMAS.core_profiles__profiles_1d;
    time0::Float64=global_time(cp1d),
    include_radiation::Bool=true,
    include_time_derivative::Bool=true
)
    exclude_indexes = Int[]
    if !include_time_derivative
        exclude_indexes = [11]
    else
        exclude_indexes = Int[]
    end
    total_source = total_sources(core_sources, cp1d; time0, fields=[:power_inside, :total_ion_power_inside], exclude_indexes)
    tot_pow_in = total_source.electrons.power_inside[end] + total_source.total_ion_power_inside[end]
    if !include_radiation
        # NOTE: this ignore radiation, since add back the radiation losses that were already subtracted
        tot_pow_in -= radiation_losses(core_sources; time0)
    end
    return tot_pow_in
end

@compat public total_power_inside
push!(document[Symbol("Physics sources")], :total_power_inside)

"""
    total_power_source(source::IMAS.core_sources__source___profiles_1d)

Returns the total power (electron + ion) for a single source
"""
function total_power_source(source::IMAS.core_sources__source___profiles_1d)
    return getproperty(source.electrons, :power_inside, [0.0])[end] + getproperty(source, :total_ion_power_inside, [0.0])[end]
end

@compat public total_power_source
push!(document[Symbol("Physics sources")], :total_power_source)

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

@compat public total_power_time
push!(document[Symbol("Physics sources")], :total_power_time)

"""
    retain_source(source::IMAS.core_sources__source, all_indexes::Vector{Int}, include_indexes::Vector{Int}, exclude_indexes::Vector{Int})

Function that decides whether a source should be kept or ignored when totaling sources
"""
function retain_source(source::IMAS.core_sources__source, all_indexes::Vector{Int}, include_indexes::Vector{Int}, exclude_indexes::Vector{Int})
    index = source.identifier.index
    if index ∈ include_indexes
        return true
    elseif !isempty(include_indexes)
        return false
    elseif index == 0
        @debug "total_sources() skipping unspecified source with index $index"
        return false
    elseif index == 1 && any(all_indexes .> 1)
        @debug "total_sources() skipping total source with index $index"
        return false
    elseif 107 >= index >= 100 && any(all_indexes .< 5)
        @debug "total_sources() skipping combination source with index $index"
        return false
    elseif index == 200 && any(300 .> all_indexes .> 200)
        @debug "total_sources() skipping total radiation source with index $index"
        return false
    elseif index ∈ exclude_indexes
        return false
    end
    return true
end

@compat public retain_source
push!(document[Symbol("Physics sources")], :retain_source)

const _core_sources_integral_value_keys = Dict{Symbol,Symbol}()
_core_sources_integral_value_keys[:power_inside] = :energy
_core_sources_integral_value_keys[:energy] = :power_inside
_core_sources_integral_value_keys[:total_ion_power_inside] = :total_ion_energy
_core_sources_integral_value_keys[:total_ion_energy] = :total_ion_power_inside
_core_sources_integral_value_keys[:particles_inside] = :particles
_core_sources_integral_value_keys[:particles] = :particles_inside
_core_sources_integral_value_keys[:current_parallel_inside] = :j_parallel
_core_sources_integral_value_keys[:j_parallel] = :current_parallel_inside
_core_sources_integral_value_keys[:torque_tor_inside] = :momentum_tor
_core_sources_integral_value_keys[:momentum_tor] = :torque_tor_inside

const _core_sources_value_keys = Dict{Symbol,Symbol}()
_core_sources_value_keys[:energy] = :power_inside
_core_sources_value_keys[:total_ion_energy] = :total_ion_power_inside
_core_sources_value_keys[:particles] = :particles_inside
_core_sources_value_keys[:j_parallel] = :current_parallel_inside
_core_sources_value_keys[:momentum_tor] = :torque_tor_inside

const _core_sources_integral_keys = Dict{Symbol,Symbol}()
_core_sources_integral_keys[:power_inside] = :energy
_core_sources_integral_keys[:total_ion_power_inside] = :total_ion_energy
_core_sources_integral_keys[:particles_inside] = :particles
_core_sources_integral_keys[:current_parallel_inside] = :j_parallel
_core_sources_integral_keys[:torque_tor_inside] = :momentum_tor

"""
    total_sources(
        core_sources::IMAS.core_sources{T},
        cp1d::IMAS.core_profiles__profiles_1d{T};
        time0::Float64;
        include_indexes::Vector{Int}=Int[],
        exclude_indexes::Vector{Int}=Int[],
        fields::Vector{Symbol}=Symbol[],
        only_positive_negative::Int=0) where {T<:Real}

Returns core_sources__source___profiles_1d with sources totals and possiblity to

  - include/exclude certain sources based on their unique index identifier
  - include only certain fields among these: [:particles_inside, :energy, :power_inside, :momentum_tor, :total_ion_power_inside, :total_ion_energy, :j_parallel, :torque_tor_inside, :current_parallel_inside, :particles]
"""
function total_sources(
    core_sources::IMAS.core_sources{T},
    cp1d::IMAS.core_profiles__profiles_1d{T};
    time0::Float64,
    include_indexes::Vector{Int}=Int[],
    exclude_indexes::Vector{Int}=Int[],
    fields::Vector{Symbol}=Symbol[],
    only_positive_negative::Int=0) where {T<:Real}

    total_source1d = core_sources__source___profiles_1d{T}()
    return total_sources!(total_source1d, core_sources, cp1d; time0, include_indexes, exclude_indexes, fields, only_positive_negative)
end

function total_sources!(
    total_source1d::IMAS.core_sources__source___profiles_1d{T},
    core_sources::IMAS.core_sources{T},
    cp1d::IMAS.core_profiles__profiles_1d{T};
    time0::Float64,
    include_indexes::Vector{Int}=Int[],
    exclude_indexes::Vector{Int}=Int[],
    fields::Vector{Symbol}=Symbol[],
    only_positive_negative::Int=0) where {T<:Real}

    total_source1d.grid.rho_tor_norm = rho = cp1d.grid.rho_tor_norm
    total_source1d.time = time0

    @assert isempty(fields) || all(field in keys(_core_sources_integral_value_keys) for field in fields) "Supported fields are $(collect(keys(_core_sources_integral_value_keys)))"

    for prop in (:volume, :area, :surface)
        value = getproperty(cp1d.grid, prop)
        if hasdata(total_source1d.grid, prop)
            getproperty(total_source1d.grid, prop) .= value
        else
            setproperty!(total_source1d.grid, prop, value)
        end
    end

    # initialize ions (get list of ions by looking both under core_profiles and core_sources)
    for ion in cp1d.ion
        resize!(total_source1d.ion, "element[1].a" => ion.element[1].z_n, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label; wipe=false)
    end
    for source in core_sources.source
        if !isempty(source.profiles_1d) && source.profiles_1d[1].time <= time0
            source1d = source.profiles_1d[1]
            for ion in source1d.ion
                resize!(total_source1d.ion, "element[1].a" => ion.element[1].z_n, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label; wipe=false)
            end
        end
    end

    # zero out total_sources
    for ids1 in [[total_source1d, total_source1d.electrons]; total_source1d.ion]
        for field in keys(ids1)
            if field in keys(_core_sources_integral_value_keys)
                if hasdata(ids1, field) && length(getfield(ids1, field)) == length(rho)
                    fill!(getproperty(ids1, field), zero(T))
                else
                    # assume coordinates have been set
                    setproperty!(ids1, field, zeros(T, size(rho)); error_on_missing_coordinates=false)
                end
            end
        end
    end

    # start accumulating
    all_indexes = [source.identifier.index for source in core_sources.source]
    for source in core_sources.source
        if !retain_source(source, all_indexes, include_indexes, exclude_indexes)
            continue
        end
        if isempty(source.profiles_1d)
            continue # skip sources that have no profiles_1d time slices
        end
        if !isempty(source.profiles_1d) && source.profiles_1d[1].time > time0
            continue # skip sources that start after time of interest
        end

        source1d = source.profiles_1d[time0]

        if ismissing(source1d.grid, :rho_tor_norm)
            continue # skip sources don't have radial coordinate, since they cannot have data
        end

        # ions that this source contributes to
        ids1_ids2 = [(total_source1d, source1d), (total_source1d.electrons, source1d.electrons)]
        for total_source1d_ion in total_source1d.ion
            for source1d_ion in source1d.ion
                if total_source1d_ion.label == source1d_ion.label
                    push!(ids1_ids2, (total_source1d_ion, source1d_ion))
                end
            end
        end

        # add to the tallies for this source
        x = source1d.grid.rho_tor_norm
        if rho == x
            rho = x
        end
        for (ids1, ids2) in ids1_ids2
            for field in keys(ids1)
                if (isempty(fields) || field ∈ fields || (field ∈ keys(_core_sources_integral_value_keys) && _core_sources_integral_value_keys[field] ∈ fields)) &&
                   field ∈ keys(_core_sources_integral_value_keys)
                    if hasdata(ids2, field) || hasdata(ids2, _core_sources_integral_value_keys[field])
                        y = getproperty(ids2, field)
                        if only_positive_negative != 0 && any((sign(yy) ∉ (0, sign(only_positive_negative)) for yy in y))
                            continue
                        end
                        if rho === x
                            interpolated_data = y
                        else
                            interpolated_data = DataInterpolations.LinearInterpolation(y, x; extrapolation=DataInterpolations.ExtrapolationType.Constant).(rho)
                        end
                        if hasdata(ids1, field)
                            getproperty(ids1, field) .= getproperty(ids1, field) .+ interpolated_data
                        else
                            setproperty!(ids1, field, getproperty(ids1, field) .+ interpolated_data)
                        end
                    end
                end
            end
        end
    end

    return total_source1d
end

"""
    total_sources(dd::IMAS.dd; time0::Float64=dd.global_time, kw...)
"""
function total_sources(dd::IMAS.dd; time0::Float64=dd.global_time, kw...)
    return total_sources(dd.core_sources, dd.core_profiles.profiles_1d[time0]; time0, kw...)
end

@compat public total_sources
push!(document[Symbol("Physics sources")], :retain_source)

function total_radiation_sources(
    core_sources::IMAS.core_sources{T},
    cp1d::IMAS.core_profiles__profiles_1d{T};
    time0::Float64,
    include_indexes::Vector{Int}=index_radiation_sources,
    exclude_indexes::Vector{Int}=Int[]) where {T<:Real}

    fields = [:power_inside, :energy]
    return total_sources(core_sources, cp1d; time0, include_indexes, exclude_indexes, fields)
end

"""
    total_radiation_sources(dd::IMAS.dd; time0::Float64=dd.global_time, kw...)
"""
function total_radiation_sources(dd::IMAS.dd; time0::Float64=dd.global_time, kw...)
    return total_radiation_sources(dd.core_sources, dd.core_profiles.profiles_1d[time0]; time0, kw...)
end

@compat public total_radiation_sources
push!(document[Symbol("Physics sources")], :total_radiation_sources)

"""
    sawteeth_source!(dd::IMAS.dd; qmin_desired::Float64=1.0)

Model sawteeth by flattening all sources where abs(q) drops below qmin_desired
"""
function sawteeth_source!(dd::IMAS.dd; qmin_desired::Float64=1.0)
    eqt1d = dd.equilibrium.time_slice[].profiles_1d

    q = abs.(eqt1d.q)
    if !any(x -> x < qmin_desired, q)
        rho0 = 0.0
    else
        rho0 = eqt1d.rho_tor_norm[findlast(qq -> abs(qq) < qmin_desired, q)]
    end

    return sawteeth_source!(dd, rho0)
end

function sawteeth_source!(dd::IMAS.dd{T}, ::Nothing) where {T<:Real}
    return sawteeth_source!(dd, 0.0)
end

function sawteeth_source!(dd::IMAS.dd{T}, i_qdes::Int) where {T<:Real}
    cp1d = dd.core_profiles.profiles_1d[]
    return sawteeth_source!(dd, cp1d.grid.rho_tor_norm[i_qdes])
end

"""
    sawteeth_source!(dd::IMAS.dd{T}, rho0::T) where {T<:Real}

Model sawteeth by flattening all sources within the inversion radius rho0
"""
function sawteeth_source!(dd::IMAS.dd{T}, rho0::T) where {T<:Real}
    @assert rho0 <= 1.0
    cp1d = dd.core_profiles.profiles_1d[]

    # fill in sawteeth
    source = resize!(dd.core_sources.source, :sawteeth, "identifier.name" => "sawteeth"; wipe=false)

    # get past value of sawteeth source
    if isempty(source.profiles_1d) || hasdata(source.profiles_1d[])
        old_source1d = total_sources(dd.core_sources, cp1d; time0=dd.global_time, include_indexes=[-10000])
    else
        source1d = source.profiles_1d[]
        itime = IMAS.index(source1d)
        if itime == 1
            old_source1d = deepcopy(source1d)
        else
            old_source1d = source.profiles_1d[itime-1]
        end
    end

    # wipe data at current time
    source1d = resize!(source.profiles_1d)

    # identify sawteeth inversion radius
    if rho0 > 0.0
        # exlude :time_dependent source (701)
        total_source1d = total_sources(dd.core_sources, cp1d; time0=dd.global_time, exclude_indexes=Int[701])
        fill!(source1d, total_source1d)
    else
        # this will return an empty source
        total_source1d = total_sources(dd.core_sources, cp1d; time0=dd.global_time, include_indexes=[-10000])
        fill!(source1d, total_source1d)
    end

    width = min(rho0 / 4, 0.05)

    α = 0.5

    # sawteeth source as difference between the total using the flattened profiles and the total using the original profiles
    for (leaf, old_leaf) in zip(IMASdd.AbstractTrees.Leaves(source1d), IMASdd.AbstractTrees.Leaves(old_source1d))
        @assert leaf.field == old_leaf.field
        if leaf.field in keys(_core_sources_value_keys)
            if leaf.field == :j_parallel
                leaf.value .= (flatten_profile!(copy(leaf.value), source1d.grid.rho_tor_norm, source1d.grid.area, rho0, width) .- leaf.value) * α .+ old_leaf.value .* (1.0 .- α)
            else
                leaf.value .= (flatten_profile!(copy(leaf.value), source1d.grid.rho_tor_norm, source1d.grid.volume, rho0, width) .- leaf.value) * α .+ old_leaf.value .* (1.0 .- α)
            end
        end
    end

    # calculate integrated quantities
    for leaf in IMASdd.AbstractTrees.Leaves(source1d)
        if leaf.field in keys(_core_sources_value_keys)
            if leaf.field == :j_parallel
                value = cumtrapz(source1d.grid.area, leaf.value)
            else
                value = cumtrapz(source1d.grid.volume, leaf.value)
            end
            setproperty!(leaf.ids, _core_sources_value_keys[leaf.field], value)
        end
    end

    return source
end

@compat public sawteeth_source!
push!(document[Symbol("Physics sources")], :sawteeth_source!)

"""
    sawteeth_profiles!(cp1d::IMAS.core_profiles__profiles_1d{T}, rho0::T) where {T<:Real}

Model sawteeth by flattening core_profiles within the inversion radius rho0
"""
function sawteeth_profiles!(cp1d::IMAS.core_profiles__profiles_1d{T}, rho0::T) where {T<:Real}

    @assert rho0 <= 1.0
    width = min(rho0 / 4, 0.05)

    α = 0.5

    for leaf in IMAS.IMASdd.AbstractTrees.Leaves(cp1d)
        if IMAS.hasdata(leaf.ids, leaf.field) && leaf.field in [:temperature, :density_thermal, :rotation_frequency_tor_sonic]
            old_leaf_value = leaf.value
            new_leaf_value = IMAS.flatten_profile!(copy(leaf.value), cp1d.grid.rho_tor_norm, cp1d.grid.volume, rho0, width) * α .+ old_leaf_value .* (1.0 .- α)
            setproperty!(leaf.ids, leaf.field, new_leaf_value)
        end
    end

    return cp1d
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
    csglbl = resize!(source.global_quantities)
    cs1d.grid.rho_tor_norm = rho
    cs1d.grid.volume = volume
    cs1d.grid.area = area

    electrons = cs1d.electrons

    if electrons_energy !== missing
        value = electrons_energy
        setproperty!(electrons, :energy, value; error_on_missing_coordinates=false)
        setproperty!(electrons, :power_inside, cumtrapz(volume, value); error_on_missing_coordinates=false)
    elseif electrons_power_inside !== missing
        value = electrons_power_inside
        setproperty!(electrons, :power_inside, value; error_on_missing_coordinates=false)
        setproperty!(electrons, :energy, gradient(volume, value); error_on_missing_coordinates=false)
    else
        setproperty!(electrons, :energy, zero(volume); error_on_missing_coordinates=false)
        setproperty!(electrons, :power_inside, zero(volume); error_on_missing_coordinates=false)
    end
    csglbl.electrons.power = electrons.power_inside[end]

    if total_ion_energy !== missing
        value = total_ion_energy
        setproperty!(cs1d, :total_ion_energy, value; error_on_missing_coordinates=false)
        setproperty!(cs1d, :total_ion_power_inside, cumtrapz(volume, value); error_on_missing_coordinates=false)
    elseif total_ion_power_inside !== missing
        value = total_ion_power_inside
        setproperty!(cs1d, :total_ion_power_inside, value; error_on_missing_coordinates=false)
        setproperty!(cs1d, :total_ion_energy, gradient(volume, value); error_on_missing_coordinates=false)
    else
        setproperty!(cs1d, :total_ion_energy, zero(volume); error_on_missing_coordinates=false)
        setproperty!(cs1d, :total_ion_power_inside, zero(volume); error_on_missing_coordinates=false)
    end
    csglbl.total_ion_power = cs1d.total_ion_power_inside[end]
    csglbl.power = csglbl.total_ion_power + csglbl.electrons.power

    if electrons_particles !== missing
        value = electrons_particles
        setproperty!(electrons, :particles, value; error_on_missing_coordinates=false)
        setproperty!(electrons, :particles_inside, cumtrapz(volume, value); error_on_missing_coordinates=false)
    elseif electrons_particles_inside !== missing
        value = electrons_particles_inside
        setproperty!(electrons, :particles_inside, value; error_on_missing_coordinates=false)
        setproperty!(electrons, :particles, gradient(volume, value); error_on_missing_coordinates=false)
    else
        setproperty!(electrons, :particles, zero(volume); error_on_missing_coordinates=false)
        setproperty!(electrons, :particles_inside, zero(volume); error_on_missing_coordinates=false)
    end
    csglbl.electrons.particles = electrons.particles_inside[end]

    if j_parallel !== missing
        value = j_parallel
        setproperty!(cs1d, :j_parallel, value; error_on_missing_coordinates=false)
        setproperty!(cs1d, :current_parallel_inside, cumtrapz(area, value); error_on_missing_coordinates=false)
    elseif current_parallel_inside !== missing
        value = current_parallel_inside
        setproperty!(cs1d, :current_parallel_inside, value; error_on_missing_coordinates=false)
        setproperty!(cs1d, :j_parallel, gradient(area, value); error_on_missing_coordinates=false)
    else
        setproperty!(cs1d, :j_parallel, zero(area); error_on_missing_coordinates=false)
        setproperty!(cs1d, :current_parallel_inside, zero(area); error_on_missing_coordinates=false)
    end
    csglbl.current_parallel = cs1d.current_parallel_inside[end]

    if momentum_tor !== missing
        value = momentum_tor
        setproperty!(cs1d, :momentum_tor, value; error_on_missing_coordinates=false)
        setproperty!(cs1d, :torque_tor_inside, cumtrapz(volume, value); error_on_missing_coordinates=false)
    elseif torque_tor_inside !== missing
        value = torque_tor_inside
        setproperty!(cs1d, :torque_tor_inside, value; error_on_missing_coordinates=false)
        setproperty!(cs1d, :momentum_tor, gradient(volume, value); error_on_missing_coordinates=false)
    else
        setproperty!(cs1d, :momentum_tor, zero(volume); error_on_missing_coordinates=false)
        setproperty!(cs1d, :torque_tor_inside, zero(volume); error_on_missing_coordinates=false)
    end
    csglbl.torque_tor = cs1d.torque_tor_inside[end]

    return source
end

@compat public new_source
push!(document[Symbol("Physics sources")], :new_source)

"""
    total_power(
        ps::Union{IMAS.pulse_schedule,IMAS.pulse_schedule__ec,IMAS.pulse_schedule__ic,IMAS.pulse_schedule__lh,IMAS.pulse_schedule__nbi},
        times::AbstractVector{Float64};
        time_smooth::Float64)

Total injected power interpolated on a given time basis
"""
function total_power(
    ps::Union{IMAS.pulse_schedule,IMAS.pulse_schedule__ec,IMAS.pulse_schedule__ic,IMAS.pulse_schedule__lh,IMAS.pulse_schedule__nbi},
    times::AbstractVector{Float64};
    tau_smooth::Float64
)
    datas = zero(times)
    for leaf in leaves(ps)
        if contains(location(leaf.ids), ".power") && hasdata(leaf.ids, leaf.field)
            time = getproperty(coordinates(leaf.ids, leaf.field)[1])
            data = getproperty(leaf.ids, leaf.field)
            sbp_itp = interp1d(time, smooth_beam_power(time, data, tau_smooth))
            @. datas += sbp_itp(times)
        end
    end
    return datas
end

@compat public total_power
push!(document[Symbol("Physics sources")], :total_power)
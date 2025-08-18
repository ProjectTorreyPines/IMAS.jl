import IMASutils: trapz, cumtrapz

function IMASdd.get_expressions(::Type{Val{:dynamic}})
    return dynamic_expressions
end

const dynamic_expressions = Dict{String,Function}()
dyexp = dynamic_expressions

#= =========== =#
# core_profiles #
#= =========== =#

#  core_profiles.profiles_1d  #

dyexp["core_profiles.profiles_1d[:].electrons.density"] =
    (; profiles_1d, electrons, _...) -> begin
        if !hasdata(electrons, :density_thermal) && !hasdata(electrons, :density_fast)
            rho_tor_norm = profiles_1d.grid.rho_tor_norm
            return zero(rho_tor_norm)
        elseif !hasdata(electrons, :density_thermal) && hasdata(electrons, :density_fast)
            return electrons.density_fast
        elseif hasdata(electrons, :density_thermal) && !hasdata(electrons, :density_fast)
            return electrons.density_thermal
        else
            return electrons.density_thermal .+ electrons.density_fast
        end
    end

dyexp["core_profiles.profiles_1d[:].electrons.pressure_thermal"] =
    (; electrons, _...) -> pressure_thermal(electrons)

dyexp["core_profiles.profiles_1d[:].electrons.pressure"] =
    (; electrons, _...) -> electrons.pressure_thermal .+ electrons.pressure_fast_parallel .+ electrons.pressure_fast_perpendicular .* 2.0


dyexp["core_profiles.profiles_1d[:].ion[:].z_ion"] =
    (; ion, _...) -> begin
        if length(ion.element) == 1
            return ion.element[1].z_n
        else
            error("z_ion expression does not yet handle multiple charge states")
        end
    end

dyexp["core_profiles.profiles_1d[:].t_i_average"] =
    (; profiles_1d, _...) -> t_i_average(profiles_1d)

dyexp["core_profiles.profiles_1d[:].ion[:].density"] =
    (; profiles_1d, ion, _...) -> begin
        if !hasdata(ion, :density_thermal) && !hasdata(ion, :density_fast)
            rho_tor_norm = profiles_1d.grid.rho_tor_norm
            return zero(rho_tor_norm)
        elseif !hasdata(ion, :density_thermal) && hasdata(ion, :density_fast)
            return ion.density_fast
        elseif hasdata(ion, :density_thermal) && !hasdata(ion, :density_fast)
            return ion.density_thermal
        else
            return ion.density_thermal .+ ion.density_fast
        end
    end

dyexp["core_profiles.profiles_1d[:].ion[:].density_fast"] =
    (; ion, _...) -> ion.density .- ion.density_thermal

dyexp["core_profiles.profiles_1d[:].ion[:].density_thermal"] =
    (; ion, _...) -> ion.density .- ion.density_fast

dyexp["core_profiles.profiles_1d[:].ion[:].pressure_thermal"] =
    (; ion, _...) -> pressure_thermal(ion)

dyexp["core_profiles.profiles_1d[:].ion[:].pressure"] =
    (; ion, _...) -> ion.pressure_thermal .+ ion.pressure_fast_parallel .+ ion.pressure_fast_perpendicular .* 2.0

dyexp["core_profiles.profiles_1d[:].pressure_ion_total"] =
    (; profiles_1d, _...) -> pressure_thermal(profiles_1d.ion)

dyexp["core_profiles.profiles_1d[:].pressure_thermal"] =
    (; profiles_1d, _...) -> pressure_thermal(profiles_1d)

dyexp["core_profiles.profiles_1d[:].pressure_parallel"] =
    (; profiles_1d, _...) -> begin
        ppar = profiles_1d.pressure_thermal ./ 3.0 .+ profiles_1d.electrons.pressure_fast_parallel
        for ion in profiles_1d.ion
            ppar .+= ion.pressure_fast_parallel
        end
        return ppar
    end

dyexp["core_profiles.profiles_1d[:].pressure_perpendicular"] =
    (; profiles_1d, _...) -> begin
        pperp = profiles_1d.pressure_thermal ./ 3.0 .+ profiles_1d.electrons.pressure_fast_perpendicular
        for ion in profiles_1d.ion
            pperp .+= ion.pressure_fast_perpendicular
        end
        return pperp
    end

dyexp["core_profiles.profiles_1d[:].pressure"] =
    (; profiles_1d, _...) -> profiles_1d.pressure_perpendicular .* 2.0 .+ profiles_1d.pressure_parallel


dyexp["core_profiles.profiles_1d[:].conductivity_parallel"] =
    (; dd, profiles_1d, _...) -> neo_conductivity(dd.equilibrium.time_slice[profiles_1d.time], profiles_1d)

dyexp["core_profiles.profiles_1d[:].j_bootstrap"] =
    (; dd, profiles_1d, _...) -> findfirst(:bootstrap_current, dd.core_sources.source).profiles_1d[profiles_1d.time].j_parallel

dyexp["core_profiles.profiles_1d[:].j_ohmic"] =
    (; profiles_1d, _...) -> profiles_1d.j_total .- profiles_1d.j_non_inductive

dyexp["core_profiles.profiles_1d[:].j_non_inductive"] =
    (; dd, profiles_1d, _...) -> total_sources(dd.core_sources, profiles_1d; time0=profiles_1d.time, exclude_indexes=[7, 409, 701], fields=[:j_parallel]).j_parallel # no ohmic, sawteeth, or time_depedent

dyexp["core_profiles.profiles_1d[:].j_total"] =
    (; dd, profiles_1d, _...) -> begin
        if hasdata(profiles_1d, :j_ohmic)
            return profiles_1d.j_ohmic .+ profiles_1d.j_non_inductive
        elseif !isempty(dd.equilibrium.time) && profiles_1d.time >= dd.equilibrium.time[1]
            rho_tor_norm = profiles_1d.grid.rho_tor_norm
            eqt1d = dd.equilibrium.time_slice[profiles_1d.time].profiles_1d
            return interp1d(eqt1d.rho_tor_norm, eqt1d.j_parallel).(rho_tor_norm)
        end
    end

dyexp["core_profiles.profiles_1d[:].j_tor"] =
    (; dd, profiles_1d, _...) -> begin
        rho_tor_norm = profiles_1d.grid.rho_tor_norm
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        Jpar_2_Jtor(rho_tor_norm, profiles_1d.j_total, true, eqt)
    end

dyexp["core_profiles.profiles_1d[:].q"] =
    (; dd, profiles_1d, _...) -> begin
        rho_tor_norm = profiles_1d.grid.rho_tor_norm
        eqt1d = dd.equilibrium.time_slice[profiles_1d.time].profiles_1d
        return interp1d(eqt1d.rho_tor_norm, eqt1d.q).(rho_tor_norm)
    end

dyexp["core_profiles.profiles_1d[:].zeff"] =
    (; dd, profiles_1d, _...) -> zeff(profiles_1d)

dyexp["core_profiles.profiles_1d[:].rotation_frequency_tor_sonic"] =
    (; dd, profiles_1d, _...) -> ωtor2sonic(profiles_1d)

dyexp["core_profiles.profiles_1d[:].ion[:].rotation_frequency_tor"] =
    (; dd, profiles_1d, ion, _...) -> sonic2ωtor(profiles_1d, ion)

#  core_profiles.global_quantities  #
dyexp["core_profiles.global_quantities.current_non_inductive"] =
    (; dd, core_profiles, _...) -> [Ip_non_inductive(cp1d, dd.equilibrium.time_slice[cp1d.time]) for cp1d in core_profiles.profiles_1d]

dyexp["core_profiles.global_quantities.current_bootstrap"] =
    (; dd, core_profiles, _...) -> [Ip_bootstrap(cp1d, dd.equilibrium.time_slice[cp1d.time]) for cp1d in core_profiles.profiles_1d]

dyexp["core_profiles.global_quantities.current_ohmic"] =
    (; dd, core_profiles, _...) -> [Ip_ohmic(cp1d, dd.equilibrium.time_slice[cp1d.time]) for cp1d in core_profiles.profiles_1d]

dyexp["core_profiles.global_quantities.ip"] =
    (; dd, core_profiles, _...) -> [Ip(cp1d, dd.equilibrium.time_slice[cp1d.time]) for cp1d in core_profiles.profiles_1d]

dyexp["core_profiles.global_quantities.beta_tor_norm"] =
    (; dd, core_profiles, _...) -> [beta_tor_norm(dd.equilibrium.time_slice[cp1d.time], cp1d) for cp1d in core_profiles.profiles_1d]

dyexp["core_profiles.global_quantities.v_loop"] =
    (; dd, _...) -> [v_loop(cp1d) for cp1d in core_profiles.profiles_1d]

#= ========= =#
# equilibrium #
#= ========= =#
dyexp["equilibrium.time_slice[:].global_quantities.energy_mhd"] =
    (; time_slice, _...) -> 3 / 2 * trapz(time_slice.profiles_1d.volume, time_slice.profiles_1d.pressure)

dyexp["equilibrium.time_slice[:].global_quantities.q_95"] =
    (; time_slice, _...) -> Interpolations.linear_interpolation(norm01(time_slice.profiles_1d.psi), time_slice.profiles_1d.q)(0.95)

dyexp["equilibrium.time_slice[:].global_quantities.q_axis"] =
    (; time_slice, _...) -> time_slice.profiles_1d.q[1]

dyexp["equilibrium.time_slice[:].global_quantities.q_min"] =
    (; time_slice, _...) -> minimum(time_slice.profiles_1d.q)

dyexp["equilibrium.time_slice[:].global_quantities.psi_axis"] =
    (; time_slice, _...) -> time_slice.profiles_1d.psi[1]

dyexp["equilibrium.time_slice[:].global_quantities.psi_boundary"] =
    (; time_slice, _...) -> time_slice.profiles_1d.psi[end]

dyexp["equilibrium.time_slice[:].global_quantities.magnetic_axis.r"] =
    (; time_slice, _...) -> time_slice.profiles_1d.geometric_axis.r[1]

dyexp["equilibrium.time_slice[:].global_quantities.magnetic_axis.z"] =
    (; time_slice, _...) -> time_slice.profiles_1d.geometric_axis.z[1]


dyexp["equilibrium.time_slice[:].global_quantities.vacuum_toroidal_field.b0"] =
    (; dd, equilibrium, time_slice, _...) -> begin
        if hasdata(dd.pulse_schedule)
            return get_time_array(dd.pulse_schedule.tf.b_field_tor_vacuum, :reference, time_slice.time, :linear)
        else
            return get_time_array(equilibrium.vacuum_toroidal_field, :b0, time_slice.time, :constant)
        end
    end

dyexp["equilibrium.time_slice[:].global_quantities.vacuum_toroidal_field.r0"] =
    (; dd, _...) -> dd.pulse_schedule.tf.r0

dyexp["equilibrium.time_slice[:].global_quantities.magnetic_axis.b_field_tor"] =
    (; time_slice, _...) -> time_slice.global_quantities.vacuum_toroidal_field.b0 * time_slice.global_quantities.vacuum_toroidal_field.r0 / time_slice.boundary.geometric_axis.r

dyexp["equilibrium.time_slice[:].profiles_1d.geometric_axis.r"] =
    (; time_slice, _...) -> (time_slice.profiles_1d.r_outboard .+ time_slice.profiles_1d.r_inboard) .* 0.5

dyexp["equilibrium.time_slice[:].profiles_1d.geometric_axis.z"] =
    (; time_slice, _...) -> zero(time_slice.profiles_1d.psi) .+ time_slice.global_quantities.magnetic_axis.z

dyexp["equilibrium.time_slice[:].boundary.geometric_axis.r"] =
    (; time_slice, _...) -> begin
        if !ismissing(time_slice.profiles_1d.geometric_axis, :r)
            return time_slice.profiles_1d.geometric_axis.r[end]
        else
            minR, maxR = extrema(time_slice.boundary.outline.r)
            return (minR + maxR) / 2
        end
    end

dyexp["equilibrium.time_slice[:].boundary.geometric_axis.z"] =
    (; time_slice, _...) -> begin
        if !ismissing(time_slice.profiles_1d.geometric_axis, :z)
            return time_slice.profiles_1d.geometric_axis.z[end]
        else
            minZ, maxZ = extrema(time_slice.boundary.outline.z)
            return (minZ + maxZ) / 2
        end
    end


dyexp["equilibrium.time_slice[:].boundary.minor_radius"] =
    (; time_slice, _...) -> (time_slice.profiles_1d.r_outboard[end] - time_slice.profiles_1d.r_inboard[end]) * 0.5

dyexp["equilibrium.time_slice[:].boundary.elongation"] =
    (; time_slice, _...) -> (time_slice.boundary.elongation_lower + time_slice.boundary.elongation_upper) * 0.5

dyexp["equilibrium.time_slice[:].boundary.elongation_lower"] =
    (; time_slice, _...) -> time_slice.profiles_1d.elongation[end] # <======= IMAS limitation

dyexp["equilibrium.time_slice[:].boundary.elongation_upper"] =
    (; time_slice, _...) -> time_slice.profiles_1d.elongation[end] # <======= IMAS limitation

dyexp["equilibrium.time_slice[:].boundary.triangularity"] =
    (; time_slice, _...) -> (time_slice.boundary.triangularity_lower + time_slice.boundary.triangularity_upper) * 0.5

dyexp["equilibrium.time_slice[:].boundary.triangularity_lower"] =
    (; time_slice, _...) -> time_slice.profiles_1d.triangularity_lower[end]

dyexp["equilibrium.time_slice[:].boundary.triangularity_upper"] =
    (; time_slice, _...) -> time_slice.profiles_1d.triangularity_upper[end]

dyexp["equilibrium.time_slice[:].boundary.squareness_lower_inner"] =
    (; time_slice, _...) -> time_slice.profiles_1d.squareness_lower_inner[end]

dyexp["equilibrium.time_slice[:].boundary.squareness_upper_inner"] =
    (; time_slice, _...) -> time_slice.profiles_1d.squareness_upper_inner[end]

dyexp["equilibrium.time_slice[:].boundary.squareness_lower_outer"] =
    (; time_slice, _...) -> time_slice.profiles_1d.squareness_lower_outer[end]

dyexp["equilibrium.time_slice[:].boundary.squareness_upper_outer"] =
    (; time_slice, _...) -> time_slice.profiles_1d.squareness_upper_outer[end]

dyexp["equilibrium.time_slice[:].boundary.squareness"] =
    (; time_slice, _...) -> begin
        mxh = MXH(time_slice.boundary.outline.r, time_slice.boundary.outline.z, 2)
        return -mxh.s[2]
    end

dyexp["equilibrium.time_slice[:].profiles_1d.j_tor"] =
    (; profiles_1d, _...) -> J_tor(profiles_1d)

dyexp["equilibrium.time_slice[:].profiles_1d.j_parallel"] =
    (; time_slice, profiles_1d, _...) -> Jtor_2_Jpar(profiles_1d.rho_tor_norm, profiles_1d.j_tor, true, time_slice)


dyexp["equilibrium.time_slice[:].profiles_1d.dpressure_dpsi"] =
    (; time_slice, profiles_1d, _...) -> gradient(profiles_1d.psi, profiles_1d.pressure)

dyexp["equilibrium.time_slice[:].profiles_1d.f_df_dpsi"] =
    (; time_slice, profiles_1d, _...) -> profiles_1d.f .* gradient(profiles_1d.psi, profiles_1d.f)


dyexp["equilibrium.time_slice[:].profiles_1d.dvolume_drho_tor"] =
    (; profiles_1d, _...) -> gradient(profiles_1d.rho_tor, profiles_1d.volume)

dyexp["equilibrium.time_slice[:].profiles_1d.dvolume_dpsi"] =
    (; profiles_1d, _...) -> gradient(profiles_1d.psi, profiles_1d.volume)

dyexp["equilibrium.time_slice[:].profiles_1d.darea_drho_tor"] =
    (; profiles_1d, _...) -> gradient(profiles_1d.rho_tor, profiles_1d.area)

dyexp["equilibrium.time_slice[:].profiles_1d.darea_dpsi"] =
    (; profiles_1d, _...) -> gradient(profiles_1d.psi, profiles_1d.area)

dyexp["equilibrium.time_slice[:].profiles_1d.dpsi_drho_tor"] =
    (; profiles_1d, _...) -> gradient(profiles_1d.rho_tor, profiles_1d.psi)

dyexp["equilibrium.time_slice[:].profiles_1d.psi_norm"] =
    (; profiles_1d, _...) -> norm01(profiles_1d.psi)

# 2D
dyexp["equilibrium.time_slice[:].profiles_2d[:].r"] =
    (; profiles_2d, _...) -> ones(length(profiles_2d.grid.dim2))' .* profiles_2d.grid.dim1

dyexp["equilibrium.time_slice[:].profiles_2d[:].z"] =
    (; profiles_2d, _...) -> profiles_2d.grid.dim2' .* ones(length(profiles_2d.grid.dim1))

dyexp["equilibrium.time_slice[:].profiles_2d[:].b_field_tor"] =
    (; time_slice, profiles_2d, _...) -> begin
        psi = time_slice.profiles_1d.psi
        fpol = time_slice.profiles_1d.f
        FPOL = extrap1d(interp1d_itp(psi, fpol, :cubic); first=:constant, last=:constant).(profiles_2d.psi)
        return FPOL ./ profiles_2d.r
    end

dyexp["equilibrium.time_slice[:].profiles_2d[:].b_field_r"] =
    (; profiles_2d, _...) -> begin
        return Br_Bz(profiles_2d).Br
    end

dyexp["equilibrium.time_slice[:].profiles_2d[:].b_field_z"] =
    (; profiles_2d, _...) -> begin
        return Br_Bz(profiles_2d).Bz
    end

dyexp["equilibrium.time_slice[:].profiles_2d[:].j_tor"] =
    (; profiles_2d, _...) -> begin
        dim1, dim2 = profiles_2d.grid.dim1, profiles_2d.grid.dim2
        dBzdR = gradient(dim1, dim2, profiles_2d.b_field_z, 1)
        dBrdZ = gradient(dim1, dim2, profiles_2d.b_field_r, 2)
        return (dBrdZ - dBzdR) ./ mks.μ_0
    end

#= ============ =#
#  core_sources  #
#= ============ =#
dyexp["core_sources.source[:].profiles_1d[:].electrons.power_inside"] =
    (; profiles_1d, _...) -> begin
        energy = profiles_1d.electrons.energy
        cumtrapz(profiles_1d.grid.volume, energy)
    end

dyexp["core_sources.source[:].profiles_1d[:].electrons.energy"] =
    (; profiles_1d, _...) -> begin
        power_inside = profiles_1d.electrons.power_inside
        gradient(profiles_1d.grid.volume, power_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].total_ion_power_inside"] =
    (; profiles_1d, _...) -> begin
        total_ion_energy = profiles_1d.total_ion_energy
        cumtrapz(profiles_1d.grid.volume, total_ion_energy)
    end

dyexp["core_sources.source[:].profiles_1d[:].total_ion_energy"] =
    (; profiles_1d, _...) -> begin
        total_ion_power_inside = profiles_1d.total_ion_power_inside
        gradient(profiles_1d.grid.volume, total_ion_power_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].electrons.particles_inside"] =
    (; profiles_1d, _...) -> begin
        particles = profiles_1d.electrons.particles
        cumtrapz(profiles_1d.grid.volume, particles)
    end

dyexp["core_sources.source[:].profiles_1d[:].electrons.particles"] =
    (; profiles_1d, _...) -> begin
        particles_inside = profiles_1d.electrons.particles_inside
        gradient(profiles_1d.grid.volume, particles_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].current_parallel_inside"] =
    (; profiles_1d, _...) -> begin
        j_parallel = profiles_1d.j_parallel
        cumtrapz(profiles_1d.grid.area, j_parallel)
    end

dyexp["core_sources.source[:].profiles_1d[:].j_parallel"] =
    (; profiles_1d, _...) -> begin
        current_parallel_inside = profiles_1d.current_parallel_inside
        gradient(profiles_1d.grid.area, current_parallel_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].torque_tor_inside"] =
    (; profiles_1d, _...) -> begin
        momentum_tor = profiles_1d.momentum_tor
        cumtrapz(profiles_1d.grid.volume, momentum_tor)
    end

dyexp["core_sources.source[:].profiles_1d[:].momentum_tor"] =
    (; profiles_1d, _...) -> begin
        torque_tor_inside = profiles_1d.torque_tor_inside
        gradient(profiles_1d.grid.volume, torque_tor_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].ion[:].particles_inside"] =
    (; profiles_1d, ion, _...) -> begin
        particles = ion.particles
        cumtrapz(profiles_1d.grid.volume, particles)
    end

dyexp["core_sources.source[:].profiles_1d[:].ion[:].particles"] =
    (; profiles_1d, ion, _...) -> begin
        particles_inside = ion.particles_inside
        gradient(profiles_1d.grid.volume, ion.particles_inside)
    end

#= ===== =#
#  build  #
#= ===== =#
dyexp["build.layer[:].identifier"] =
    (; build, layer, _...) -> begin
        plasma_index = index(get_build_layer(build.layer; type=_plasma_))
        in_len = length(get_build_layers(build.layer; fs=_in_))
        return -(index(layer) - plasma_index) * layer.side
    end

dyexp["build.layer[:].outline.r"] =
    (; layer, _...) -> opposite_side_layer(layer).outline.r

dyexp["build.layer[:].outline.z"] =
    (; layer, _...) -> opposite_side_layer(layer).outline.z

dyexp["build.layer[:].shape"] =
    (; layer, _...) -> opposite_side_layer(layer).shape

dyexp["build.layer[:].shape_parameters"] =
    (; layer, _...) -> opposite_side_layer(layer).shape_parameters

dyexp["build.layer[:].start_radius"] =
    (; build, layer_index, _...) -> begin
        if build.layer[layer_index].side == Int(_out_)
            return 0.0
        else
            build_radii(build.layer)[1:end-1][layer_index]
        end
    end

dyexp["build.layer[:].end_radius"] =
    (; build, layer_index, _...) -> build_radii(build.layer)[2:end][layer_index]

dyexp["build.layer[:].area"] =
    (; layer, _...) -> area(layer)

dyexp["build.layer[:].volume"] =
    (; layer, _...) -> volume(layer)

dyexp["build.tf.ripple"] =
    (; build, _...) -> tf_ripple(get_build_layer(build.layer; type=_plasma_).end_radius, get_build_layer(build.layer; type=_tf_, fs=_lfs_).start_radius, build.tf.coils_n)

dyexp["build.tf.wedge_thickness"] =
    (; build, _...) -> 2π * get_build_layer(build.layer; type=_tf_, fs=_hfs_).end_radius / build.tf.coils_n

#= ======= =#
#  costing  #
#= ======= =#
dyexp["costing.cost_direct_capital.system[:].cost"] =
    (; system, _...) -> isempty(system.subsystem) ? error("no subsystem") : sum(sub.cost for sub in system.subsystem if !ismissing(sub, :cost))

dyexp["costing.cost_direct_capital.cost"] =
    (; cost_direct_capital, _...) -> isempty(cost_direct_capital.system) ? error("no system") : sum(sys.cost for sys in cost_direct_capital.system if !ismissing(sys, :cost))

dyexp["costing.cost_operations.system[:].yearly_cost"] =
    (; system, _...) -> isempty(system.subsystem) ? error("no subsystem") : sum(sub.yearly_cost for sub in system.subsystem if !ismissing(sub, :yearly_cost))

dyexp["costing.cost_operations.yearly_cost"] =
    (; cost_operations, _...) -> isempty(cost_operations.system) ? error("no system") : sum(sys.yearly_cost for sys in cost_operations.system if !ismissing(sys, :yearly_cost))

dyexp["costing.cost_decommissioning.system[:].cost"] =
    (; system, _...) -> isempty(system.subsystem) ? error("no subsystem") : sum(sub.cost for sub in system.subsystem if !ismissing(sub, :cost))

dyexp["costing.cost_decommissioning.cost"] =
    (; cost_decommissioning, _...) -> isempty(cost_decommissioning.system) ? error("no system") : sum(sys.cost for sys in cost_decommissioning.system if !ismissing(sys, :cost))

#= ============== =#
#  BalanceOfPlant  #
#= ============== =#
dyexp["balance_of_plant.Q_plant"] =
    (; balance_of_plant, _...) -> balance_of_plant.power_plant.power_electric_generated ./ balance_of_plant.power_electric_plant_operation.total_power

dyexp["balance_of_plant.power_electric_net"] =
    (; balance_of_plant, _...) -> balance_of_plant.power_plant.power_electric_generated .- balance_of_plant.power_electric_plant_operation.total_power

dyexp["balance_of_plant.power_electric_plant_operation.total_power"] =
    (; power_electric_plant_operation, _...) -> sum(sys.power for sys in power_electric_plant_operation.system)

#= ============== =#
#  pulse_schedule  #
#= ============== =#
dyexp["pulse_schedule.tf.b_field_tor_vacuum_r.reference"] =
    (; tf, _...) -> tf.r0 .* tf.b_field_tor_vacuum.reference

dyexp["pulse_schedule.tf.r0"] =
    (; dd, _...) -> dd.equilibrium.vacuum_toroidal_field.r0

dyexp["pulse_schedule.tf.b_field_tor_vacuum.reference"] =
    (; dd, _...) -> dd.equilibrium.vacuum_toroidal_field.b0

dyexp["pulse_schedule.tf.time"] =
    (; dd, _...) -> dd.equilibrium.time

dyexp["pulse_schedule.nbi.power.reference"] =
    (; nbi, _...) -> sum(unit.power.reference for unit in nbi.unit)

dyexp["pulse_schedule.ec.power_launched.reference"] =
    (; ec, _...) -> sum(beam.power_launched.reference for beam in ec.unit)

dyexp["pulse_schedule.ic.power.reference"] =
    (; ic, _...) -> sum(antenna.power.reference for antenna in ic.antenna)

dyexp["pulse_schedule.lh.power.reference"] =
    (; lh, _...) -> sum(antenna.power.reference for antenna in lh.antenna)

#= ====== =#
#  limits  #
#= ====== =#
dyexp["limits.model[:].cleared"] =
    (; model, _...) -> Int.(model.fraction .< 1.0)

dyexp["limits.all_cleared"] =
    (; limits, _...) -> begin
        time = limits.time
        all_cleared = ones(Int, length(time))
        for model in limits.model
            all_cleared .= all_cleared .* model.cleared
        end
        return all_cleared
    end

#= ============== =#
#  interferometer  #
#= ============== =#
dyexp["interferometer.channel[:].n_e_line_average.data"] = (; channel, _...) -> channel.n_e_line.data ./ lenght_line_of_sight(channel.line_of_sight)
dyexp["interferometer.channel[:].n_e_line_average.time"] = (; channel, _...) -> channel.n_e_line.time

#= ======= =#
#  summary  #
#= ======= =#
dyexp["summary.fusion.power.value"] = # NOTE: This is the fusion power that is coupled to the plasma
    (; dd, summary, _...) -> [fusion_plasma_power(dd.core_profiles.profiles_1d[time0]) for time0 in summary.time]

dyexp["summary.global_quantities.ip.value"] =
    (; dd, summary, _...) -> [dd.equilibrium.time_slice[time0].global_quantities.ip for time0 in summary.time]

dyexp["summary.global_quantities.current_bootstrap.value"] =
    (; dd, summary, _...) -> begin
        tmp = eltype(summary)[]
        for time0 in summary.time
            cp1d = dd.core_profiles.profiles_1d[time0]
            eqt = dd.equilibrium.time_slice[time0]
            push!(tmp, Ip_bootstrap(cp1d, eqt))
        end
        return tmp
    end

dyexp["summary.global_quantities.current_non_inductive.value"] =
    (; dd, summary, _...) -> begin
        tmp = eltype(summary)[]
        for time0 in summary.time
            cp1d = dd.core_profiles.profiles_1d[time0]
            eqt = dd.equilibrium.time_slice[time0]
            push!(tmp, Ip_non_inductive(cp1d, eqt))
        end
        return tmp
    end

dyexp["summary.global_quantities.current_ohm.value"] =
    (; dd, summary, _...) -> begin
        tmp = eltype(summary)[]
        for time0 in summary.time
            cp1d = dd.core_profiles.profiles_1d[time0]
            eqt = dd.equilibrium.time_slice[time0]
            push!(tmp, Ip_ohmic(cp1d, eqt))
        end
        return tmp
    end


dyexp["summary.global_quantities.beta_pol_mhd.value"] =
    (; dd, summary, _...) -> [dd.equilibrium.time_slice[time0].global_quantities.beta_pol for time0 in summary.time]

dyexp["summary.global_quantities.beta_tor.value"] =
    (; dd, summary, _...) -> [beta_tor(dd.equilibrium.time_slice[time0], dd.core_profiles.profiles_1d[time0]) for time0 in summary.time]

dyexp["summary.global_quantities.beta_tor_mhd.value"] =
    (; dd, summary, _...) -> [dd.equilibrium.time_slice[time0].global_quantities.beta_tor for time0 in summary.time]

dyexp["summary.global_quantities.beta_tor_norm_mhd.value"] =
    (; dd, summary, _...) -> [dd.equilibrium.time_slice[time0].global_quantities.beta_normal for time0 in summary.time]

dyexp["summary.global_quantities.beta_tor_norm.value"] =
    (; dd, summary, _...) -> [beta_tor_norm(dd.equilibrium.time_slice[time0], dd.core_profiles.profiles_1d[time0]) for time0 in summary.time]

dyexp["summary.global_quantities.beta_tor_thermal_norm.value"] =
    (; dd, summary, _...) -> [beta_tor_thermal_norm(dd.equilibrium.time_slice[time0], dd.core_profiles.profiles_1d[time0]) for time0 in summary.time]

dyexp["summary.global_quantities.energy_thermal.value"] =
    (; dd, summary, _...) -> [energy_thermal(dd.core_profiles.profiles_1d[time0]) for time0 in summary.time]

dyexp["summary.global_quantities.tau_energy.value"] =
    (; dd, summary, _...) -> [tau_e_thermal(dd; time0, include_radiation=true) for time0 in summary.time]

dyexp["summary.global_quantities.tau_energy_98.value"] =
    (; dd, summary, _...) -> [tau_e_h98(dd; time0, include_radiation=true) for time0 in summary.time]

dyexp["summary.global_quantities.h_98.value"] =
    (; dd, summary, _...) -> summary.global_quantities.tau_energy.value ./ summary.global_quantities.tau_energy_98.value


dyexp["summary.heating_current_drive.power_launched_ec.value"] =
    (; dd, summary, _...) -> sum(interp1d(beam.power_launched.time, beam.power_launched.data, :constant).(summary.time) for beam in dd.ec_launchers.beam)

dyexp["summary.heating_current_drive.power_launched_ic.value"] =
    (; dd, summary, _...) -> sum(interp1d(antenna.power_launched.time, antenna.power_launched.data, :constant).(summary.time) for antenna in dd.ic_antennas.antenna)

dyexp["summary.heating_current_drive.power_launched_lh.value"] =
    (; dd, summary, _...) -> sum(interp1d(antenna.power_launched.time, antenna.power_launched.data, :constant).(summary.time) for antenna in dd.lh_antennas.antenna)

dyexp["summary.heating_current_drive.power_launched_nbi.value"] =
    (; dd, summary, _...) -> sum(interp1d(unit.power_launched.time, unit.power_launched.data, :constant).(summary.time) for unit in dd.nbi.unit)

dyexp["summary.heating_current_drive.power_launched_total.value"] =
    (; dd, summary, _...) ->
        getproperty(dd.summary.heating_current_drive.power_launched_nbi, :value, zeros(length(summary.time))) .+
        getproperty(dd.summary.heating_current_drive.power_launched_ec, :value, zeros(length(summary.time))) .+
        getproperty(dd.summary.heating_current_drive.power_launched_ic, :value, zeros(length(summary.time))) .+
        getproperty(dd.summary.heating_current_drive.power_launched_lh, :value, zeros(length(summary.time)))


dyexp["summary.local.magnetic_axis.t_e.value"] =
    (; dd, summary, _...) -> [dd.core_profiles.profiles_1d[time0].electrons.temperature[1] for time0 in summary.time]

dyexp["summary.local.magnetic_axis.n_e.value"] =
    (; dd, summary, _...) -> [dd.core_profiles.profiles_1d[time0].electrons.density[1] for time0 in summary.time]

dyexp["summary.local.magnetic_axis.t_i_average.value"] =
    (; dd, summary, _...) -> [dd.core_profiles.profiles_1d[time0].t_i_average[1] for time0 in summary.time]

dyexp["summary.local.magnetic_axis.zeff.value"] =
    (; dd, summary, _...) -> [dd.core_profiles.profiles_1d[time0].zeff[1] for time0 in summary.time]


dyexp["summary.local.separatrix.t_e.value"] =
    (; dd, summary, _...) -> [dd.core_profiles.profiles_1d[time0].electrons.temperature[end] for time0 in summary.time]

dyexp["summary.local.separatrix.n_e.value"] =
    (; dd, summary, _...) -> [dd.core_profiles.profiles_1d[time0].electrons.density[end] for time0 in summary.time]

dyexp["summary.local.separatrix.t_i_average.value"] =
    (; dd, summary, _...) -> [dd.core_profiles.profiles_1d[time0].t_i_average[end] for time0 in summary.time]

dyexp["summary.local.separatrix.zeff.value"] =
    (; dd, summary, _...) -> [dd.core_profiles.profiles_1d[time0].zeff[end] for time0 in summary.time]


dyexp["summary.volume_average.t_e.value"] =
    (; dd, summary, _...) -> begin
        tmp = eltype(summary)[]
        for time0 in summary.time
            cp1d = dd.core_profiles.profiles_1d[time0]
            push!(tmp, trapz(cp1d.grid.volume, cp1d.electrons.temperature) / cp1d.grid.volume[end])
        end
        return tmp
    end

dyexp["summary.volume_average.n_e.value"] =
    (; dd, summary, _...) -> begin
        tmp = eltype(summary)[]
        for time0 in summary.time
            cp1d = dd.core_profiles.profiles_1d[time0]
            push!(tmp, trapz(cp1d.grid.volume, cp1d.electrons.density) / cp1d.grid.volume[end])
        end
        return tmp
    end

dyexp["summary.volume_average.t_i_average.value"] =
    (; dd, summary, _...) -> begin
        tmp = eltype(summary)[]
        for time0 in summary.time
            cp1d = dd.core_profiles.profiles_1d[time0]
            push!(tmp, trapz(cp1d.grid.volume, cp1d.t_i_average) / cp1d.grid.volume[end])
        end
        return tmp
    end

dyexp["summary.volume_average.zeff.value"] =
    (; dd, summary, _...) -> begin
        tmp = eltype(summary)[]
        for time0 in summary.time
            cp1d = dd.core_profiles.profiles_1d[time0]
            push!(tmp, trapz(cp1d.grid.volume, cp1d.zeff) / cp1d.grid.volume[end])
        end
        return tmp
    end

# ============ #

Base.Docs.@doc """
    dynamic_expressions = Dict{String,Function}()

Expressions
* `$(join(sort!(collect(keys(dynamic_expressions))),"`\n* `"))`
""" dynamic_expressions

push!(document[:Expressions], :dynamic_expressions)
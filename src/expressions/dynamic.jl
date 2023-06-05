import NumericalIntegration: integrate, cumul_integrate

function IMASDD.get_expressions(::Val{:dynamic})
    return dynamic_expressions
end

const dynamic_expressions = dyexp = Dict{String,Function}()

# NOTE: make sure that expressions accept as argument (not keyword argument)
# the coordinates of the quantitiy you are writing the expression of
# 
# For example, this will FAIL:
#    dyexp["core_profiles.profiles_1d[:].electrons.pressure"] =
#         (; electrons, _...) -> electrons.temperature .* electrons.density * 1.60218e-19
#
# This is GOOD:
#    dyexp["core_profiles.profiles_1d[:].electrons.pressure"] =
#         (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density * 1.60218e-19

#= =========== =#
# core_profiles #
#= =========== =#
dyexp["core_profiles.profiles_1d[:].electrons.density"] =
    (rho_tor_norm; electrons, _...) -> electrons.density_thermal .+ electrons.density_fast

dyexp["core_profiles.profiles_1d[:].electrons.density_fast"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

dyexp["core_profiles.profiles_1d[:].electrons.pressure_thermal"] =
    (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density .* constants.e

dyexp["core_profiles.profiles_1d[:].electrons.pressure_fast_parallel"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

dyexp["core_profiles.profiles_1d[:].electrons.pressure_fast_perpendicular"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

dyexp["core_profiles.profiles_1d[:].electrons.pressure"] =
    (rho_tor_norm; electrons, _...) -> electrons.pressure_thermal .+ electrons.pressure_fast_parallel .+ electrons.pressure_fast_perpendicular .* 2.0


dyexp["core_profiles.profiles_1d[:].ion[:].z_ion"] =
    (; ion, _...) -> begin
        if length(ion.element) == 1
            return ion.element[1].z_n
        else
            error("z_ion expression does not yet handle multiple charge states")
        end
    end

dyexp["core_profiles.profiles_1d[:].t_i_average"] =
    (rho_tor_norm; profiles_1d, _...) -> t_i_average(profiles_1d)

dyexp["core_profiles.profiles_1d[:].ion[:].density"] =
    (rho_tor_norm; ion, _...) -> ion.density_thermal .+ ion.density_fast

dyexp["core_profiles.profiles_1d[:].ion[:].density_fast"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

dyexp["core_profiles.profiles_1d[:].ion[:].pressure_thermal"] =
    (rho_tor_norm; ion, _...) -> ion.temperature .* ion.density .* constants.e

dyexp["core_profiles.profiles_1d[:].ion[:].pressure_fast_parallel"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

dyexp["core_profiles.profiles_1d[:].ion[:].pressure_fast_perpendicular"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

dyexp["core_profiles.profiles_1d[:]..ion[:].pressure"] =
    (rho_tor_norm; ion, _...) -> ion.pressure_thermal .+ ion.pressure_fast_parallel .+ ion.pressure_fast_perpendicular .* 2.0


dyexp["core_profiles.profiles_1d[:].pressure_ion_total"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        tmp = zero(rho_tor_norm)
        for ion in profiles_1d.ion
            tmp .+= ion.pressure_thermal
        end
        return tmp
    end

dyexp["core_profiles.profiles_1d[:].pressure_thermal"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.electrons.pressure_thermal .+ profiles_1d.pressure_ion_total

dyexp["core_profiles.profiles_1d[:].pressure_parallel"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.pressure_thermal ./ 3.0 .+ profiles_1d.electrons.pressure_fast_parallel .+ sum(ion.pressure_fast_parallel for ion in profiles_1d.ion)

dyexp["core_profiles.profiles_1d[:].pressure_perpendicular"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.pressure_thermal ./ 3.0 .+ profiles_1d.electrons.pressure_fast_perpendicular .+ sum(ion.pressure_fast_perpendicular for ion in profiles_1d.ion)

dyexp["core_profiles.profiles_1d[:].pressure"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.pressure_perpendicular .* 2.0 .+ profiles_1d.pressure_parallel


dyexp["core_profiles.profiles_1d[:].conductivity_parallel"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> nclass_conductivity(dd.equilibrium.time_slice[Float64(profiles_1d.time)], profiles_1d)

dyexp["core_profiles.profiles_1d[:].j_bootstrap"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> Sauter_neo2021_bootstrap(dd.equilibrium.time_slice[Float64(profiles_1d.time)], profiles_1d)

dyexp["core_profiles.profiles_1d[:].j_ohmic"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.j_total .- profiles_1d.j_non_inductive

dyexp["core_profiles.profiles_1d[:].j_non_inductive"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> total_sources(dd.core_sources, profiles_1d; exclude_indexes=[7, 13]).j_parallel .+ profiles_1d.j_bootstrap

dyexp["core_profiles.profiles_1d[:].j_total"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> profiles_1d.j_ohmic .+ profiles_1d.j_non_inductive

dyexp["core_profiles.profiles_1d[:].j_tor"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        Jpar_2_Jtor(rho_tor_norm, profiles_1d.j_total, true, eqt)
    end

dyexp["core_profiles.profiles_1d[:].time"] =
    (; core_profiles, profiles_1d_index, _...) -> core_profiles.time[profiles_1d_index]

dyexp["core_profiles.vacuum_toroidal_field.b0"] =
    (time; dd, core_profiles, _...) -> interp1d(dd.equilibrium.time, dd.equilibrium.vacuum_toroidal_field.b0, :constant).(core_profiles.time)

dyexp["core_profiles.vacuum_toroidal_field.r0"] =
    (; dd, _...) -> dd.equilibrium.vacuum_toroidal_field.r0

#= ============ =#
# core_transport #
#= ============ =#
dyexp["core_profiles.profiles_1d[:].time"] =
    (; core_profiles, profiles_1d_index, _...) -> core_profiles.time[profiles_1d_index]

dyexp["core_profiles.vacuum_toroidal_field.b0"] =
    (time; dd, core_profiles, _...) -> interp1d(dd.equilibrium.time, dd.equilibrium.vacuum_toroidal_field.b0, :constant).(core_profiles.time)

dyexp["core_profiles.vacuum_toroidal_field.r0"] =
    (; dd, _...) -> dd.equilibrium.vacuum_toroidal_field.r0

#= ========= =#
# equilibrium #
#= ========= =#
# IMAS does not hold B0 information in a given time slice, but we can get that info from `B0=f/R0`
# This trick propagates the B0 information to a time_slice even when that time_slice has not been initialized with profiles_1d data
dyexp["equilibrium.time_slice[:].profiles_1d.f"] =
    (psi; equilibrium, time_slice_index, _...) -> (psi === missing ? [1] : ones(size(psi))) .* (equilibrium.vacuum_toroidal_field.b0[time_slice_index] * equilibrium.vacuum_toroidal_field.r0)

dyexp["equilibrium.time_slice[:].global_quantities.energy_mhd"] =
    (; time_slice, _...) -> 3 / 2 * integrate(time_slice.profiles_1d.volume, time_slice.profiles_1d.pressure)

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

dyexp["equilibrium.time_slice[:].global_quantities.magnetic_axis.b_field_tor"] =
    (; equilibrium, time_slice_index, _...) -> equilibrium.vacuum_toroidal_field.b0[time_slice_index] * equilibrium.vacuum_toroidal_field.r0 / equilibrium.time_slice[time_slice_index].boundary.geometric_axis.r

dyexp["equilibrium.time_slice[:].profiles_1d.geometric_axis.r"] =
    (psi; time_slice, _...) -> (time_slice.profiles_1d.r_outboard .+ time_slice.profiles_1d.r_inboard) .* 0.5

dyexp["equilibrium.time_slice[:].profiles_1d.geometric_axis.z"] =
    (psi; time_slice, _...) -> psi .* 0.0 .+ time_slice.global_quantities.magnetic_axis.z

dyexp["equilibrium.time_slice[:].boundary.geometric_axis.r"] =
    (; time_slice, _...) -> time_slice.profiles_1d.geometric_axis.r[end]

dyexp["equilibrium.time_slice[:].boundary.geometric_axis.z"] =
    (; time_slice, _...) -> time_slice.profiles_1d.geometric_axis.z[end]

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
    (psi; profiles_1d, _...) -> begin
        j_parallel = profiles_1d.j_parallel
        Jpar_2_Jtor(profiles_1d.rho_tor_norm, j_parallel, true, time_slice)
    end

dyexp["equilibrium.time_slice[:].profiles_1d.j_parallel"] =
    (psi; time_slice, profiles_1d, _...) -> begin
        j_tor = profiles_1d.j_tor
        Jtor_2_Jpar(profiles_1d.rho_tor_norm, j_tor, true, time_slice)
    end

dyexp["equilibrium.time_slice[:].profiles_1d.dpressure_dpsi"] =
    (psi; time_slice, profiles_1d, _...) -> gradient(psi, profiles_1d.pressure)

dyexp["equilibrium.time_slice[:].time"] =
    (; equilibrium, time_slice_index, _...) -> equilibrium.time[time_slice_index]

dyexp["equilibrium.time_slice[:].profiles_1d.psi_norm"] =
    (psi; _...) -> norm01(psi)


dyexp["equilibrium.time_slice[:].profiles_2d[:].r"] =
    (dim1, dim2; _...) -> ones(length(dim2))' .* dim1

dyexp["equilibrium.time_slice[:].profiles_2d[:].z"] =
    (dim1, dim2; _...) -> dim2' .* ones(length(dim1))

dyexp["equilibrium.time_slice[:].profiles_2d[:].b_field_tor"] =
    (dim1, dim2; time_slice, profiles_2d, _...) -> begin
        psi = time_slice.profiles_1d.psi
        fpol = time_slice.profiles_1d.f
        FPOL = extrap1d(interp1d(psi, fpol, :cubic); first=:flat, last=:flat).(profiles_2d.psi)
        return FPOL ./ profiles_2d.r
    end

dyexp["equilibrium.time_slice[:].profiles_2d[:].b_field_r"] =
    (dim1, dim2; profiles_2d, _...) -> begin
        return Br_Bz(profiles_2d)[1]
    end

dyexp["equilibrium.time_slice[:].profiles_2d[:].b_field_z"] =
    (dim1, dim2; profiles_2d, _...) -> begin
        return Br_Bz(profiles_2d)[2]
    end

dyexp["equilibrium.time_slice[:].profiles_2d[:].j_tor"] =
    (dim1, dim2; profiles_2d, _...) -> begin
        dBrdZ = gradient(dim1, dim2, profiles_2d.b_field_r; dim=2)
        dBzdR = gradient(dim1, dim2, profiles_2d.b_field_z, dim=1)
        return (dBrdZ - dBzdR) ./ constants.μ_0
    end

#= ============ =#
#  core_sources  #
#= ============ =#
dyexp["core_sources.source[:].profiles_1d[:].electrons.power_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        energy = profiles_1d.electrons.energy
        cumul_integrate(profiles_1d.grid.volume, energy)
    end

dyexp["core_sources.source[:].profiles_1d[:].electrons.energy"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        power_inside = profiles_1d.electrons.power_inside
        gradient(profiles_1d.grid.volume, power_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].total_ion_power_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        total_ion_energy = profiles_1d.total_ion_energy
        cumul_integrate(profiles_1d.grid.volume, total_ion_energy)
    end

dyexp["core_sources.source[:].profiles_1d[:].total_ion_energy"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        total_ion_power_inside = profiles_1d.total_ion_power_inside
        gradient(profiles_1d.grid.volume, total_ion_power_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].electrons.particles_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        particles = profiles_1d.electrons.particles
        cumul_integrate(profiles_1d.grid.volume, particles)
    end

dyexp["core_sources.source[:].profiles_1d[:].electrons.particles"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        particles_inside = profiles_1d.electrons.particles_inside
        gradient(profiles_1d.grid.volume, particles_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].current_parallel_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        j_parallel = profiles_1d.j_parallel
        cumul_integrate(profiles_1d.grid.area, j_parallel)
    end

dyexp["core_sources.source[:].profiles_1d[:].j_parallel"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        current_parallel_inside = profiles_1d.current_parallel_inside
        gradient(profiles_1d.grid.area, current_parallel_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].torque_tor_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        momentum_tor = profiles_1d.momentum_tor
        cumul_integrate(profiles_1d.grid.volume, momentum_tor)
    end

dyexp["core_sources.source[:].profiles_1d[:].momentum_tor"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        torque_tor_inside = profiles_1d.torque_tor_inside
        gradient(profiles_1d.grid.volume, torque_tor_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].ion[:].particles_inside"] =
    (rho_tor_norm; profiles_1d, ion, _...) -> begin
        particles = ion.particles
        cumul_integrate(profiles_1d.grid.volume, particles)
    end

dyexp["core_sources.source[:].profiles_1d[:].ion[:].particles"] =
    (rho_tor_norm; profiles_1d, ion, _...) -> begin
        particles_inside = ion.particles_inside
        gradient(profiles_1d.grid.volume, ion.particles_inside)
    end


dyexp["core_sources.source[:].profiles_1d[:].time"] =
    (; core_sources, profiles_1d_index, _...) -> core_sources.time[profiles_1d_index]

dyexp["core_sources.source[:].global_quantities[:].time"] =
    (; core_sources, global_quantities_index, _...) -> core_sources.time[global_quantities_index]


dyexp["core_sources.vacuum_toroidal_field.b0"] =
    (time; dd, core_sources, _...) -> interp1d(dd.equilibrium.time, dd.equilibrium.vacuum_toroidal_field.b0, :constant).(core_sources.time)

dyexp["core_sources.vacuum_toroidal_field.r0"] =
    (; dd, _...) -> dd.equilibrium.vacuum_toroidal_field.r0

#= ===== =#
#  build  #
#= ===== =#
dyexp["build.layer[:].outline.r"] =
    (x; build, layer, _...) -> get_build_layer(build.layer; identifier=layer.identifier, fs=(layer.fs == Int(_lfs_)) ? _hfs_ : _lfs_).outline.r

dyexp["build.layer[:].outline.z"] =
    (x; build, layer, _...) -> get_build_layer(build.layer; identifier=layer.identifier, fs=(layer.fs == Int(_lfs_)) ? _hfs_ : _lfs_).outline.z

dyexp["build.layer[:].shape"] =
    (; build, layer, _...) -> get_build_layer(build.layer; identifier=layer.identifier, fs=(layer.fs == Int(_lfs_)) ? _hfs_ : _lfs_).shape

dyexp["build.layer[:].shape_parameters"] =
    (; build, layer, _...) -> get_build_layer(build.layer; identifier=layer.identifier, fs=(layer.fs == Int(_lfs_)) ? _hfs_ : _lfs_).shape_parameters

dyexp["build.layer[:].start_radius"] =
    (; build, layer_index, _...) -> build_radii(build)[1:end-1][layer_index]

dyexp["build.layer[:].end_radius"] =
    (; build, layer_index, _...) -> build_radii(build)[2:end][layer_index]

dyexp["build.layer[:].area"] =
    (; layer, _...) -> area(layer)

dyexp["build.layer[:].volume"] =
    (; layer, _...) -> volume(layer)

dyexp["build.tf.ripple"] =
    (; build, _...) -> tf_ripple(get_build_layer(build.layer, type=_plasma_).end_radius, get_build_layer(build.layer, type=_tf_, fs=_lfs_).start_radius, build.tf.coils_n)

dyexp["build.tf.wedge_thickness"] =
    (; build, _...) -> 2π * get_build_layer(build.layer, type=_tf_, fs=_hfs_).end_radius / build.tf.coils_n

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
    (; cost_operations, _...) -> isempty(cost_operations.system) ?  error("no system") : sum(sys.yearly_cost for sys in cost_operations.system if !ismissing(sys, :yearly_cost))

dyexp["costing.cost_decommissioning.system[:].cost"] =
    (; system, _...) -> isempty(system.subsystem) ? error("no subsystem") : sum(sub.cost for sub in system.subsystem if !ismissing(sub, :cost))

dyexp["costing.cost_decommissioning.cost"] =
    (; cost_decommissioning, _...) -> isempty(cost_decommissioning.system) ? error("no system") : sum(sys.cost for sys in cost_decommissioning.system if !ismissing(sys, :cost))

#= ============== =#
#  BalanceOfPlant  #
#= ============== =#
dyexp["balance_of_plant.Q_plant"] =
    (time; balance_of_plant, _...) -> balance_of_plant.thermal_cycle.power_electric_generated ./ balance_of_plant.power_electric_plant_operation.total_power

dyexp["balance_of_plant.power_electric_net"] =
    (time; balance_of_plant, _...) -> balance_of_plant.thermal_cycle.power_electric_generated .- balance_of_plant.power_electric_plant_operation.total_power

dyexp["balance_of_plant.power_electric_plant_operation.total_power"] =
    (time; power_electric_plant_operation, _...) -> sum(sys.power for sys in power_electric_plant_operation.system)

dyexp["balance_of_plant.thermal_cycle.total_useful_heat_power"] =
    (time; balance_of_plant, _...) -> balance_of_plant.heat_transfer.wall.heat_delivered .+ balance_of_plant.heat_transfer.divertor.heat_delivered .+ balance_of_plant.heat_transfer.breeder.heat_delivered

dyexp["balance_of_plant.thermal_cycle.power_electric_generated"] =
    (time; balance_of_plant, thermal_cycle, _...) -> thermal_cycle.net_work .* thermal_cycle.generator_conversion_efficiency

#= ========= =#
#  divertors  #
#= ========= =#
dyexp["divertors.divertor[:].power_black_body.time"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_black_body)[1]

dyexp["divertors.divertor[:].power_black_body.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_black_body)[2]


dyexp["divertors.divertor[:].power_conducted.time"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_conducted)[1]

dyexp["divertors.divertor[:].power_conducted.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_conducted)[2]


dyexp["divertors.divertor[:].power_convected.time"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_convected)[1]

dyexp["divertors.divertor[:].power_convected.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_convected)[2]


dyexp["divertors.divertor[:].power_currents.time"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_currents)[1]

dyexp["divertors.divertor[:].power_currents.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_currents)[2]


dyexp["divertors.divertor[:].power_incident.time"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_incident)[1]

dyexp["divertors.divertor[:].power_incident.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_incident)[2]


dyexp["divertors.divertor[:].power_neutrals.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_neutrals)[1]

dyexp["divertors.divertor[:].power_neutrals.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_neutrals)[2]


dyexp["divertors.divertor[:].power_radiated.time"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_radiated)[1]

dyexp["divertors.divertor[:].power_radiated.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_radiated)[2]


dyexp["divertors.divertor[:].power_recombination_neutrals.time"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_recombination_neutrals)[1]

dyexp["divertors.divertor[:].power_recombination_neutrals.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_recombination_neutrals)[2]


dyexp["divertors.divertor[:].power_recombination_plasma.time"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_recombination_plasma)[1]

dyexp["divertors.divertor[:].power_recombination_plasma.data"] =
    (time; divertor, _...) -> divertor_totals_from_targets(divertor, :power_recombination_plasma)[2]


#= ========= =#
#  stability  #
#= ========= =#
dyexp["stability.model[:].cleared"] =
    (time; model, _...) -> Int.(model.fraction .<= 1.0)

dyexp["stability.all_cleared"] =
    (time; stability, _...) -> begin
        all_cleared = ones(Int, length(time))
        for model in stability.model
            all_cleared .= all_cleared .* model.cleared
        end
        return all_cleared
    end

#= ======= =#
#  summary  #
#= ======= =#
dyexp["summary.fusion.power.value"] = # NOTE: This is the fusion power that is coupled to the plasma
    (time; dd, summary, _...) -> begin
        type = typeof(summary).parameters[1]
        tmp = type[]
        for time in summary.time
            push!(tmp, fusion_plasma_power(dd.core_profiles.profiles_1d[Float64(time)]))
        end
        return tmp
    end

dyexp["summary.global_quantities.ip.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.ip for time in summary.time]

dyexp["summary.global_quantities.b0.value"] =
    (time; dd, summary, _...) -> interp1d(dd.equilibrium.time, dd.equilibrium.vacuum_toroidal_field.b0, :constant).(summary.time)

dyexp["summary.global_quantities.beta_tor_norm.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.beta_normal for time in summary.time]

dyexp["summary.global_quantities.r0.value"] =
    (; dd, summary, _...) -> dd.equilibrium.vacuum_toroidal_field.r0

dyexp["summary.global_quantities.current_bootstrap.value"] =
    (time; dd, summary, _...) -> begin
        type = typeof(summary).parameters[1]
        tmp = type[]
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.area, cp1d.j_bootstrap))
        end
        return tmp
    end

dyexp["summary.global_quantities.current_non_inductive.value"] =
    (time; dd, summary, _...) -> begin
        type = typeof(summary).parameters[1]
        tmp = type[]
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.area, cp1d.j_non_inductive))
        end
        return tmp
    end

dyexp["summary.global_quantities.current_ohm.value"] =
    (time; dd, summary, _...) -> begin
        type = typeof(summary).parameters[1]
        tmp = type[]
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.area, cp1d.j_ohmic))
        end
        return tmp
    end


dyexp["summary.global_quantities.beta_pol_mhd.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.beta_pol for time in summary.time]

dyexp["summary.global_quantities.beta_tor.value"] =
    (time; dd, summary, _...) -> [beta_tor(dd.equilibrium, dd.core_profiles.profiles_1d[Float64(time)]) for time in summary.time]

dyexp["summary.global_quantities.beta_tor_mhd.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.beta_tor for time in summary.time]

dyexp["summary.global_quantities.beta_tor_norm_mhd.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.beta_normal for time in summary.time]

dyexp["summary.global_quantities.beta_tor_thermal_norm.value"] =
    (time; dd, summary, _...) -> [beta_tor_thermal_norm(dd.equilibrium, dd.core_profiles.profiles_1d[Float64(time)]) for time in summary.time]

dyexp["summary.global_quantities.energy_thermal.value"] =
    (time; dd, summary, _...) -> [energy_thermal(dd.core_profiles.profiles_1d[Float64(time)]) for time in summary.time]

dyexp["summary.global_quantities.tau_energy.value"] =
    (time; dd, summary, _...) -> [tau_e_thermal(dd.core_profiles.profiles_1d[Float64(time)], dd.core_sources) for time in summary.time]

dyexp["summary.global_quantities.tau_energy_98.value"] =
    (time; dd, summary, _...) -> [tau_e_h98(dd, time=time) for time in summary.time]

dyexp["summary.global_quantities.h_98.value"] =
    (time; dd, summary, _...) -> summary.global_quantities.tau_energy.value ./ summary.global_quantities.tau_energy_98.value


dyexp["summary.heating_current_drive.power_launched_ec.value"] =
    (time; dd, summary, _...) -> sum(interp1d(beam.power_launched.time, beam.power_launched.data, :constant).(summary.time) for beam in dd.ec_launchers.beam)

dyexp["summary.heating_current_drive.power_launched_ic.value"] =
    (time; dd, summary, _...) -> sum(interp1d(antenna.power_launched.time, antenna.power_launched.data, :constant).(summary.time) for antenna in dd.ic_antennas.antenna)

dyexp["summary.heating_current_drive.power_launched_lh.value"] =
    (time; dd, summary, _...) -> sum(interp1d(antenna.power_launched.time, antenna.power_launched.data, :constant).(summary.time) for antenna in dd.lh_antennas.antenna)

dyexp["summary.heating_current_drive.power_launched_nbi.value"] =
    (time; dd, summary, _...) -> sum(interp1d(unit.power_launched.time, unit.power_launched.data, :constant).(summary.time) for unit in dd.nbi.unit)

dyexp["summary.heating_current_drive.power_launched_total.value"] =
    (time; dd, summary, _...) -> getproperty(dd.summary.heating_current_drive.power_launched_nbi, :value, zeros(length(summary.time))) .+ getproperty(dd.summary.heating_current_drive.power_launched_ec, :value, zeros(length(summary.time))) .+ getproperty(dd.summary.heating_current_drive.power_launched_ic, :value, zeros(length(summary.time))) .+ getproperty(dd.summary.heating_current_drive.power_launched_lh, :value, zeros(length(summary.time)))


dyexp["summary.local.magnetic_axis.t_e.value"] =
    (time; dd, summary, _...) -> [dd.core_profiles.profiles_1d[Float64(time)].electrons.temperature[1] for time in summary.time]

dyexp["summary.local.magnetic_axis.n_e.value"] =
    (time; dd, summary, _...) -> [dd.core_profiles.profiles_1d[Float64(time)].electrons.density[1] for time in summary.time]

dyexp["summary.local.magnetic_axis.t_i_average.value"] =
    (time; dd, summary, _...) -> [dd.core_profiles.profiles_1d[Float64(time)].t_i_average[1] for time in summary.time]

dyexp["summary.local.magnetic_axis.zeff.value"] =
    (time; dd, summary, _...) -> [dd.core_profiles.profiles_1d[Float64(time)].zeff[1] for time in summary.time]


dyexp["summary.local.separatrix.t_e.value"] =
    (time; dd, summary, _...) -> [dd.core_profiles.profiles_1d[Float64(time)].electrons.temperature[end] for time in summary.time]

dyexp["summary.local.separatrix.n_e.value"] =
    (time; dd, summary, _...) -> [dd.core_profiles.profiles_1d[Float64(time)].electrons.density[end] for time in summary.time]

dyexp["summary.local.separatrix.t_i_average.value"] =
    (time; dd, summary, _...) -> [dd.core_profiles.profiles_1d[Float64(time)].t_i_average[end] for time in summary.time]

dyexp["summary.local.separatrix.zeff.value"] =
    (time; dd, summary, _...) -> [dd.core_profiles.profiles_1d[Float64(time)].zeff[end] for time in summary.time]


dyexp["summary.volume_average.t_e.value"] =
    (time; dd, summary, _...) -> begin
        type = typeof(summary).parameters[1]
        tmp = type[]
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.volume, cp1d.electrons.temperature) / cp1d.grid.volume[end])
        end
        return tmp
    end

dyexp["summary.volume_average.n_e.value"] =
    (time; dd, summary, _...) -> begin
        type = typeof(summary).parameters[1]
        tmp = type[]
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.volume, cp1d.electrons.density) / cp1d.grid.volume[end])
        end
        return tmp
    end

dyexp["summary.volume_average.t_i_average.value"] =
    (time; dd, summary, _...) -> begin
        type = typeof(summary).parameters[1]
        tmp = type[]
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.volume, cp1d.t_i_average) / cp1d.grid.volume[end])
        end
        return tmp
    end

dyexp["summary.volume_average.zeff.value"] =
    (time; dd, summary, _...) -> begin
        type = typeof(summary).parameters[1]
        tmp = type[]
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.volume, cp1d.zeff) / cp1d.grid.volume[end])
        end
        return tmp
    end

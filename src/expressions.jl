import NumericalIntegration: integrate, cumul_integrate

function IMASDD.get_expressions(::Val{:IMAS})
    return expressions
end

const expressions = Dict{String,Function}()

# NOTE: make sure that expressions accept as argument (not keyword argument)
# the coordinates of the quantitiy you are writing the expression of
# 
# For example, this will FAIL:
#    expressions["core_profiles.profiles_1d[:].electrons.pressure"] =
#               (; electrons, _...) -> electrons.temperature .* electrons.density * 1.60218e-19
#
# This is GOOD:
#    expressions["core_profiles.profiles_1d[:].electrons.pressure"] =
#               (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density * 1.60218e-19

#= =========== =#
# core_profiles #
#= =========== =#
expressions["core_profiles.profiles_1d[:].electrons.density"] =
    (rho_tor_norm; electrons, _...) -> electrons.density_thermal .+ electrons.density_fast

expressions["core_profiles.profiles_1d[:].electrons.density_fast"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

expressions["core_profiles.profiles_1d[:].electrons.pressure_thermal"] =
    (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density .* constants.e

expressions["core_profiles.profiles_1d[:].electrons.pressure_fast_parallel"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

expressions["core_profiles.profiles_1d[:].electrons.pressure_fast_perpendicular"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

expressions["core_profiles.profiles_1d[:].electrons.pressure"] =
    (rho_tor_norm; electrons, _...) -> electrons.pressure_thermal .+ electrons.pressure_fast_parallel .+ electrons.pressure_fast_perpendicular .* 2.0


expressions["core_profiles.profiles_1d[:].ion[:].z_ion"] =
    (; ion, _...) -> begin
        if length(ion.element) == 1
            return ion.element[1].z_n
        else
            error("z_ion expression does not yet handle multiple charge states")
        end
    end

expressions["core_profiles.profiles_1d[:].ion[:].density"] =
    (rho_tor_norm; ion, _...) -> ion.density_thermal .+ ion.density_fast

expressions["core_profiles.profiles_1d[:].ion[:].density_fast"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

expressions["core_profiles.profiles_1d[:].ion[:].pressure_thermal"] =
    (rho_tor_norm; ion, _...) -> ion.temperature .* ion.density .* constants.e

expressions["core_profiles.profiles_1d[:].ion[:].pressure_fast_parallel"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

expressions["core_profiles.profiles_1d[:].ion[:].pressure_fast_perpendicular"] =
    (rho_tor_norm; _...) -> zero(rho_tor_norm)

expressions["core_profiles.profiles_1d[:]..ion[:].pressure"] =
    (rho_tor_norm; ion, _...) -> ion.pressure_thermal .+ ion.pressure_fast_parallel .+ ion.pressure_fast_perpendicular .* 2.0


expressions["core_profiles.profiles_1d[:].pressure_ion_total"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
        tmp = zero(rho_tor_norm)
        for ion in profiles_1d.ion
            tmp .+= ion.pressure_thermal
        end
        return tmp
    end

expressions["core_profiles.profiles_1d[:].pressure_thermal"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.electrons.pressure_thermal .+ profiles_1d.pressure_ion_total

expressions["core_profiles.profiles_1d[:].pressure_parallel"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.pressure_thermal ./ 3.0 .+ profiles_1d.electrons.pressure_fast_parallel .+ sum([ion.pressure_fast_parallel for ion in profiles_1d.ion])

expressions["core_profiles.profiles_1d[:].pressure_perpendicular"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.pressure_thermal ./ 3.0 .+ profiles_1d.electrons.pressure_fast_perpendicular .+ sum([ion.pressure_fast_perpendicular for ion in profiles_1d.ion])

expressions["core_profiles.profiles_1d[:].pressure"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.pressure_perpendicular .* 2.0 .+ profiles_1d.pressure_parallel


expressions["core_profiles.profiles_1d[:].conductivity_parallel"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> nclass_conductivity(dd.equilibrium.time_slice[Float64(profiles_1d.time)], profiles_1d)

expressions["core_profiles.profiles_1d[:].j_bootstrap"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> Sauter_neo2021_bootstrap(dd.equilibrium.time_slice[Float64(profiles_1d.time)], profiles_1d)

expressions["core_profiles.profiles_1d[:].j_ohmic"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.j_total .- profiles_1d.j_non_inductive

expressions["core_profiles.profiles_1d[:].j_non_inductive"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> total_sources(dd.core_sources, profiles_1d; exclude_indexes=[7, 13]).j_parallel .+ profiles_1d.j_bootstrap

expressions["core_profiles.profiles_1d[:].j_total"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> profiles_1d.j_ohmic .+ profiles_1d.j_non_inductive

expressions["core_profiles.profiles_1d[:].j_tor"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        Jpar_2_Jtor(rho_tor_norm, profiles_1d.j_total, true, eqt)
    end

expressions["core_profiles.profiles_1d[:].grid.psi_norm"] =
    (rho_tor_norm; grid, _...) -> norm01(grid.psi)

expressions["core_profiles.profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, volume, :cubic).(rho_tor_norm)
    end

expressions["core_profiles.profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, area, :cubic).(rho_tor_norm)
    end

expressions["core_profiles.profiles_1d[:].grid.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end

expressions["core_profiles.profiles_1d[:].grid.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        return interp1d(eqt.profiles_1d.rho_tor_norm, psi, :cubic).(rho_tor_norm)
    end

expressions["core_profiles.profiles_1d[:].time"] =
    (; core_profiles, profiles_1d_index, _...) -> core_profiles.time[profiles_1d_index]

expressions["core_profiles.vacuum_toroidal_field.b0"] =
    (time; dd, core_profiles, _...) -> interp1d(dd.equilibrium.time, dd.equilibrium.vacuum_toroidal_field.b0, :constant).(core_profiles.time)

expressions["core_profiles.vacuum_toroidal_field.r0"] =
    (; dd, _...) -> dd.equilibrium.vacuum_toroidal_field.r0

#= ========= =#
# equilibrium #
#= ========= =#
# IMAS does not hold B0 information in a given time slice, but we can get that info from `B0=f/R0`
# This trick propagates the B0 information to a time_slice even when that time_slice has not been initialized with profiles_1d data
expressions["equilibrium.time_slice[:].profiles_1d.f"] =
    (psi; equilibrium, time_slice_index, _...) -> (psi === missing ? [1] : ones(size(psi))) .* (equilibrium.vacuum_toroidal_field.b0[time_slice_index] * equilibrium.vacuum_toroidal_field.r0)

expressions["equilibrium.time_slice[:].global_quantities.energy_mhd"] =
    (; time_slice, _...) -> 3 / 2 * integrate(time_slice.profiles_1d.volume, time_slice.profiles_1d.pressure)

expressions["equilibrium.time_slice[:].global_quantities.q_95"] =
    (; time_slice, _...) -> Interpolations.linear_interpolation(norm01(time_slice.profiles_1d.psi), time_slice.profiles_1d.q)(0.95)

expressions["equilibrium.time_slice[:].global_quantities.q_axis"] =
    (; time_slice, _...) -> time_slice.profiles_1d.q[1]

expressions["equilibrium.time_slice[:].global_quantities.q_min"] =
    (; time_slice, _...) -> minimum(time_slice.profiles_1d.q)

expressions["equilibrium.time_slice[:].global_quantities.psi_axis"] =
    (; time_slice, _...) -> time_slice.profiles_1d.psi[1]

expressions["equilibrium.time_slice[:].global_quantities.psi_boundary"] =
    (; time_slice, _...) -> time_slice.profiles_1d.psi[end]

expressions["equilibrium.time_slice[:].global_quantities.magnetic_axis.r"] =
    (; time_slice, _...) -> time_slice.profiles_1d.geometric_axis.r[1]

expressions["equilibrium.time_slice[:].global_quantities.magnetic_axis.z"] =
    (; time_slice, _...) -> time_slice.profiles_1d.geometric_axis.z[1]

expressions["equilibrium.time_slice[:].global_quantities.magnetic_axis.b_field_tor"] =
    (; equilibrium, time_slice_index, _...) -> equilibrium.vacuum_toroidal_field.b0[time_slice_index] * equilibrium.vacuum_toroidal_field.r0 / equilibrium.time_slice[time_slice_index].boundary.geometric_axis.r

expressions["equilibrium.time_slice[:].profiles_1d.geometric_axis.r"] =
    (psi; time_slice, _...) -> (time_slice.profiles_1d.r_outboard .+ time_slice.profiles_1d.r_inboard) .* 0.5

expressions["equilibrium.time_slice[:].profiles_1d.geometric_axis.z"] =
    (psi; time_slice, _...) -> psi .* 0.0 .+ time_slice.global_quantities.magnetic_axis.z

expressions["equilibrium.time_slice[:].boundary.geometric_axis.r"] =
    (; time_slice, _...) -> time_slice.profiles_1d.geometric_axis.r[end]

expressions["equilibrium.time_slice[:].boundary.geometric_axis.z"] =
    (; time_slice, _...) -> time_slice.profiles_1d.geometric_axis.z[end]

expressions["equilibrium.time_slice[:].boundary.minor_radius"] =
    (; time_slice, _...) -> (time_slice.profiles_1d.r_outboard[end] - time_slice.profiles_1d.r_inboard[end]) * 0.5

expressions["equilibrium.time_slice[:].boundary.elongation"] =
    (; time_slice, _...) -> (time_slice.boundary.elongation_lower + time_slice.boundary.elongation_upper) * 0.5

expressions["equilibrium.time_slice[:].boundary.elongation_lower"] =
    (; time_slice, _...) -> time_slice.profiles_1d.elongation[end] # <======= IMAS limitation

expressions["equilibrium.time_slice[:].boundary.elongation_upper"] =
    (; time_slice, _...) -> time_slice.profiles_1d.elongation[end] # <======= IMAS limitation

expressions["equilibrium.time_slice[:].boundary.triangularity"] =
    (; time_slice, _...) -> (time_slice.boundary.triangularity_lower + time_slice.boundary.triangularity_upper) * 0.5

expressions["equilibrium.time_slice[:].boundary.triangularity_lower"] =
    (; time_slice, _...) -> time_slice.profiles_1d.triangularity_lower[end]

expressions["equilibrium.time_slice[:].boundary.triangularity_upper"] =
    (; time_slice, _...) -> time_slice.profiles_1d.triangularity_upper[end]

expressions["equilibrium.time_slice[:].boundary.squareness_lower_inner"] =
    (; time_slice, _...) -> time_slice.profiles_1d.squareness_lower_inner[end]

expressions["equilibrium.time_slice[:].boundary.squareness_upper_inner"] =
    (; time_slice, _...) -> time_slice.profiles_1d.squareness_upper_inner[end]

expressions["equilibrium.time_slice[:].boundary.squareness_lower_outer"] =
    (; time_slice, _...) -> time_slice.profiles_1d.squareness_lower_outer[end]

expressions["equilibrium.time_slice[:].boundary.squareness_upper_outer"] =
    (; time_slice, _...) -> time_slice.profiles_1d.squareness_upper_outer[end]

expressions["equilibrium.time_slice[:].boundary.squareness"] =
    (; time_slice, _...) -> begin
        mxh = MXH(time_slice.boundary.outline.r, time_slice.boundary.outline.z, 2)
        return -mxh.s[2]
    end

expressions["equilibrium.time_slice[:].profiles_1d.j_tor"] =
    (psi; profiles_1d, _...) -> begin
    j_parallel = profiles_1d.j_parallel
    Jpar_2_Jtor(profiles_1d.rho_tor_norm, j_parallel, true, time_slice)
end

expressions["equilibrium.time_slice[:].profiles_1d.j_parallel"] =
    (psi; time_slice, profiles_1d, _...) -> begin
    j_tor = profiles_1d.j_tor
    Jtor_2_Jpar(profiles_1d.rho_tor_norm, j_tor, true, time_slice)
end

expressions["equilibrium.time_slice[:].time"] =
    (; equilibrium, time_slice_index, _...) -> equilibrium.time[time_slice_index]

expressions["equilibrium.time_slice[:].profiles_1d.psi_norm"] =
    (psi; _...) -> norm01(psi)

#= ============ =#
#  core_sources  #
#= ============ =#
expressions["core_sources.source[:].profiles_1d[:].electrons.power_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    energy = profiles_1d.electrons.energy
    cumul_integrate(profiles_1d.grid.volume, energy)
end

expressions["core_sources.source[:].profiles_1d[:].electrons.energy"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    power_inside = profiles_1d.electrons.power_inside
    gradient(profiles_1d.grid.volume, power_inside)
end


expressions["core_sources.source[:].profiles_1d[:].total_ion_power_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    total_ion_energy = profiles_1d.total_ion_energy
    cumul_integrate(profiles_1d.grid.volume, total_ion_energy)
end

expressions["core_sources.source[:].profiles_1d[:].total_ion_energy"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    total_ion_power_inside = profiles_1d.total_ion_power_inside
    gradient(profiles_1d.grid.volume, total_ion_power_inside)
end


expressions["core_sources.source[:].profiles_1d[:].electrons.particles_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    particles = profiles_1d.electrons.particles
    cumul_integrate(profiles_1d.grid.volume, particles)
end

expressions["core_sources.source[:].profiles_1d[:].electrons.particles"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    particles_inside = profiles_1d.electrons.particles_inside
    gradient(profiles_1d.grid.volume, particles_inside)
end


expressions["core_sources.source[:].profiles_1d[:].current_parallel_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    j_parallel = profiles_1d.j_parallel
    cumul_integrate(profiles_1d.grid.area, j_parallel)
end

expressions["core_sources.source[:].profiles_1d[:].j_parallel"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    current_parallel_inside = profiles_1d.current_parallel_inside
    gradient(profiles_1d.grid.area, current_parallel_inside)
end


expressions["core_sources.source[:].profiles_1d[:].torque_tor_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    momentum_tor = profiles_1d.momentum_tor
    cumul_integrate(profiles_1d.grid.volume, momentum_tor)
end

expressions["core_sources.source[:].profiles_1d[:].momentum_tor"] =
    (rho_tor_norm; profiles_1d, _...) -> begin
    torque_tor_inside = profiles_1d.torque_tor_inside
    gradient(profiles_1d.grid.volume, torque_tor_inside)
end


expressions["core_sources.source[:].profiles_1d[:].grid.psi_norm"] =
    (rho_tor_norm; grid, _...) -> norm01(grid.psi)

expressions["core_sources.source[:].profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, volume, :cubic).(rho_tor_norm)
    end

expressions["core_sources.source[:].profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, area, :cubic).(rho_tor_norm)
    end

expressions["core_sources.source[:].profiles_1d[:].grid.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end

expressions["core_sources.source[:].profiles_1d[:].grid.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        return interp1d(eqt.profiles_1d.rho_tor_norm, psi, :cubic).(rho_tor_norm)
    end

expressions["core_sources.source[:].profiles_1d[:].time"] =
    (; core_sources, profiles_1d_index, _...) -> core_sources.time[profiles_1d_index]

expressions["core_sources.source[:].global_quantities[:].time"] =
    (; core_sources, global_quantities_index, _...) -> core_sources.time[global_quantities_index]

expressions["core_sources.vacuum_toroidal_field.b0"] =
    (time; dd, core_sources, _...) -> interp1d(dd.equilibrium.time, dd.equilibrium.vacuum_toroidal_field.b0, :constant).(core_sources.time)

expressions["core_sources.vacuum_toroidal_field.r0"] =
    (; dd, _...) -> dd.equilibrium.vacuum_toroidal_field.r0

#= ===== =#
#  build  #
#= ===== =#
expressions["build.layer[:].outline.r"] =
    (x; build, layer, _...) -> get_build(build.layer; identifier=layer.identifier, fs=(layer.fs == Int(_lfs_)) ? _hfs_ : _lfs_).outline.r

expressions["build.layer[:].outline.z"] =
    (x; build, layer, _...) -> get_build(build.layer; identifier=layer.identifier, fs=(layer.fs == Int(_lfs_)) ? _hfs_ : _lfs_).outline.z

expressions["build.layer[:].shape"] =
    (; build, layer, _...) -> get_build(build.layer; identifier=layer.identifier, fs=(layer.fs == Int(_lfs_)) ? _hfs_ : _lfs_).shape

expressions["build.layer[:].shape_parameters"] =
    (; build, layer, _...) -> get_build(build.layer; identifier=layer.identifier, fs=(layer.fs == Int(_lfs_)) ? _hfs_ : _lfs_).shape_parameters

expressions["build.layer[:].start_radius"] =
    (; build, layer_index, _...) -> build_radii(build)[1:end-1][layer_index]

expressions["build.layer[:].end_radius"] =
    (; build, layer_index, _...) -> build_radii(build)[2:end][layer_index]

expressions["build.layer[:].area"] =
    (; layer, _...) -> area(layer)

expressions["build.layer[:].volume"] =
    (; layer, _...) -> volume(layer)

expressions["build.tf.ripple"] =
    (; build, _...) -> tf_ripple(get_build(build, type=_plasma_).end_radius, get_build(build, type=_tf_, fs=_lfs_).start_radius, build.tf.coils_n)

#= ======= =#
#  costing  #
#= ======= =#
expressions["costing.cost_direct_capital.system[:].cost"] =
    (; system, _...) -> isempty(system.subsystem) ? 0.0 : sum([sub.cost for sub in system.subsystem])

expressions["costing.cost_direct_capital.cost"] =
    (; cost_direct_capital, _...) -> isempty(cost_direct_capital.system) ? 0.0 : sum([sys.cost for sys in cost_direct_capital.system])

expressions["costing.cost_operations.system[:].yearly_cost"] =
    (; system, _...) -> isempty(system.subsystem) ? 0.0 : sum([sub.yearly_cost for sub in system.subsystem])

expressions["costing.cost_operations.yearly_cost"] =
    (; cost_operations, _...) -> isempty(cost_operations.system) ? 0.0 : sum([sys.yearly_cost for sys in cost_operations.system])

expressions["costing.cost_decommissioning.system[:].cost"] =
    (; system, _...) -> isempty(system.subsystem) ? 0.0 : sum([sub.cost for sub in system.subsystem])

expressions["costing.cost_decommissioning.cost"] =
    (; cost_decommissioning, _...) -> isempty(cost_decommissioning.system) ? 0.0 : sum([sys.cost for sys in cost_decommissioning.system])

#= ============== =#
#  BalanceOfPlant  #
#= ============== =#
expressions["balance_of_plant.Q_plant"] =
    (time; balance_of_plant, _...) -> balance_of_plant.thermal_cycle.power_electric_generated ./ balance_of_plant.power_electric_plant_operation.total_power

expressions["balance_of_plant.power_electric_net"] =
    (time; balance_of_plant, _...) -> balance_of_plant.thermal_cycle.power_electric_generated .- balance_of_plant.power_electric_plant_operation.total_power

expressions["balance_of_plant.power_electric_plant_operation.total_power"] =
    (time; power_electric_plant_operation, _...) -> sum([sys.power for sys in power_electric_plant_operation.system])

expressions["balance_of_plant.thermal_cycle.power_thermal_convertable_total"] =
    (time; thermal_cycle, _...) -> sum([sys.power_in for sys in thermal_cycle.system])

expressions["balance_of_plant.thermal_cycle.power_electric_generated"] =
    (time; thermal_cycle, _...) -> thermal_cycle.thermal_electric_conversion_efficiency .* sum([sys.power_in for sys in thermal_cycle.system])

#= ======= =#
#  summary  #
#= ======= =#
expressions["summary.global_quantities.ip.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.ip for time in summary.time]

expressions["summary.global_quantities.beta_tor_mhd.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.beta_tor for time in summary.time]

expressions["summary.global_quantities.beta_pol_mhd.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.beta_pol for time in summary.time]

expressions["summary.global_quantities.beta_tor_norm_mhd.value"] =
    (time; dd, summary, _...) -> [dd.equilibrium.time_slice[Float64(time)].global_quantities.beta_normal for time in summary.time]

expressions["summary.global_quantities.current_bootstrap.value"] =
    (time; dd, summary, _...) -> begin
        tmp = []
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.area, cp1d.j_bootstrap))
        end
        return tmp
    end

expressions["summary.global_quantities.current_non_inductive.value"] =
    (time; dd, summary, _...) -> begin
        tmp = []
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.area, cp1d.j_non_inductive))
        end
        return tmp
    end

expressions["summary.global_quantities.current_ohm.value"] =
    (time; dd, summary, _...) -> begin
        tmp = []
        for time in summary.time
            cp1d = dd.core_profiles.profiles_1d[Float64(time)]
            push!(tmp, integrate(cp1d.grid.area, cp1d.j_ohmic))
        end
        return tmp
    end

expressions["summary.global_quantities.beta_tor_thermal_norm.value"] =
    (time; dd, summary, _...) -> [calc_beta_thermal_norm(dd.equilibrium, dd.core_profiles.profiles_1d[Float64(time)]) for time in summary.time]

expressions["summary.global_quantities.energy_thermal.value"] =
    (time; dd, summary, _...) -> [energy_thermal(dd.core_profiles.profiles_1d[Float64(time)]) for time in summary.time]

expressions["summary.global_quantities.tau_energy.value"] =
    (time; dd, summary, _...) -> [tau_e_thermal(dd.core_profiles.profiles_1d[Float64(time)], dd.core_sources) for time in summary.time]

expressions["summary.global_quantities.tau_energy_98.value"] =
    (time; dd, summary, _...) -> [tau_e_h98(dd, time=time) for time in summary.time]

expressions["summary.global_quantities.h_98.value"] =
    (time; dd, summary, _...) -> summary.global_quantities.tau_energy.value ./ summary.global_quantities.tau_energy_98.value
import NumericalIntegration: integrate, cumul_integrate
expressions = Dict{String,Function}()

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
# Core Profiles #
#= =========== =#
expressions["core_profiles.profiles_1d[:].electrons.pressure"] =
    (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density * constants.e

expressions["core_profiles.profiles_1d[:].electrons.density"] =
    (rho_tor_norm; electrons, _...) -> electrons.pressure ./ (electrons.temperature * constants.e)

expressions["core_profiles.profiles_1d[:].electrons.temperature"] =
    (rho_tor_norm; electrons, _...) -> electrons.pressure ./ (electrons.density * constants.e)

expressions["dd.summary.global_quantities.beta_tor_thermal_norm"] =
    (;dd, _...) -> calc_beta_thermal_norm!(dd.summary, dd.equilibrium, dd.core_profiles)

expressions["core_profiles.profiles_1d[:].pressure_thermal"] =
    (rho_tor_norm; core_profiles, _...) -> total_pressure_thermal!(core_profiles)

expressions["core_profiles.profiles_1d[:].conductivity_parallel"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> nclass_conductivity!(dd; time=profiles_1d.time)

expressions["core_profiles.profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        return interp(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume)(rho_tor_norm)
    end

expressions["core_profiles.source[:].profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        return interp(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area)(rho_tor_norm)
    end

#= ========= =#
# Equilibrium #
#= ========= =#
# IMAS does not hold B0 information in a given time slice, but we can get that info from `B0=f/R0`
# This trick propagates the B0 information to a time_slice even when that time_slice has not been initialized with profiles_1d data
expressions["equilibrium.time_slice[:].profiles_1d.f"] = (psi; equilibrium, time_slice_index, _...) -> (psi === missing ? [1] : ones(size(psi))) .* (equilibrium.vacuum_toroidal_field.b0[time_slice_index] * equilibrium.vacuum_toroidal_field.r0)

expressions["equilibrium.time_slice[:].global_quantities.energy_mhd"] =
    (;time_slice, _...) -> 3 / 2 * integrate(time_slice.profiles_1d.volume, time_slice.profiles_1d.pressure)

expressions["equilibrium.time_slice[:].global_quantities.q_95"] =
    (;time_slice, _...) -> Interpolations.LinearInterpolation(time_slice.profiles_1d.rho_tor_norm, time_slice.profiles_1d.q)(0.95)

expressions["equilibrium.time_slice[:].global_quantities.q_axis"] =
    (;time_slice, _...) -> time_slice.profiles_1d.q[1]

expressions["equilibrium.time_slice[:].global_quantities.q_min"] =
    (;time_slice, _...) -> minimum(time_slice.profiles_1d.q)

expressions["equilibrium.time_slice[:].global_quantities.psi_axis"] =
    (;time_slice, _...) -> time_slice.profiles_1d.psi[1]

expressions["equilibrium.time_slice[:].global_quantities.psi_boundary"] =
    (;time_slice, _...) -> time_slice.profiles_1d.psi[end]

expressions["equilibrium.time_slice[:].global_quantities.magnetic_axis.r"] =
    (;time_slice, _...) -> time_slice.profiles_1d.geometric_axis.r[1]

expressions["equilibrium.time_slice[:].global_quantities.magnetic_axis.z"] =
    (;time_slice, _...) -> time_slice.profiles_1d.geometric_axis.z[1]

expressions["equilibrium.time_slice[:].global_quantities.magnetic_axis.b_field_tor"] =
    (; equilibrium, _...) ->  @ddtime(equilibrium.vacuum_toroidal_field.b0) * equilibrium.vacuum_toroidal_field.r0 / equilibrium.boundary.geometric_axis.r

expressions["equilibrium.time_slice[:].profiles_1d.geometric_axis.r"] =
    (psi; time_slice, _...) -> (time_slice.profiles_1d.r_outboard .+ time_slice.profiles_1d.r_inboard) .* 0.5

expressions["equilibrium.time_slice[:].profiles_1d.geometric_axis.z"] =
    (psi; time_slice, _...) -> psi .* 0.0 .+ time_slice.global_quantities.magnetic_axis.z

expressions["equilibrium.time_slice[:].boundary.geometric_axis.r"] =
    (;time_slice, _...) -> time_slice.profiles_1d.geometric_axis.r[end]

expressions["equilibrium.time_slice[:].boundary.geometric_axis.z"] =
    (;time_slice, _...) -> time_slice.profiles_1d.geometric_axis.z[end]

expressions["equilibrium.time_slice[:].boundary.minor_radius"] =
    (;time_slice, _...) -> (time_slice.profiles_1d.r_outboard[end] - time_slice.profiles_1d.r_inboard[end]) * 0.5

expressions["equilibrium.time_slice[:].boundary.elongation"] =
    (;time_slice, _...) -> (time_slice.boundary.elongation_lower + time_slice.boundary.elongation_upper) * 0.5

expressions["equilibrium.time_slice[:].boundary.elongation_lower"] =
    (;time_slice, _...) -> time_slice.profiles_1d.elongation[end] # <======= IMAS 3.30.0 limitation

expressions["equilibrium.time_slice[:].boundary.elongation_upper"] =
    (;time_slice, _...) -> time_slice.profiles_1d.elongation[end] # <======= IMAS 3.30.0 limitation

expressions["equilibrium.time_slice[:].boundary.triangularity"] =
    (;time_slice, _...) -> (time_slice.boundary.triangularity_lower + time_slice.boundary.triangularity_upper) * 0.5

expressions["equilibrium.time_slice[:].boundary.triangularity_lower"] =
    (;time_slice, _...) -> time_slice.profiles_1d.triangularity_lower[end]

expressions["equilibrium.time_slice[:].boundary.triangularity_upper"] =
    (;time_slice, _...) -> time_slice.profiles_1d.triangularity_upper[end]

expressions["equilibrium.time_slice[:].boundary.squareness_lower_inner"] =
    (;time_slice, _...) -> time_slice.profiles_1d.squareness_lower_inner[end]

expressions["equilibrium.time_slice[:].boundary.squareness_upper_inner"] =
    (;time_slice, _...) -> time_slice.profiles_1d.squareness_upper_inner[end]

expressions["equilibrium.time_slice[:].boundary.squareness_lower_outer"] =
    (;time_slice, _...) -> time_slice.profiles_1d.squareness_lower_outer[end]

expressions["equilibrium.time_slice[:].boundary.squareness_upper_outer"] =
    (;time_slice, _...) -> time_slice.profiles_1d.squareness_upper_outer[end]

#= ============ =#
#  core_sources  #
#= ============ =#
expressions["core_sources.source[:].profiles_1d[:].electrons.power_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> cumul_integrate(profiles_1d.grid.volume, profiles_1d.electrons.energy)

expressions["core_sources.source[:].profiles_1d[:].electrons.energy"] =
    (rho_tor_norm; profiles_1d, _...) -> gradient(profiles_1d.grid.volume, profiles_1d.electrons.power_inside)


expressions["core_sources.source[:].profiles_1d[:].total_ion_power_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> cumul_integrate(profiles_1d.grid.volume, profiles_1d.total_ion_energy)

expressions["core_sources.source[:].profiles_1d[:].total_ion_energy"] =
    (rho_tor_norm; profiles_1d, _...) -> gradient(profiles_1d.grid.volume, profiles_1d.total_ion_power_inside)


expressions["core_sources.source[:].profiles_1d[:].electrons.particles_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> cumul_integrate(profiles_1d.grid.volume, profiles_1d.electrons.particles)

expressions["core_sources.source[:].profiles_1d[:].electrons.particles"] =
    (rho_tor_norm; profiles_1d, _...) -> gradient(profiles_1d.grid.volume, profiles_1d.electrons.particles_inside)


expressions["core_sources.source[:].profiles_1d[:].current_parallel_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> cumul_integrate(profiles_1d.grid.volume, profiles_1d.j_parallel)

expressions["core_sources.source[:].profiles_1d[:]. "] =
    (rho_tor_norm; profiles_1d, _...) -> gradient(profiles_1d.grid.volume, profiles_1d.current_parallel_inside)


expressions["core_sources.source[:].profiles_1d[:].torque_tor_inside"] =
    (rho_tor_norm; profiles_1d, _...) -> cumul_integrate(profiles_1d.grid.volume, profiles_1d.momentum_tor)

expressions["core_sources.source[:].profiles_1d[:].momentum_tor"] =
    (rho_tor_norm; profiles_1d, _...) -> gradient(profiles_1d.grid.volume, profiles_1d.torque_tor_inside)


expressions["core_sources.source[:].profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        return interp(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume)(rho_tor_norm)
    end

expressions["core_sources.source[:].profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        return interp(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area)(rho_tor_norm)
    end


#= ===== =#
#  Build  #
#= ===== =#
expressions["build.layer[:].start_radius"] =
    (;build, layer_index, _...) -> build_radii(build)[1:end - 1][layer_index]

expressions["build.layer[:].end_radius"] =
    (;build, layer_index, _...) -> build_radii(build)[2:end][layer_index]

document[:Expressions] = Symbol[]

function IMASdd.get_expressions(::Type{Val{:onetime}})
    return onetime_expressions
end

const onetime_expressions = Dict{String,Function}()
otexp = onetime_expressions

# These expressions are frozen the first time they are accessed.
# This is necessary to ensure that core_profiles, core_sources, and core_transport grids do not change after changing the equilibrium.
# The idea is that we want to freeze in the core_profiles, core_sources, and core_transport grids the rho, psi, volume, area, ... info that were used when those IDSs were filled.
# While this is generally ok, this is not desirable when iterating the equilibrium solver with other actors.
# In this case, at each iteration we want core_profiles, core_sources, and core_transport to take the grids from the latest equilibrium.
#
# NOTE: make sure that expressions accept as argument (not keyword argument)
# the coordinates of the quantitiy you are writing the expression of
#
# For example, this will FAIL:
#    otexp["core_profiles.profiles_1d[:].electrons.pressure_thermal"] =
#         (; electrons, _...) -> electrons.temperature .* electrons.density_thermal * mks.e
#
# This is GOOD:
#    otexp["core_profiles.profiles_1d[:].electrons.pressure_thermal"] =
#         (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density_thermal * mks.e

#= =========== =#
# core_profiles #
#= =========== =#
otexp["core_profiles.profiles_1d[:].grid.psi_norm"] =
    (rho_tor_norm; grid, _...) -> norm01(grid.psi)

otexp["core_profiles.profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(volume), :cubic).(rho_tor_norm) .^ 2
    end

otexp["core_profiles.profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(area), :cubic).(rho_tor_norm) .^ 2
    end

otexp["core_profiles.profiles_1d[:].grid.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end


otexp["core_profiles.profiles_1d[:].grid.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        sign_psi = sign(psi[end] - psi[1])
        return sign_psi .* (interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(abs.(psi .- psi[1])), :cubic).(rho_tor_norm) .^ 2) .+ psi[1]
    end

#= ============ =#
# core_transport #
#= ============ =#
otexp["core_transport.model[:].profiles_1d[:].grid_flux.rho_tor_norm"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.grid_d.rho_tor_norm

otexp["core_transport.model[:].profiles_1d[:].grid_d.rho_tor_norm"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.grid_flux.rho_tor_norm

otexp["core_transport.model[:].profiles_1d[:].grid_flux.psi_norm"] =
    (rho_tor_norm; grid_flux, _...) -> norm01(grid_flux.psi)

otexp["core_transport.model[:].profiles_1d[:].grid_flux.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, volume, :cubic).(rho_tor_norm)
    end

otexp["core_transport.model[:].profiles_1d[:].grid_flux.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, area, :cubic).(rho_tor_norm)
    end

otexp["core_transport.model[:].profiles_1d[:].grid_flux.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end

otexp["core_transport.model[:].profiles_1d[:].grid_flux.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        return interp1d(eqt.profiles_1d.rho_tor_norm, psi, :cubic).(rho_tor_norm)
    end

#= ============ =#
#  core_sources  #
#= ============ =#
otexp["core_sources.source[:].profiles_1d[:].grid.psi_norm"] =
    (rho_tor_norm; grid, _...) -> norm01(grid.psi)

otexp["core_sources.source[:].profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, volume, :cubic).(rho_tor_norm)
    end

otexp["core_sources.source[:].profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, area, :cubic).(rho_tor_norm)
    end

otexp["core_sources.source[:].profiles_1d[:].grid.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end

otexp["core_sources.source[:].profiles_1d[:].grid.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        return interp1d(eqt.profiles_1d.rho_tor_norm, psi, :cubic).(rho_tor_norm)
    end

#= ===== =#
#  waves  #
#= ===== =#
otexp["waves.coherent_wave[:].profiles_1d[:].grid.psi_norm"] =
    (rho_tor_norm; grid, _...) -> norm01(grid.psi)

otexp["waves.coherent_wave[:].profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, volume, :cubic).(rho_tor_norm)
    end

otexp["waves.coherent_wave[:].profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, area, :cubic).(rho_tor_norm)
    end

otexp["waves.coherent_wave[:].profiles_1d[:].grid.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end

otexp["waves.coherent_wave[:].profiles_1d[:].grid.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        return interp1d(eqt.profiles_1d.rho_tor_norm, psi, :cubic).(rho_tor_norm)
    end

# ============ #

Base.Docs.@doc """
    onetime_expressions = Dict{String,Function}()

Expressions that are frozen after first evaluation
* `$(join(sort!(collect(keys(onetime_expressions))),"`\n* `"))`
""" onetime_expressions

push!(document[:Expressions], :onetime_expressions)
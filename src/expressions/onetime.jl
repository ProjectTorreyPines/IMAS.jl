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

#= =========== =#
# core_profiles #
#= =========== =#
otexp["core_profiles.profiles_1d[:].grid.psi_norm"] =
    (; grid, _...) -> norm01(grid.psi)

otexp["core_profiles.profiles_1d[:].grid.volume"] =
    (; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        volume = eqt.profiles_1d.volume
        vitp =  cubic_interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(volume))
        return @. vitp(profiles_1d.grid.rho_tor_norm) ^ 2
    end

otexp["core_profiles.profiles_1d[:].grid.area"] =
    (; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        area = eqt.profiles_1d.area
        aitp = cubic_interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(area))
        return @. aitp(profiles_1d.grid.rho_tor_norm) ^ 2
    end

otexp["core_profiles.profiles_1d[:].grid.surface"] =
    (; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        surface = eqt.profiles_1d.surface
        sitp = cubic_interp1d(eqt.profiles_1d.rho_tor_norm, surface)
        return sitp.(profiles_1d.grid.rho_tor_norm)
    end


otexp["core_profiles.profiles_1d[:].grid.psi"] =
    (; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        psi = eqt.profiles_1d.psi
        sign_psi = sign(psi[end] - psi[1])
        pitp = cubic_interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(abs.(psi .- psi[1])))
        return @. sign_psi * (pitp(profiles_1d.grid.rho_tor_norm) ^ 2) + psi[1]
    end


otexp["core_profiles.profiles_1d[:].electrons.pressure_fast_parallel"] =
    (; profiles_1d, _...) -> zero(profiles_1d.grid.rho_tor_norm)

otexp["core_profiles.profiles_1d[:].electrons.pressure_fast_perpendicular"] =
    (; profiles_1d, _...) -> zero(profiles_1d.grid.rho_tor_norm)

otexp["core_profiles.profiles_1d[:].ion[:].pressure_fast_parallel"] =
    (; profiles_1d, _...) -> zero(profiles_1d.grid.rho_tor_norm)

otexp["core_profiles.profiles_1d[:].ion[:].pressure_fast_perpendicular"] =
    (; profiles_1d, _...) -> zero(profiles_1d.grid.rho_tor_norm)

#= ============ =#
# core_transport #
#= ============ =#
otexp["core_transport.model[:].profiles_1d[:].grid_flux.rho_tor_norm"] =
    (; profiles_1d, _...) -> profiles_1d.grid_d.rho_tor_norm

otexp["core_transport.model[:].profiles_1d[:].grid_d.rho_tor_norm"] =
    (; profiles_1d, _...) -> profiles_1d.grid_flux.rho_tor_norm

otexp["core_transport.model[:].profiles_1d[:].grid_flux.psi_norm"] =
    (; grid_flux, _...) -> norm01(grid_flux.psi)

otexp["core_transport.model[:].profiles_1d[:].grid_flux.volume"] =
    (; dd, profiles_1d, grid_flux, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        volume = eqt.profiles_1d.volume
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, volume).(grid_flux.rho_tor_norm)
    end

otexp["core_transport.model[:].profiles_1d[:].grid_flux.area"] =
    (; dd, profiles_1d, grid_flux, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        area = eqt.profiles_1d.area
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, area).(grid_flux.rho_tor_norm)
    end

otexp["core_transport.model[:].profiles_1d[:].grid_flux.surface"] =
    (; dd, profiles_1d, grid_flux, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        surface = eqt.profiles_1d.surface
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, surface).(grid_flux.rho_tor_norm)
    end

otexp["core_transport.model[:].profiles_1d[:].grid_flux.psi"] =
    (; dd, profiles_1d, grid_flux, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        psi = eqt.profiles_1d.psi
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, psi).(grid_flux.rho_tor_norm)
    end

#= ============ =#
#  core_sources  #
#= ============ =#
otexp["core_sources.source[:].profiles_1d[:].grid.psi_norm"] =
    (; grid, _...) -> norm01(grid.psi)

otexp["core_sources.source[:].profiles_1d[:].grid.volume"] =
    (; dd, profiles_1d, grid, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        volume = eqt.profiles_1d.volume
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, volume).(grid.rho_tor_norm)
    end

otexp["core_sources.source[:].profiles_1d[:].grid.area"] =
    (; dd, profiles_1d, grid, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        area = eqt.profiles_1d.area
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, area).(grid.rho_tor_norm)
    end

otexp["core_sources.source[:].profiles_1d[:].grid.surface"] =
    (; dd, profiles_1d, grid, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        surface = eqt.profiles_1d.surface
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, surface).(grid.rho_tor_norm)
    end

otexp["core_sources.source[:].profiles_1d[:].grid.psi"] =
    (; dd, profiles_1d, grid, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        psi = eqt.profiles_1d.psi
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, psi).(grid.rho_tor_norm)
    end

#= ===== =#
#  waves  #
#= ===== =#
otexp["waves.coherent_wave[:].profiles_1d[:].grid.psi_norm"] =
    (; grid, _...) -> norm01(grid.psi)

otexp["waves.coherent_wave[:].profiles_1d[:].grid.volume"] =
    (; dd, profiles_1d, grid, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        volume = eqt.profiles_1d.volume
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, volume).(grid.rho_tor_norm)
    end

otexp["waves.coherent_wave[:].profiles_1d[:].grid.area"] =
    (; dd, profiles_1d, grid, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        area = eqt.profiles_1d.area
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, area).(grid.rho_tor_norm)
    end

otexp["waves.coherent_wave[:].profiles_1d[:].grid.surface"] =
    (; dd, profiles_1d, grid, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        surface = eqt.profiles_1d.surface
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, surface).(grid.rho_tor_norm)
    end

otexp["waves.coherent_wave[:].profiles_1d[:].grid.psi"] =
    (; dd, profiles_1d, grid, _...) -> begin
        eqt = dd.equilibrium.time_slice[profiles_1d.time]
        psi = eqt.profiles_1d.psi
        return cubic_interp1d(eqt.profiles_1d.rho_tor_norm, psi).(grid.rho_tor_norm)
    end

# ============ #

Base.Docs.@doc """
    onetime_expressions = Dict{String,Function}()

Expressions that are frozen after first evaluation
* `$(join(sort!(collect(keys(onetime_expressions))),"`\n* `"))`
""" onetime_expressions

push!(document[:Expressions], :onetime_expressions)

using Plots

document[Symbol("Benchmark")] = Symbol[]

"""
    Benchmark{T}

Struct for benchmarking time-dependent dd data structures
"""
struct Benchmark{T}
    dd::IMAS.dd{T}
    data::Dict{String,Vector{T}}
end

@compat public Benchmark
push!(document[Symbol("Benchmark")], :Benchmark)

"""
    Benchmark(dd::IMAS.dd{T}) where {T<:Real}

Construct Benchmark from frozen dd, extracting all leaf data
"""
function Benchmark(dd::IMAS.dd{T}) where {T<:Real}
    @assert IMAS.isfrozen(dd)

    data = Dict{String,Vector{T}}()
    leaves = collect(IMAS.IMASdd.leaves(dd))
    for leaf in leaves
        if typeof(leaf) <: IMASnodeRepr
            tloc = IMAS.IMASdd.utlocation(leaf.ids, leaf.field)
            if tloc ∉ keys(data)
                data[tloc] = Float64[]
            end
            append!(data[tloc], leaf.value)
        end
    end

    # keep only one time array
    tmp = []
    for item in collect(keys(data))
        if endswith(item, ".time")
            push!(tmp, pop!(data, item))
        end
    end
    if !isempty(tmp)
        data["time"] = tmp[1]
    end

    return Benchmark(dd, data)
end

"""
    benchmark(dd_ben::IMAS.dd{T}, times::AbstractVector{Float64}) where {T<:Real}

Evaluate self-consistency of a given IDS (flux-matching errors, etc.)
"""
function benchmark(dd_ben::IMAS.dd{T}, times::AbstractVector{Float64}) where {T<:Real}
    dd_cst = IMAS.dd(; frozen=true)

    for time0 in times
        # core transport: compare total fluxes vs fluxes calculated total sources
        if !isempty(dd_ben.core_transport.model) && !isempty(dd_ben.core_sources.source)
            # Get total fluxes from benchmark
            ct1d_ben_total = IMAS.total_fluxes(dd_ben; time0)

            # Get total sources from reference and calculate fluxes
            cs1d_ben_total = IMAS.total_sources(dd_ben; time0)
            resize!(dd_cst.core_transport.model, 1; wipe=false)
            ct1d_cst = resize!(dd_cst.core_transport.model[1].profiles_1d, time0)
            resize!(ct1d_cst.ion, length(ct1d_ben_total.ion))

            # Set up grids
            rho_ben = ct1d_ben_total.grid_flux.rho_tor_norm
            rho_flux = ct1d_ben_total.grid_flux.rho_tor_norm
            rho_source = cs1d_ben_total.grid.rho_tor_norm

            surface_interp = IMAS.interp1d(rho_source, cs1d_ben_total.grid.surface)

            # electrons energy flux
            if !ismissing(ct1d_ben_total.electrons, :energy) && !ismissing(cs1d_ben_total.electrons, :energy)
                rho0 = rho_ben[ct1d_ben_total.electrons.energy.flux.!=0.0]
                flux_from_source = IMAS.interp1d(rho_source, cs1d_ben_total.electrons.power_inside).(rho0) ./ surface_interp.(rho0)
                ct1d_cst.electrons.energy.flux = [benchmark(
                    rho_flux, ct1d_ben_total.electrons.energy.flux,
                    rho0, flux_from_source)]
            end

            # electrons particles flux
            if !ismissing(ct1d_ben_total.electrons, :particles) && !ismissing(cs1d_ben_total.electrons, :particles)
                rho0 = rho_ben[ct1d_ben_total.electrons.particles.flux.!=0.0]
                flux_from_source = IMAS.interp1d(rho_source, cs1d_ben_total.electrons.particles_inside).(rho0) ./ surface_interp.(rho0)
                ct1d_cst.electrons.particles.flux = [benchmark(
                    rho_flux, ct1d_ben_total.electrons.particles.flux,
                    rho0, flux_from_source)]
            end

            # total ion energy flux
            if !ismissing(ct1d_ben_total, :total_ion_energy) && !ismissing(cs1d_ben_total, :total_ion_energy)
                rho0 = rho_ben[ct1d_ben_total.total_ion_energy.flux.!=0.0]
                flux_from_source = IMAS.interp1d(rho_source, cs1d_ben_total.total_ion_power_inside).(rho0) ./ surface_interp.(rho0)
                ct1d_cst.total_ion_energy.flux = [benchmark(
                    rho_flux, ct1d_ben_total.total_ion_energy.flux,
                    rho0, flux_from_source)]
            end

            # momentum torque flux
            if !ismissing(ct1d_ben_total, :momentum_tor) && !ismissing(cs1d_ben_total, :momentum_tor)
                rho0 = rho_ben[ct1d_ben_total.momentum_tor.flux.!=0.0]
                flux_from_source = IMAS.interp1d(rho_source, cs1d_ben_total.torque_tor_inside).(rho0) ./ surface_interp.(rho0)
                ct1d_cst.momentum_tor.flux = [benchmark(
                    rho_flux, ct1d_ben_total.momentum_tor.flux,
                    rho0, flux_from_source)]
            end

            # ion particles flux for each ion species
            for kion in 1:length(ct1d_ben_total.ion)
                if !ismissing(ct1d_ben_total.ion[kion], :particles) && !ismissing(cs1d_ben_total.ion[kion], :particles)
                    rho0 = rho_ben[ct1d_ben_total.ion[kion].particles.flux.!=0.0]
                    flux_from_source = IMAS.interp1d(rho_source, cs1d_ben_total.ion[kion].particles_inside).(rho0) ./ surface_interp.(rho0)
                    ct1d_cst.ion[kion].particles.flux = [benchmark(
                        rho_flux, ct1d_ben_total.ion[kion].particles.flux,
                        rho0, flux_from_source)]
                end
            end
        end
    end

    return Benchmark(dd_cst)
end

"""
    benchmark(dd_ben::IMAS.dd{T}, dd_ref::IMAS.dd{T}, times::AbstractVector{Float64}; self_benchmark::Bool=true) where {T<:Real}

Benchmark two dd structures at specified times
"""
function benchmark(dd_ben::IMAS.dd{T}, dd_ref::IMAS.dd{T}, times::AbstractVector{Float64}; self_benchmark::Bool=true) where {T<:Real}
    if self_benchmark
        dd_cst = benchmark(dd_ben, times).dd
    else
        dd_cst = IMAS.dd(; frozen=true)
    end

    for time0 in times
        # equilibrium
        eqt_cst = resize!(dd_cst.equilibrium.time_slice, time0)
        eqt1d_ben = dd_ben.equilibrium.time_slice[time0].profiles_1d
        eqt1d_ref = dd_ref.equilibrium.time_slice[time0].profiles_1d
        eqt1d_cst = eqt_cst.profiles_1d

        eqt1d_cst.pressure = [benchmark(
            eqt1d_ben.psi_norm, eqt1d_ben.pressure,
            eqt1d_ref.psi_norm, eqt1d_ref.pressure)]
        eqt1d_cst.psi = [benchmark(
            eqt1d_ben.psi_norm, eqt1d_ben.psi,
            eqt1d_ref.psi_norm, eqt1d_ref.psi)]
        eqt1d_cst.j_tor = [benchmark(
            eqt1d_ben.psi_norm, eqt1d_ben.j_tor,
            eqt1d_ref.psi_norm, eqt1d_ref.j_tor)]

        # core profiles
        cp1d_ben = dd_ben.core_profiles.profiles_1d[time0]
        cp1d_ref = dd_ref.core_profiles.profiles_1d[time0]
        cp1d_cst = resize!(dd_cst.core_profiles.profiles_1d, time0)

        cp1d_cst.electrons.density_thermal = [benchmark(
            cp1d_ben.grid.rho_tor_norm, cp1d_ben.electrons.density_thermal,
            cp1d_ref.grid.rho_tor_norm, cp1d_ref.electrons.density_thermal)]
        cp1d_cst.electrons.temperature = [benchmark(
            cp1d_ben.grid.rho_tor_norm, cp1d_ben.electrons.temperature,
            cp1d_ref.grid.rho_tor_norm, cp1d_ref.electrons.temperature)]
        cp1d_cst.t_i_average = [benchmark(
            cp1d_ben.grid.rho_tor_norm, cp1d_ben.t_i_average,
            cp1d_ref.grid.rho_tor_norm, cp1d_ref.t_i_average)]
        cp1d_cst.rotation_frequency_tor_sonic = [benchmark(
            cp1d_ben.grid.rho_tor_norm, cp1d_ben.rotation_frequency_tor_sonic,
            cp1d_ref.grid.rho_tor_norm, cp1d_ref.rotation_frequency_tor_sonic)]
        resize!(cp1d_cst.ion, length(cp1d_ref.ion))
        for kion in 1:length(cp1d_ref.ion)
            cp1d_cst.ion[kion].density_thermal = [benchmark(
                cp1d_ben.grid.rho_tor_norm, cp1d_ben.ion[kion].density_thermal,
                cp1d_ref.grid.rho_tor_norm, cp1d_ref.ion[kion].density_thermal)]
        end
        cp1d_cst.pressure = [benchmark(
            cp1d_ben.grid.rho_tor_norm, cp1d_ben.pressure,
            cp1d_ref.grid.rho_tor_norm, cp1d_ref.pressure)]
    end

    return Benchmark(dd_cst)
end

"""
    benchmark(x::AbstractVector{T}, y::AbstractVector{T}, x_ref::AbstractVector{T}, y_ref::AbstractVector{T}) where {T<:Real}

evaluate error between two 1D arrays
"""
function benchmark(x::AbstractVector{T}, y::AbstractVector{T}, x_ref::AbstractVector{T}, y_ref::AbstractVector{T}) where {T<:Real}
    return norm((y_ref .- IMAS.interp1d(x, y).(x_ref)) ./ norm(y))
end

@compat public benchmark
push!(document[Symbol("Benchmark")], :benchmark)

# plot Benchmark structure
@recipe function plot_benchmark(bnch::Benchmark)
    items = sort!(collect(filter(x -> x != "time", keys(bnch.data))))
    time = bnch.data["time"]
    layout := length(items)
    size --> (1000, 1000)
    titlefontsize --> 5
    for (k, utloc) in enumerate(items)
        @series begin
            subplot := k
            label := ""
            time, bnch.data[utloc]
        end
        @series begin
            primary := false
            title := utloc
            subplot := k
            seriestype := :hline
            color := :black
            [1.0]
        end
    end
end

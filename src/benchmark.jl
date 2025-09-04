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
    for leaf in IMAS.IMASdd.leaves(dd)
        uloc = IMAS.ulocation(leaf.ids, leaf.field)
        if uloc ∉ keys(data)
            # println(uloc)
            data[uloc] = Float64[]
        end
        append!(data[uloc], leaf.value)
    end

    return Benchmark(dd, data)
end

@recipe function plot_benchmark(bnch::Benchmark)
    items = collect(filter(x -> !endswith(x, ".time"), keys(bnch.data)))
    layout := length(items)
    size --> (1000, 1000)
    titlefontsize --> 8
    for (k, uloc) in enumerate(items)
        @series begin
            subplot := k
            label := ""
            bnch.data["equilibrium.time"], bnch.data[uloc]
        end
        @series begin
            primary := false
            title := uloc
            subplot := k
            seriestype := :hline
            color := :black
            [1.0]
        end
    end
end

"""
    benchmark(dd_ben::IMAS.dd{T}, dd_ref::IMAS.dd{T}, times::AbstractVector{Float64}) where {T<:Real}

Benchmark two dd structures at specified times
"""
function benchmark(dd_ben::IMAS.dd{T}, dd_ref::IMAS.dd{T}, times::AbstractVector{Float64}) where {T<:Real}
    dd_cst = IMAS.dd(; frozen=true)

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
        cp1d_cst.pressure = [benchmark(
            cp1d_ben.grid.rho_tor_norm, cp1d_ben.pressure,
            cp1d_ref.grid.rho_tor_norm, cp1d_ref.pressure)]
    end

    return Benchmark(dd_cst)
end

@compat public benchmark
push!(document[Symbol("Benchmark")], :benchmark)

function benchmark(x::AbstractVector{T}, y::AbstractVector{T}, x_ref::AbstractVector{T}, y_ref::AbstractVector{T}) where {T<:Real}
    return norm((y_ref .- IMAS.interp1d(x, y).(x_ref)) ./ norm(y_ref))
end

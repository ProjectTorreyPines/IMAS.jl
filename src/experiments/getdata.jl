"""
    getdata(
        what_val::Union{Val{:t_i},Val{:n_i_over_n_e},Val{:zeff},Val{:n_imp},Val{:ω_tor}},
        dd::IMAS.dd{T},
        time0::Union{Nothing,Float64}=nothing,
        time_averaging::Float64=0.0
    ) where {T<:Real}

Extract charge exchange spectroscopy data for ion temperature, impurity to electron density ratio, Zeff or impurity density

  - time0: optional specific time for extraction

  - time_averaging: time averaging window size (required if time0 is specified)

Returns NamedTuple with (:time, :rho, :data, :weights, :units)

Data is always of type Measurements.Measurement{T}
"""
function getdata(
    what_val::Union{Val{:t_i},Val{:n_i_over_n_e},Val{:zeff},Val{:n_imp},Val{:ω_tor}},
    dd::IMAS.dd{T},
    time0::Union{Nothing,Float64}=nothing,
    time_averaging::Float64=0.0
) where {T<:Real}

    what = typeof(what_val).parameters[1]
    cer = dd.charge_exchange

    if time0 !== nothing
        @assert time_averaging > 0.0
    end

    random_seed = 0
    rng = Random.MersenneTwister(random_seed)
    data = Measurements.Measurement{T}[]
    weights = Float64[]
    times = Float64[]
    chr = T[]
    chz = T[]
    for ch in cer.channel
        if what == :n_imp
            ch_data = ch.ion[1].n_i_over_n_e
        elseif what ∈ (:v_tor, :ω_tor)
            ch_data = ch.ion[1].velocity_tor
        elseif what == :zeff
            ch_data = ch.zeff
        else
            ch_data = getproperty(ch.ion[1], what)
        end

        if hasdata(ch_data, :data)
            position_r = interp1d(ch.position.r.time, ch.position.r.data)
            position_z = interp1d(ch.position.z.time, ch.position.z.data)
            if time0 !== nothing
                selection = select_time_window(ch_data, :data, time0; window_size=time_averaging)
                v = selection.data
                if hasdata(ch_data, :data_σ)
                    σ = getproperty(ch_data, :data_σ)[selection.idx_range]
                else
                    σ = v .* 0.0
                end
                _data = Measurements.measurement.(v, σ)
                if what == :ω_tor
                    _data .= _data ./ position_r.(selection.time)
                end
                append!(data, _data)
                append!(times, selection.time)
                append!(weights, selection.weights)
                append!(chr, position_r.(selection.time))
                append!(chz, position_z.(selection.time))
            else
                time = getproperty(ch_data, :time)
                v = getproperty(ch_data, :data)
                if hasdata(ch_data, :data_σ)
                    σ = getproperty(ch_data, :data_σ)
                else
                    σ = v .* 0.0
                end
                _data = Measurements.measurement.(v, σ)
                if what == :ω_tor
                    _data .= _data ./ position_r.(time)
                end
                n = length(_data)
                append!(data, _data)
                append!(times, time)
                append!(chr, position_r.(time) .+ randn(rng, n) .* 1E-6)
                append!(chz, position_z.(time) .+ randn(rng, n) .* 1E-6)
                # note we add a small random variation on top of the R and Z channels
                # because in some situations (eg. DIII-D shot 200000) different CER
                # channels can have the same exact spatial location, which messes
                # up the DelaunayTriangulation in the NaturalNeighbours interpolator.
            end
        end
    end

    rho = chz .* 0.0
    for time0 in unique(times)
        i = nearest_causal_time(dd.equilibrium.time, time0; bounds_error=false).index
        eqt = dd.equilibrium.time_slice[i]
        r, z, RHO_interpolant = ρ_interpolant(eqt)
        index = (times .== time0)
        rho[index] .= RHO_interpolant.(chr[index], chz[index])
    end

    if what == :n_imp
        cp1d = top_dd(cer).core_profiles.profiles_1d[time0]
        n_e = interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(rho)
        data = data .* n_e
    end

    if what == :t_i
        units = "eV"
    elseif what ∈ (:n_i_over_n_e, :zeff)
        units = ""
    elseif what == :n_imp
        units = "m^-3"
    elseif what == :ω_tor
        units = "rad s^-1"
    else
        error("getdata() should not be here")
    end

    return (time=times, rho=rho, data=data, weights=weights, units=units)
end

"""
    getdata(what_val::Union{Val{:t_e},Val{:n_e}}, dd::IMAS.dd{T}, time0::Union{Nothing,Float64}=nothing, time_averaging::Float64=0.0) where {T<:Real}

Extract Thomson scattering data for electron temperature or density.

  - time0: optional specific time for extraction

  - time_averaging: time averaging window size (required if time0 is specified)

Returns NamedTuple with (:time, :rho, :data, :weights)

Data is always of type Measurements.Measurement{T}
"""
function getdata(what_val::Union{Val{:t_e},Val{:n_e}}, dd::IMAS.dd{T}, time0::Union{Nothing,Float64}=nothing, time_averaging::Float64=0.0) where {T<:Real}
    what = typeof(what_val).parameters[1]
    ts = dd.thomson_scattering

    if time0 !== nothing
        @assert time_averaging > 0.0
    end

    data = Measurements.Measurement{T}[]
    weights = Float64[]
    times = Float64[]
    chr = T[]
    chz = T[]
    for ch in ts.channel
        ch_data = getproperty(ch, what)
        if time0 !== nothing
            selection = select_time_window(ch_data, :data, time0; window_size=time_averaging)
            if hasdata(ch_data, :data_σ)
                append!(data, Measurements.measurement.(selection.data, getproperty(ch_data, :data_σ)[selection.idx_range]))
            else
                append!(data, Measurements.measurement.(selection.data, selection.data .* 0.0))
            end
            append!(times, selection.time)
            append!(weights, selection.weights)
            append!(chr, fill(ch.position.r, length(selection.time)))
            append!(chz, fill(ch.position.z, length(selection.time)))
        else
            _data = getproperty(ch_data, :data)
            if hasdata(ch_data, :data_σ)
                _data = Measurements.measurement.(_data, getproperty(ch_data, :data_σ))
            else
                _data = Measurements.measurement.(_data, _data .* 0.0)
            end
            append!(data, _data)
            append!(times, getproperty(ch_data, :time))
            append!(chr, fill(ch.position.r, size(_data)))
            append!(chz, fill(ch.position.z, size(_data)))
        end
    end

    rho = chz .* 0.0
    for t in unique(times)
        i = nearest_causal_time(dd.equilibrium.time, t; bounds_error=false).index
        eqt = dd.equilibrium.time_slice[i]
        r, z, RHO_interpolant = ρ_interpolant(eqt)
        index = (times .== t)
        rho[index] = RHO_interpolant.(chr[index], chz[index])
    end

    if what == :t_e
        units = "eV"
    elseif what == :n_e
        units = "m^-3"
    else
        error("getdata() should not be here")
    end

    return (time=times, rho=rho, data=data, weights=weights, units=units)
end

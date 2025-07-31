"""
    getdata(
        what_val::Union{Val{:t_i},Val{:n_i_over_n_e},Val{:zeff},Val{:n_imp}},
        dd::IMAS.dd{T},
        time0::Union{Nothing,Float64}=nothing,
        time_averaging::Float64=0.0
    ) where {T<:Real}

Extract charge exchange spectroscopy data for ion temperature, impurity to electron density ratio, Zeff or impurity density

  - time0: optional specific time for extraction

  - time_averaging: time averaging window size (required if time0 is specified)

Returns NamedTuple with (:time, :rho, :data, :weights)

Data is always of type Measurements.Measurement
"""
function getdata(
    what_val::Union{Val{:t_i},Val{:n_i_over_n_e},Val{:zeff},Val{:n_imp}},
    dd::IMAS.dd{T},
    time0::Union{Nothing,Float64}=nothing,
    time_averaging::Float64=0.0
) where {T<:Real}

    what = typeof(what_val).parameters[1]
    cer = dd.charge_exchange

    if time0 !== nothing
        @assert time_averaging > 0.0
    end

    data = Measurements.Measurement[]
    weights = Float64[]
    time = Float64[]
    chr = T[]
    chz = T[]
    for ch in cer.channel
        if what == :n_imp
            ch_data = getproperty(ch.ion[1], :n_i_over_n_e)
        elseif what == :zeff
            ch_data = getproperty(ch, what)
        else
            ch_data = getproperty(ch.ion[1], what)
        end
        if time0 !== nothing
            selection = select_time_window(ch_data, :data, time0; window_size=time_averaging)
            if hasdata(ch_data, :data_σ)
                append!(data, Measurements.measurement.(selection.data, getproperty(ch_data, :data_σ)[selection.idx_range]))
            else
                append!(data, Measurements.measurement.(selection.data, selection.data .* 0.0))
            end
            append!(time, selection.time)
            append!(weights, selection.weights)
            append!(chr, ch.position.r.data[selection.idx_range])
            append!(chz, ch.position.z.data[selection.idx_range])
        else
            _data = getproperty(ch_data, :data)
            if hasdata(ch_data, :data_σ)
                _data = Measurements.measurement.(_data, getproperty(ch_data, :data_σ))
            else
                _data = Measurements.measurement.(_data, _data .* 0.0)
            end
            append!(data, _data)
            append!(time, getproperty(ch_data, :time))
            append!(chr, ch.position.r.data)
            append!(chz, ch.position.z.data)
        end
    end

    rho = chz .* 0.0
    for time0 in unique(time)
        i = nearest_causal_time(dd.equilibrium.time, time0; bounds_error=false).index
        eqt = dd.equilibrium.time_slice[i]
        r, z, RHO_interpolant = ρ_interpolant(eqt)
        index = time .== time0
        rho[index] = RHO_interpolant.(chr[index], chz[index])
    end

    if what == :n_imp
        cp1d = top_dd(cer).core_profiles.profiles_1d[time0]
        n_e = interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(rho)
        data = data .* n_e
    end

    return (time=time, rho=rho, data=data, weights=weights)
end

"""
    getdata(what_val::Union{Val{:t_e},Val{:n_e}}, dd::IMAS.dd{T}, time0::Union{Nothing,Float64}=nothing, time_averaging::Float64=0.0) where {T<:Real}

Extract Thomson scattering data for electron temperature or density.

  - time0: optional specific time for extraction

  - time_averaging: time averaging window size (required if time0 is specified)

Returns NamedTuple with (:time, :rho, :data, :weights)

Data is always of type Measurements.Measurement
"""
function getdata(what_val::Union{Val{:t_e},Val{:n_e}}, dd::IMAS.dd{T}, time0::Union{Nothing,Float64}=nothing, time_averaging::Float64=0.0) where {T<:Real}
    what = typeof(what_val).parameters[1]
    ts = dd.thomson_scattering

    if time0 !== nothing
        @assert time_averaging > 0.0
    end

    data = Measurements.Measurement[]
    weights = Float64[]
    time = Float64[]
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
            append!(time, selection.time)
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
            append!(time, getproperty(ch_data, :time))
            append!(chr, fill(ch.position.r, size(_data)))
            append!(chz, fill(ch.position.z, size(_data)))
        end
    end

    rho = chz .* 0.0
    for t in unique(time)
        i = nearest_causal_time(dd.equilibrium.time, t; bounds_error=false).index
        eqt = dd.equilibrium.time_slice[i]
        r, z, RHO_interpolant = ρ_interpolant(eqt)
        index = time .== t
        rho[index] = RHO_interpolant.(chr[index], chz[index])
    end

    return (time=time, rho=rho, data=data, weights=weights)
end

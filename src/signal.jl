using DSP

document[Symbol("Signal")] = Symbol[]

"""
    step(t::Float64)

Unitary step triggered at t=0
"""
function step(t::Float64)
    return Float64(t >= 0.0)
end

"""
    step(t::Float64, t_start::Float64)

Unitary step triggered at t=t_start
"""
function step(t::Float64, t_start::Float64)
    return step(t - t_start)
end

@compat public step
push!(document[Symbol("Signal")], :step)

"""
    pulse(t::Float64)

Unitary pulse with width of 1, starting at t=0
"""
function pulse(t::Float64)
    a = step(t)
    b = step(-t + 1)
    return a * b
end

"""
    pulse(t::Float64, t_start::Float64, Δt::Float64)

Unitary pulse with given width Δt, starting at t=t_start
"""
function pulse(t::Float64, t_start::Float64, Δt::Float64)
    return pulse((t - t_start) / Δt)
end

@compat public pulse
push!(document[Symbol("Signal")], :pulse)

"""
    ramp(t::Float64)

Unitary ramp from t=0 to t=1
"""
function ramp(t::Float64)
    a = t * (t < 1) * (t > 0)
    b = t >= 1
    return a + b
end

"""
    ramp(t::Float64, ramp_fraction::Float64)

Unitary ramp

The `ramp_fraction` defines the fraction of ramp with respect to 1.0 and must be between [0.0,1.0]

NOTE: This function is designed as is to be able to switch between `ramp(t, ramp_fraction)` and `trap(t, ramp_fraction)`.
"""
function ramp(t::Float64, ramp_fraction::Float64)
    @assert 0 <= ramp_fraction <= 1 "ramp `ramp_fraction` must be between [0.0,1.0]"
    k = 1.0 / ramp_fraction
    if ramp_fraction == 0.0
        return pulse(t)
    else
        return ramp(t * k)
    end
end

"""
    ramp(t::Float64, t_start::Float64, Δt::Float64)

Unitary ramp of duration Δt, starting at t=t_start
"""
function ramp(t::Float64, t_start::Float64, Δt::Float64)
    return ramp((t - t_start) / Δt)
end

@compat public ramp
push!(document[Symbol("Signal")], :ramp)

"""
    trap(t::Float64, ramp_fraction::Float64)

Unitary trapezoid

The `ramp_fraction` defines the fraction of ramp with respect to flattop and must be between [0.0,0.5]
"""
function trap(t::Float64, ramp_fraction::Float64)
    @assert 0 <= ramp_fraction <= 0.5 "trap `ramp_fraction` must be between [0.0,0.5]"
    k = 1.0 / ramp_fraction
    if ramp_fraction == 0.0
        return pulse(t)
    else
        a = ramp(t * k) * (t < 0.5)
        b = ramp(-t * k + k) * (t >= 0.5)
        return a + b
    end
end

"""
    trap(t::Float64, t_start::Float64, Δt::Float64, ramp_fraction::Float64)

Unitary trapezoid of duration Δt, starting at t=t_start

The `ramp_fraction` defines the fraction of ramp with respect to flattop and must be between [0.0,0.5]
"""
function trap(t::Float64, t_start::Float64, Δt::Float64, ramp_fraction::Float64)
    return trap((t - t_start) / Δt, ramp_fraction)
end

@compat public trap
push!(document[Symbol("Signal")], :trap)

"""
    gaus(t::Float64, order::Float64=1.0)

Unitary gaussian
"""
function gaus(t::Float64, order::Float64=1.0)
    return exp(-(t^2 / 2.0)^order)
end

"""
    gaus(t::Float64, t_start::Float64, Δt::Float64, order::Float64=1.0)

Unitary gaussian centered at t_start and with standard deviation Δt
"""
function gaus(t::Float64, t_start::Float64, Δt::Float64, order::Float64=1.0)
    return gaus((t - t_start) / Δt, order)
end

@compat public gaus
push!(document[Symbol("Signal")], :gaus)

"""
    beta(t::Float64, mode::Float64)

Unitary beta distribution

The `mode` [-1.0, 1.0] defines how skewed the distribution is
"""
function beta(t::Float64, mode::Float64)
    @assert -1 <= mode <= 1 "beta `mode` must be between [-1.0,1.0]"
    if t < 0 || t > 1
        return 0  # Outside the support
    end

    α = 2.0  # keeping alpha constant

    # mode expressed from 0 to 1
    mode = (mode / 2.0) + 0.5

    # Special conditions where the mode is at the boundaries
    if mode == 0.0
        return t == 0.0 ? 1.0 : 0.0  # spike at 0
    elseif mode == 1.0
        return t == 1.0 ? 1.0 : 0.0  # spike at 1
    end

    if mode > 0.5
        t = -t + 1.0
        mode = -(mode - 0.5) + 0.5
    end

    # For other cases, we find the corresponding β from the mode
    β = ((α - 1) / mode) - (α - 2)

    # Calculate the value of the beta distribution at its mode
    peak_value = (mode^(α - 1) * (1 - mode)^(β - 1))

    # normalize
    return (t^(α - 1) * (1 - t)^(β - 1)) / peak_value
end

"""
    beta(t::Float64, t_start::Float64, Δt::Float64, mode::Float64)

Unitary beta distribution of duration Δt, starting at t=t_start

The `mode` [-1.0, 1.0] defines how skewed the distribution is
"""
function beta(t::Float64, t_start::Float64, Δt::Float64, mode::Float64)
    return beta((t - t_start) / Δt, mode)
end

@compat public beta
push!(document[Symbol("Signal")], :beta)

"""
    sequence(t::Float64, t_y_sequence::Vector{Tuple{Float64,Float64}}; scheme::Symbol=:linear)

returns interpolated data given a sequence (tuple) of time/value points
"""
function sequence(t::Float64, t_y_sequence::Vector{Tuple{Float64,Float64}}; scheme::Symbol=:linear)
    tt = [t0 for (t0, y0) in t_y_sequence]
    yy = [y0 for (t0, y0) in t_y_sequence]
    return IMAS.extrap1d(IMAS.interp1d_itp(tt, yy, scheme); first=:flat, last=:flat)(t)
end

@compat public sequence
push!(document[Symbol("Signal")], :sequence)

"""
    moving_average(data::Vector{<:Real}, window_size::Int)

Calculate the moving average of a data vector using a specified window size.
The window size is always rounded up to the closest odd number to maintain symmetry around each data point.
"""
function moving_average(data::Vector{<:Real}, window_size::Int)
    smoothed_data = copy(data)
    window_size = isodd(window_size) ? window_size : window_size + 1
    pad_size = div(window_size - 1, 2)
    if pad_size < 1
        return smoothed_data
    end
    for i in eachindex(data)
        window_start = max(1, i - pad_size)
        window_end = min(length(data), i + pad_size)
        smoothed_data[i] = sum(data[window_start:window_end]) / (window_end - window_start + 1)
    end
    return smoothed_data
end

"""
    moving_average(t::AbstractVector, data::AbstractVector, new_time::AbstractVector; causal::Bool)

Calculate the moving average of a data vector on a new time basis
"""
function moving_average!(result::Ref{Float64}, w::Vector{Float64}, time::AbstractVector{Float64}, data::AbstractVector{T}, t0::Float64, width::Float64, interp, causal::Bool) where {T<:Real}
    if causal
        @. w = pulse(-time, -t0 - width, width)
    else
        @. w = pulse(-time, -t0 - width / 2, width)
    end
    norm = sum(w)
    if norm == 0.0
        result[] = interp(t0)
    else
        acc = 0.0
        for i in eachindex(w)
            acc += data[i] * (w[i] / norm)
        end
        result[] = acc
    end
end

function moving_average(time::AbstractVector{Float64}, data::AbstractVector{T}, new_time::AbstractVector{Float64}; causal::Bool=false) where {T<:Real}
    new_data = similar(new_time, Float64)
    width = new_time[2] - new_time[1]
    w = zeros(Float64, length(time))
    result = Ref{Float64}()
    interp = interp1d(time, data)

    for k in eachindex(new_time)
        moving_average!(result, w, time, data, new_time[k], width, interp, causal)
        new_data[k] = result[]
    end
    return new_data
end


@compat public moving_average
push!(document[Symbol("Signal")], :moving_average)

"""
    highpassfilter(signals, fs, cutoff, order=4)

Apply butterworth high pass filter of given order to signals sampled with frequency `fs`
"""
function highpassfilter(signals::AbstractArray{<:Real}, fs::Real, cutoff::Real, order=4)
    wdo = 2.0 * cutoff / fs
    filth = DSP.digitalfilter(DSP.Highpass(wdo), DSP.Butterworth(order))
    return DSP.filtfilt(filth, signals)
end


@compat public highpassfilter
push!(document[Symbol("Signal")], :highpassfilter)

"""
    lowpassfilter(signals, fs, cutoff, order=4)

Apply butterworth low pass filter of given order to signals sampled with frequency `fs`
"""
function lowpassfilter(signals::AbstractArray{<:Real}, fs::Real, cutoff::Real, order=4)
    wdo = 2.0 * cutoff / fs
    filth = DSP.digitalfilter(DSP.Lowpass(wdo), DSP.Butterworth(order))
    return DSP.filtfilt(filth, signals)
end

@compat public lowpassfilter
push!(document[Symbol("Signal")], :lowpassfilter)
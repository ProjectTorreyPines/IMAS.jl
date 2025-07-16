using DSP

document[Symbol("Signal")] = Symbol[]

"""
    heaviside(t::Float64)

Heaviside step triggered at t=0
"""
function heaviside(t::Float64)
    return Float64(t >= 0.0)
end

"""
    heaviside(t::Float64, t_start::Float64)

Heaviside step triggered at t=t_start
"""
function heaviside(t::Float64, t_start::Float64)
    return heaviside(t - t_start)
end

@compat public heaviside
push!(document[Symbol("Signal")], :heaviside)

"""
    pulse(t::Float64)

Unitary pulse with width of 1, starting at t=0
"""
function pulse(t::Float64)
    a = heaviside(t)
    b = heaviside(-t + 1)
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
    return IMAS.extrap1d(IMAS.interp1d_itp(tt, yy, scheme); first=:constant, last=:constant)(t)
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
function moving_average!(
    result::Ref{Float64},
    w::Vector{Float64},
    time::AbstractVector{Float64},
    data::AbstractVector{T},
    t0::Float64,
    width::Float64,
    interp,
    causal::Bool
) where {T<:Real}
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

# Window function implementations that work on the full weights array
# These now write directly to weights[idx] for each idx in idx_range
function window_weights_range!(weights, ::Val{:gaussian}, xi, idx_range, x_target, window_size)
    @inbounds for idx in idx_range
        x_rel = xi[idx] - x_target
        weights[idx] = exp(-0.5 * (4.0 * x_rel / window_size)^2) / (window_size / 4.0) / sqrt(2π)
    end
end

function window_weights_range!(weights, ::Val{:hanning}, xi, idx_range, x_target, window_size)
    a = window_size
    @inbounds for idx in idx_range
        x_rel = xi[idx] - x_target
        if abs(x_rel) <= a / 2.0
            weights[idx] = (1 - cos(2π * (x_rel + a / 2.0) / a)) / a
        else
            weights[idx] = 0.0
        end
    end
end

function window_weights_range!(weights, ::Val{:bartlett}, xi, idx_range, x_target, window_size)
    a = window_size
    @inbounds for idx in idx_range
        x_rel = xi[idx] - x_target
        if abs(x_rel) <= a / 2.0
            weights[idx] = 1 - abs(x_rel) * 2.0 / a
        else
            weights[idx] = 0.0
        end
    end
end

function window_weights_range!(weights, ::Val{:blackman}, xi, idx_range, x_target, window_size)
    a = window_size
    @inbounds for idx in idx_range
        x_rel = xi[idx] - x_target
        if abs(x_rel) <= a / 2.0
            weights[idx] = 0.42 - 0.5 * cos(2π * (x_rel + a / 2.0) / a) + 0.08 * cos(4π * (x_rel + a / 2.0) / a)
        else
            weights[idx] = 0.0
        end
    end
end

function window_weights_range!(weights, ::Val{:boxcar}, xi, idx_range, x_target, window_size)
    a = window_size
    @inbounds for idx in idx_range
        x_rel = xi[idx] - x_target
        if abs(x_rel) <= a / 2.0
            weights[idx] = 1.0 / a
        else
            weights[idx] = 0.0
        end
    end
end

function window_weights_range!(weights, ::Val{:triangle}, xi, idx_range, x_target, window_size)
    a = window_size
    @inbounds for idx in idx_range
        x_rel = xi[idx] - x_target
        if abs(x_rel) <= a
            weights[idx] = 1 - abs(x_rel) / a
        else
            weights[idx] = 0.0
        end
    end
end

"""
    select_by_window!(weights, xi, x_target, window_size, window_function; causal=false)

Fill weights array with window weights around x_target. Returns the index range of non-zero weights.
The weights array must be the same length as xi.

Returns:

  - idx_range: UnitRange of indices in xi that have non-zero weights
"""
function select_by_window!(
    weights::AbstractVector{Float64},
    xi::AbstractVector{T},
    x_target::Real,
    window_size::Real,
    window_function::Union{Symbol,Val};
    causal::Bool=false
) where {T<:Real}

    @assert issorted(xi)
    @assert length(weights) == length(xi)

    # Clear all weights first
    fill!(weights, 0.0)

    # Convert symbol to Val if needed
    wf_val = window_function isa Symbol ? Val(window_function) : window_function

    # Determine window extent (internal implementation detail)
    win_mult = if wf_val isa Val{:gaussian} || wf_val isa Val{:triangle}
        2.0
    else
        1.0
    end

    win_extent = window_size * win_mult
    half_extent = win_extent / 2.0

    # Find indices within window using binary search
    if causal
        idx_range = searchsortedfirst(xi, x_target):searchsortedlast(xi, x_target + win_extent)
    else
        idx_range = searchsortedfirst(xi, x_target - half_extent):searchsortedlast(xi, x_target + half_extent)
    end

    if isempty(idx_range)
        return 1:0
    end

    # Calculate weights in-place for the relevant range
    window_weights_range!(weights, wf_val, xi, idx_range, x_target, window_size)

    # Apply causal constraint if needed
    if causal
        @inbounds for idx in idx_range
            x_rel = xi[idx] - x_target
            if x_rel < 0
                weights[idx] = 0.0
            end
        end
    end

    # Normalize weights in the range
    norm = 0.0
    @inbounds for idx in idx_range
        norm += weights[idx]
    end

    if norm == 0
        @inbounds for idx in idx_range
            weights[idx] = NaN
        end
    else
        inv_norm = 1.0 / norm
        @inbounds for idx in idx_range
            weights[idx] *= inv_norm
        end
    end

    return idx_range
end

function select_time_window(ids::IDS, field::Symbol, time0::Float64; window_size::Float64=0.1, window_function::Symbol=:gaussian, causal::Bool=false)
    @assert time_coordinate_index(ids, field; error_if_not_time_dependent=true) == 1
    data = getproperty(ids, field)
    time = getproperty(time_coordinate(ids, field))
    weights = Vector{Float64}(undef, length(time))
    idx_range = select_by_window!(weights, time, time0, window_size, window_function; causal)
    return (time=time[idx_range], data=data[idx_range], weights=weights[idx_range], idx_range=idx_range)
end

"""
    smooth_by_convolution(
        yi::AbstractVector{T};
        xi=nothing,
        xo=nothing,
        window_size=nothing,
        window_function::Symbol=:gaussian,
        causal::Bool=false,
        interpolate::Int=0
    ) where {T<:Real}

Smooth non-uniformly sampled data with a chosen convolution window function. Supports error propagation.

The output values are nan where no points are found in finite windows (weight is zero).
The gaussian window is infinite in extent, and thus returns values for all xo.

:param yi: Values of input array

:param xi: Original grid points of input array (default y indicies)

:param xo: Output grid points of convolution array (default xi)

:param window_size: float.
Width of passed to window function (default maximum xi step).
For the Gaussian, sigma=window_size/4. and the convolution is integrated across +/-4.*sigma.

:param window_function: Symbol
Accepted strings are :hanning, :bartlett, :blackman, :gaussian, or :boxcar

:param axis: int. Axis of y along which convolution is performed

:param causal: int. Forces fw(x>0) = 0.

:param interpolate: Int > 0
Paramter indicating to interpolate data so that there are`interpolate`
number of data points within a time window. This is useful in presence of sparse
data, which would result in stair-case output if not interpolated.
The integer value sets the # of points per window size.

:return: convolved array on xo
"""
function smooth_by_convolution(
    yi::AbstractVector{T};
    xi=nothing,
    xo=nothing,
    window_size=nothing,
    window_function::Symbol=:gaussian,
    causal::Bool=false,
    interpolate::Int=0
) where {T<:Real}

    xi = xi === nothing ? collect(1:length(yi)) : xi
    xo = xo === nothing ? xi : xo

    @assert issorted(xi)

    if length(xi) < 2
        return zeros(T, length(xo)) .* NaN
    end

    # Remove any NaN
    if any(isnan, yi) || any(isnan, xi)
        index = .!isnan.(yi) .&& .!isnan.(xi)
        yi = yi[index]
        xi = xi[index]
    end

    # Convert symbol to Val for dispatch
    wf_val = Val(window_function)

    if window_size === nothing
        dx = diff(xi)  # xi is already sorted
        window_size = maximum(filter(x -> x > 0.0, dx))
    end

    # Handle interpolation
    if interpolate > 1 && length(xi) > 1
        n_interp = round(Int, (maximum(xi) - minimum(xi)) / window_size * interpolate, RoundUp)
        interp_y = DataInterpolations.LinearInterpolation(yi, xi)
        xi = range(minimum(xi); stop=maximum(xi), length=n_interp)
        yi = interp_y.(xi)
    end

    N = length(xo)
    yo = Vector{eltype(yi)}(undef, N)

    # Pre-allocate weights array once
    weights = Vector{Float64}(undef, length(xi))

    # Main smoothing loop using select_by_window!
    for k in 1:N
        idx_range = select_by_window!(weights, xi, xo[k], window_size, wf_val; causal=causal)

        if isempty(idx_range) || all(isnan, (@view weights[idx_range]))
            yo[k] = NaN
        else
            yo[k] = sum((@view weights[idx_range]) .* (@view yi[idx_range]))
        end
    end

    return yo
end

function smooth_by_convolution(ids::IDS, field::Symbol, time0::Vector{Float64}; window_size::Float64=0.1, window_function::Symbol=:gaussian, causal::Bool=false, interpolate::Int=0)
    @assert time_coordinate_index(ids, field; error_if_not_time_dependent=true) == 1
    data = getproperty(ids, field)
    field_σ = Symbol("$(field)_σ")
    if hasdata(ids, field_σ)
        data = Measurements.measurement.(data, getproperty(ids, field_σ))
    end
    time = getproperty(time_coordinate(ids, field))
    return smooth_by_convolution(data; xi=time, xo=time0, window_size, window_function, causal, interpolate)
end

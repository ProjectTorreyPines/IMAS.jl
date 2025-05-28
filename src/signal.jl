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

"""
    smooth_by_convolution(yi::AbstractVector; xi=nothing, xo=nothing, window_size=nothing, window_function="gaussian", causal=false, interpolate=false, std_dev=2)

Smooth non-uniformly sampled data with a chosen convolution window function. Supports error propagation.

The output values are nan where no points are found in finite windows (weight is zero).
The gaussian window is infinite in extent, and thus returns values for all xo.

:param yi: Values of input array

:param xi: Original grid points of input array (default y indicies)

:param xo: Output grid points of convolution array (default xi)

:param window_size: float.
    Width of passed to window function (default maximum xi step).
    For the Gaussian, sigma=window_size/4. and the convolution is integrated across +/-4.*sigma.

:param window_function: str/function.
    Accepted strings are 'hanning','bartlett','blackman','gaussian', or 'boxcar'.
    Function should accept x and window_size as arguments and return a corresponding weight.

:param axis: int. Axis of y along which convolution is performed

:param causal: int. Forces f(x>0) = 0.

:param interpolate: False or integer number > 0
    Paramter indicating to interpolate data so that there are`interpolate`
    number of data points within a time window. This is useful in presence of sparse
    data, which would result in stair-case output if not interpolated.
    The integer value sets the # of points per window size.

:param std_dev: str/int
    Accepted strings are 'none', 'propagate', 'population', 'expand', 'deviation', 'variance'.
    Only 'population' and 'none' are valid if yi is not an uncertainties array (i.e. std_devs(yi) is all zeros).
    Setting to an integer will convolve the error uncertainties to the std_dev power before taking the std_dev root.
    std_dev = 'propagate' is true propagation of errors (slow if not interpolating)
    std_dev = 'population' is the weighted "standard deviation" of the points themselves (strictly correct for the boxcar window)
    std_dev = 'expand' is propagation of errors weighted by w~1/window_function
    std_dev = 'deviation' is equivalent to std_dev=1
    std_dev = 'variance' is equivalent to std_dev=2

:return: convolved array on xo
"""
function smooth_by_convolution(yi::AbstractVector; xi=nothing, xo=nothing, window_size=nothing, window_function="gaussian", causal=false, interpolate=false, std_dev=2)

    xi = xi === nothing ? collect(1:length(yi)) : xi
    xo = xo === nothing ? xi : xo

    function get_window_function(name::String)
        if name == "gaussian"
            return (x, a) -> exp.(-0.5 * (4.0 * x / a).^2) ./ (a / 4.0) ./ sqrt(2π), 2.0
        elseif name == "hanning"
            return (x, a) -> (abs.(x) .<= a / 2.0) .* (1 .- cos.(2π * (x .- a / 2.0) ./ a)) ./ a, 1.0
        elseif name == "bartlett"
            return (x, a) -> (abs.(x) .<= a / 2.0) .* (1 .- abs.(x) * 2.0 / a), 1.0
        elseif name == "blackman"
            return (x, a) -> (abs.(x) .<= a / 2.0) .* (0.42 .- 0.5 * cos.(2π * (x .+ a / 2.0) ./ a) .+ 0.08 * cos.(4π * (x .+ a / 2.0) ./ a)), 1.0
        elseif name == "boxcar"
            return (x, a) -> (abs.(x) .<= a / 2.0) ./ a, 1.0
        elseif name == "triangle"
            return (x, a) -> (abs.(x) .<= a) .* (1 .- abs.(x) ./ a), 2.0
        else
            error("Unknown window function: $name")
        end
    end

    f, win_mult = typeof(window_function) == String ? get_window_function(window_function) : (window_function, 1.0)

    if window_size === nothing
        dx = diff(sort(xi))
        window_size = maximum(dx)
    end

    win_extent = window_size * win_mult
    half_extent = win_extent / 2.0

    if interpolate != false && length(xi) > 1
        interp_points = interpolate === true ? 10 : interpolate
        n_interp = ceil(Int, (maximum(xi) - minimum(xi)) / win_extent * interp_points)
        xi_new = range(minimum(xi), stop=maximum(xi), length=n_interp)
        interp_y = DataInterpolations.LinearInterpolation(xi, yi)
        yi = [interp_y(x) for x in xi_new]
        xi = xi_new
    end

    N = length(xo)
    yo = Vector{eltype(yi)}(undef, N)

    for k in 1:N
        xdiff = xo[k] .- xi
        mask = causal ? ((xdiff .>= 0) .& (xdiff .<= win_extent)) : (abs.(xdiff) .<= half_extent)

        if !any(mask)
            yo[k] = NaN
            continue
        end

        x_sel = xdiff[mask]
        y_sel = yi[mask]
        w = f(x_sel, window_size)
        if causal
            w .= w .* (x_sel .>= 0)
        end

        norm = sum(w)
        if norm == 0
            yo[k] = NaN
            continue
        end

        w ./= norm

        # === Handle std_dev cases ===
        if std_dev == "none"
            # Suppress uncertainty, treat values as plain reals
            yo[k] = sum(w .* y_sel)  # works for Measurement and Float64
        elseif std_dev == "population"
            μ = sum(w .* y_sel)
            w2 = sum(w .^ 2)
            denom = 1 - w2
            variance = denom > 0 ? sum(w .* (y_sel .- μ).^2) / denom : 0.0
            yo[k] = μ + sqrt(variance) * zero(eltype(yi))  # adds zero uncertainty if `yi` is Float64
        elseif std_dev == "propagate"
            yo[k] = sum(w .* y_sel)  # Measurement will propagate automatically
        elseif std_dev == "expand"
            max_w = maximum(w)
            # expand uncertainty (works if y_sel has it, otherwise does nothing)
            yo[k] = sum(w .* (max_w .* y_sel))
        elseif std_dev == "deviation" || std_dev == "variance" || isa(std_dev, Int)
            p = std_dev == "deviation" ? 1 :
                std_dev == "variance" ? 2 : std_dev
            err = (sum(abs.(w .* y_sel).^p)) ^ (1 / p)
            yo[k] = err  # not necessarily a measurement unless input was
        else
            error("Invalid std_dev: $std_dev")
        end
    end

    return yo
end

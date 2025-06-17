import NaturalNeighbours

"""
    fit2d(what::Val, dd::IMAS.dd{T}; transform::F=x -> x) where {T<:Real, F<:Function}

Create 2D interpolation of experimental data over time and radial coordinate.

Optional transformation function applied to data before interpolation

Returns interpolation function for 2D data (rho, time) -> transformed_data

NOTE: what is a Val{<:Symbol} that gets passed to the `getrawdata(what, dd)`
"""
function fit2d(what::Val, dd::IMAS.dd{T}; transform::F=x -> x) where {T<:Real,F<:Function}
    # get data
    time, rho, data = getdata(what, dd)

    @assert eltype(data) <: Measurements.Measurement

    # remove any NaN
    index = .!isnan.(rho) .&& .!isnan.(data)
    if sum(.!index) != length(data)
        time = @views time[index]
        rho = @views rho[index]
        data = @views data[index]
    end

    return NaturalNeighbours.interpolate(rho, time, transform.(data))
end

"""
    fit1d(rho::AbstractVector{T1}, data::AbstractVector{T2}, rho_tor_norm::AbstractVector{T3}; smooth1::Float64, smooth2::Float64) where {T1<:Real,T2<:Real,T3<:Real}

Fit 1D profile data with spatial smoothing using linearized coordinate transformation.

  - smooth1: smoothing parameter for space linearization

  - smooth2: smoothing parameter for linearized space

Returns NamedTuple with fitted data and coordinate information (:rho, :data, :fit, :rho_tor_norm, :rho_linearized)
"""
function fit1d(rho::AbstractVector{T1}, data::AbstractVector{T2}, rho_tor_norm::AbstractVector{T3}; smooth1::Float64, smooth2::Float64) where {T1<:Real,T2<:Real,T3<:Real}
    # linearize space
    result = smooth_by_convolution(Measurements.value.(data); xi=rho, xo=rho_tor_norm, window_size=smooth1)
    g = abs.(gradient(rho_tor_norm, result) ./ result)
    g .= cumtrapz(rho_tor_norm, g)
    rho_inverse = @. (g - g[1]) / (g[end] - g[1]) * (rho_tor_norm[end] - rho_tor_norm[1]) + rho_tor_norm[1]
    rho_linearized = interp1d(rho_tor_norm, rho_inverse).(rho)

    result = smooth_by_convolution(data; xi=rho_linearized, xo=rho_inverse, window_size=smooth2)

    return (rho=rho, data=data, fit=result, rho_tor_norm=rho_tor_norm, rho_linearized=rho_linearized)
end
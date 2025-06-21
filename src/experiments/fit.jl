import NaturalNeighbours

struct NaturalNeighboursMeasurementInterpolator
    itp
    itp_σ
end

function (nnmi::NaturalNeighboursMeasurementInterpolator)(rho::Real, time::Real)
    return Measurements.measurement(nnmi.itp(rho, time), nnmi.itp_σ(rho, time))
end

function (nnmi::NaturalNeighboursMeasurementInterpolator)(rho::AbstractVector{<:Real}, time::AbstractVector{<:Real})
    return [Measurements.measurement(v, e) for (v, e) in zip(nnmi.itp(rho, time), nnmi.itp_σ(rho, time))]
end

"""
    fit2d(what::Val, dd::IMAS.dd{T}; transform::F=x -> x) where {T<:Real, F<:Function}

Create 2D interpolation of experimental data over time and radial coordinate.

Optional transformation function applied to data before interpolation

Returns interpolation function for 2D data (rho, time) -> transformed_data

NOTE: what is a Val{<:Symbol} that gets passed to the `getrawdata(what, dd)`
"""
function fit2d(what::Val, dd::IMAS.dd{T}; transform::F=x -> x) where {T<:Real,F<:Function}
    # get data
    time, rho, data_measurement = getdata(what, dd)

    @assert eltype(data_measurement) <: Measurements.Measurement
    @assert eltype(rho) <: T
    @assert eltype(time) <: Float64

    data = [d.val for d in data_measurement]
    #data_σ = [d.err for d in data_measurement]

    # remove any NaN
    index = .!isnan.(rho) .&& .!isnan.(data)
    if sum(.!index) != length(data)
        time = @views time[index]
        rho = @views rho[index]
        data = @views data[index]
        #data_σ = @views data_σ[index]
    end

    itp = NaturalNeighbours.interpolate(rho, time, transform.(data))
    # itp_σ = NaturalNeighbours.interpolate(rho, time, transform.(data_σ))

    return itp #NaturalNeighboursMeasurementInterpolator(itp, itp_σ)
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
    g ./= maximum(g)
    g .= g .+ rho_tor_norm
    rho_inverse = @. (g - g[1]) / (g[end] - g[1]) * (rho_tor_norm[end] - rho_tor_norm[1]) + rho_tor_norm[1]
    rho_linearized = interp1d(rho_tor_norm, rho_inverse).(rho)

    result = smooth_by_convolution(data; xi=rho_linearized, xo=rho_inverse, window_size=smooth2)

    return (rho=rho, data=data, fit=result, rho_tor_norm=rho_tor_norm, rho_linearized=rho_linearized)
end
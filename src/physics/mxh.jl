import Interact

mutable struct MXH
    R0::Real          # Major Radius
    Z0::Real          # Elevation
    ϵ::Real           # Inverse aspect ratio a/R0
    κ::Real           # Elongation
    c0::Real          # Tilt
    c::Vector{<:Real} # Cosine coefficients acos.([ovality,...])
    s::Vector{<:Real} # Sine coefficients asin.([triangularity,-squareness,...]
end

function MXH(R0::Real, n_coeffs::Integer)
    return MXH(R0, 0.0, 0.3, 1.0, 0.0, zeros(n_coeffs), zeros(n_coeffs))
end

function flat_coeffs(mxh::MXH)
    return vcat(mxh.R0, mxh.Z0, mxh.ϵ, mxh.κ, mxh.c0, mxh.c, mxh.s)
end

function MXH(flat_coeffs::Vector{<:Real})
    R0, Z0, ϵ, κ, c0 = flat_coeffs[1:5]
    L = Int((length(flat_coeffs) - 5) / 2)
    c = flat_coeffs[6:5+L]
    s = flat_coeffs[7+L:5+L*2]
    return MXH(R0, Z0, ϵ, κ, c0, c, s)
end

"""
    MXH_moment(f, w, d)

This does Int[f.w]/Int[w.w]
If w is a pure Fourier mode, this gives the Fourier coefficient
"""
function MXH_moment(f, w, d)
    # Could probably be replaced by some Julia trapz
    s0 = sum((f[1:end-1] .* w[1:end-1] .+ f[2:end] .* w[2:end]) .* d[1:end-1])
    s1 = sum((w[1:end-1] .* w[1:end-1] .+ w[2:end] .* w[2:end]) .* d[1:end-1])
    res = s0 / s1
    return res
end

"""
    MXH(pr::Vector{T}, pz::Vector{T}, MXH_modes::Integer) where {T<:Real}

Compute Fourier coefficients for Miller-extended-harmonic representation:

    R(r,θ) = R(r) + a(r)*cos(θᵣ(r,θ)) where θᵣ(r,θ) = θ + C₀(r) + sum[Cᵢ(r)*cos(i*θ) + Sᵢ(r)*sin(i*θ)]
    Z(r,θ) = Z(r) + κ(r)*a(r)*sin(θ)

Where pr,pz are the flux surface coordinates and MXH_modes is the number of modes
"""
function MXH(pr::Vector{T}, pz::Vector{T}, MXH_modes::Integer=5) where {T<:Real}
    R0 = 0.5 * (maximum(pr) + minimum(pr))
    Z0 = 0.5 * (maximum(pz) + minimum(pz))
    a = 0.5 * (maximum(pr) - minimum(pr))
    b = 0.5 * (maximum(pz) - minimum(pz))
    return MXH(pr, pz, R0, Z0, a, b, MXH_modes)
end

using Plots

function MXH(pr::Vector{T}, pz::Vector{T}, R0::T, Z0::T, a::T, b::T, MXH_modes::Integer) where {T<:Real}
    L = length(pr)
    θ = zeros(L)
    θᵣ = zeros(L)

    # Calculate angles with proper branches
    th = 0.0
    thr = 0.0
    for j in 1:L
        th_old = th
        thr_old = thr
        aa = (pz[j] - Z0) / b
        aa = max(-1, min(1, aa))
        th = asin(aa)
        bb = (pr[j] - R0) / a
        bb = max(-1, min(1, bb))
        thr = acos(bb)
        if (j == 1) || ((th > th_old) && (thr > thr_old))
            θ[j] = th
            θᵣ[j] = thr
        elseif (th < th_old) && (thr > thr_old)
            θ[j] = π - th
            θᵣ[j] = thr
        elseif (th < th_old) && (thr < thr_old)
            θ[j] = π - th
            θᵣ[j] = 2π - thr
        elseif (th > th_old) && (thr < thr_old)
            θ[j] = 2π + th
            θᵣ[j] = 2π - thr
        end
    end

    dθ = zeros(L)
    for j in 1:L-1
        dθ[j] = θ[j+1] - θ[j]
        if dθ[j] < 0
            dθ[j] += 2π
        end
    end
    dθ[end] = dθ[1]

    tilt = MXH_moment(θᵣ - θ, ones(L), dθ)

    sin_coeffs = zeros(MXH_modes)
    cos_coeffs = zeros(MXH_modes)
    for m in 1:MXH_modes
        sin_coeffs[m] = MXH_moment(θᵣ - θ, sin.(m * θ), dθ)
        cos_coeffs[m] = MXH_moment(θᵣ - θ, cos.(m * θ), dθ)
    end

    return MXH(R0, Z0, a / R0, b / a, tilt, cos_coeffs, sin_coeffs)
end

function (mxh::MXH)(adaptive_grid_N::Integer=100)
    step = mxh.R0 / adaptive_grid_N
    a = mxh.ϵ * mxh.R0
    N = Int(ceil(2π * a * mxh.κ / step))
    Θ = LinRange(0, 2π, N)
    return mxh.(Θ)
end

function (mxh::MXH)(θ::Real)
    a = mxh.ϵ * mxh.R0
    if length(mxh.c) > 0
        c_sum = sum(mxh.c[n] * cos(n * θ) for n in 1:length(mxh.c))
    else
        c_sum = 0.0
    end
    if length(mxh.s) > 0
        s_sum = sum(mxh.s[n] * sin(n * θ) for n in 1:length(mxh.s))
    else
        s_sum = 0.0
    end
    θ_R = θ + mxh.c0 + c_sum + s_sum
    r = mxh.R0 + a * cos(θ_R)
    z = mxh.Z0 + mxh.κ * a * sin(θ)
    return r, z
end

@recipe function plot_mxh(mxh::MXH; adaptive_grid_N=100)
    @series begin
        aspect_ratio --> :equal
        label --> ""
        mxh(adaptive_grid_N)
    end
end

function Base.show(io::IO, mxh::MXH)
    println(io, "R0: $(mxh.R0)")
    println(io, "Z0: $(mxh.Z0)")
    println(io, "ϵ: $(mxh.ϵ)")
    println(io, "κ: $(mxh.κ)")
    println(io, "c0: $(mxh.c0)")
    println(io, "c: $(mxh.c)")
    println(io, "s: $(mxh.s)")
end

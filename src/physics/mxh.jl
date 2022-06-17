"""
    MXH_moment(f, w, d)

This does Int[f.w]/Int[w.w]
If w is a pure Fourier mode, this gives the Fourier coefficient
"""
function MXH_moment(f, w, d)
    # Could probably be replaced by some Julia trapz
    s0 = sum((f[1:end-1] .* w[1:end-1] .+ f[2:end] .* w[2:end]) .* d[1:end-1])
    s1 = sum((w[1:end-1] .* w[1:end-1] .+ w[2:end] .* w[2:end]) .* d[1:end-1])
    return s0 / s1
end

"""
    miller_extended_harmonic(pr::Vector{T}, pz::Vector{T}, Rm::T, Zm::T, a::T, b::T, MXH_modes::Integer) where {T<:Real}

Compute Fourier coefficients for Miller-extended-harmonic representation:

    R(r,θ) = R(r) + a(r)*cos(θᵣ(r,θ)) where θᵣ(r,θ) = θ + C₀(r) + sum[Cᵢ(r)*cos(i*θ) + Sᵢ(r)*sin(i*θ)]
    Z(r,θ) = Z(r) + κ(r)*a(r)*sin(θ)

Where pr,pz are the flux surface coordinates around the Rm,Zm magnetic axis.
a and b are the minor radius and the major radius of the flux surface ellipse,
and MXH_modes is the number of modes
"""
function miller_extended_harmonic(pr::Vector{T}, pz::Vector{T}, Rm::T, Zm::T, a::T, b::T, MXH_modes::Integer) where {T<:Real}
    MXH_amplitudes = zeros(2 * MXH_modes + 1)
    th = 0.0
    th_old = 0.0
    thr = 0.0
    thr_old = 0.0

    L = length(pr)
    θ = zeros(L)
    θᵣ = zeros(L)

    # Calculate angles with proper branches
    for j in 1:L
        th_old = th
        thr_old = thr
        aa = (pz[j] - Zm) / b
        th = asin(min(-1, max(1, aa)))
        bb = (pr[j] - Rm) / a
        thr = acos(min(-1, max(1, bb)))
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
        dθ[j] < 0 && (dθ[j] += 2π)
    end
    dθ[end] = dθ[1]

    MXH_amplitudes[1] = MXH_moment(θᵣ - θ, ones(L), dθ)
    for m in 1:MXH_modes
        MXH_amplitudes[2m] = MXH_moment(θᵣ - θ, sin.(m * θ), dθ)
        MXH_amplitudes[2m+1] = MXH_moment(θᵣ - θ, cos.(m * θ), dθ)
    end
    return MXH_amplitudes
end

function Θr(θ, x)
    θr = θ + x[5]
    M = (length(x) - 5) ÷ 2
    for m in 1:M
        S, C = sincos(m * θ)
        θr += x[4+2m] * S + x[5+2m] * C
    end
    return θr
end

function R_MXH(θ, x)
    return x[1] + x[3] * cos(Θr(θ, x))
end

function Z_MXH(θ, x)
    return x[2] - x[3] * x[4] * sin(θ)
end
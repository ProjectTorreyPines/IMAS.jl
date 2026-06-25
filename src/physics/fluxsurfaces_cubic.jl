# Cubic-interpolant-based flux-surface geometry.
#
# Flux-surface geometry computed directly from the ψ cubic interpolant (analytic
# gradient/Hessian via FastInterpolations, aliased `FI` in IMAS.jl) with Newton
# iteration — as opposed to the Contour-based tracing in fluxsurfaces.jl.
#
# Currently: X-point-aware refinement of a traced surface's geometric extrema
# (max_r/min_r/max_z/min_z). Intended to also host a general cubic-interpolant-
# based flux-surface tracer. All helpers here are private (no public API yet).

# Split form of the critical-point system ∇ψ = 0 (Jacobian = Hessian) for the globalized
# [`_damped_newton2d`](@ref): `residual` returns the residual (∇ψ) and the first Hessian row,
# `hessrow` the second — same contract as [`_extremum_eqs`](@ref), so the one safe solver serves
# both systems. (Critical-point finding is rare, so the probe recomputing the Hessian is not hot.)
function _critical_eqs(H::AbstractMatrix, itp)
    residual(R, Z) = begin
        g = _gradient(itp, R, Z)
        _hessian!(H, itp, R, Z)
        return (g[1], g[2], H[1, 1], H[1, 2])
    end
    hessrow(R, Z) = begin
        _hessian!(H, itp, R, Z)
        return (H[2, 1], H[2, 2])
    end
    return (; residual, hessrow)
end

# Split form of the extremum system for the globalized [`_damped_newton2d`](@ref): `residual`
# (one `value_gradient`) returns the residual and the ψ=target Jacobian row (∇ψ, free from the
# same eval); `hessrow` (one `hessian!`) returns the ∂ψ/∂n Jacobian row. The line search probes
# with `residual` alone and pays `hessrow` only on accepted steps — halving the interpolant work
# per Newton iterate, with a bit-identical iterate sequence.
function _extremum_eqs(H::AbstractMatrix, itp, target_psi::Real, daxis::Int)
    residual(R, Z) = begin
        val, g = _value_gradient(itp, R, Z)
        return (val - target_psi, g[daxis], g[1], g[2])
    end
    hessrow(R, Z) = begin
        _hessian!(H, itp, R, Z)
        return (H[daxis, 1], H[daxis, 2])
    end
    return (; residual, hessrow)
end

"""
    _damped_newton2d(eqs, seed; tol=1e-11, maxit=50, αmin=1e-3)

Globalized 2×2 Newton with a backtracking line search on `‖F‖` and a gradient-descent fallback
when the Jacobian is near-singular — robust where pure Newton diverges (e.g. near the separatrix,
where `det(J) = ψ_R·ψ_ZZ → 0`, small poloidal curvature, blows up the plain step). Returns
`((R, Z), converged::Bool)`.

`eqs` is the split extremum system from [`_extremum_eqs`](@ref): `eqs.residual(R,Z)` returns
`(F1, F2, J11, J12)` — the residual plus the ψ=target Jacobian row (∇ψ), both from one
`value_gradient` — and `eqs.hessrow(R,Z)` returns `(J21, J22)`, the `∂ψ/∂n` row (a `hessian`).
The line search probes with `residual` alone; `hessrow` runs only on accepted steps, halving the
interpolant work per iterate versus a combined residual+Jacobian eval (the iterate sequence is
bit-identical: the `‖F‖` accept test uses only the residual, and `hessrow` runs at the exact
accepted point a combined eval would have).
"""
function _damped_newton2d(eqs, seed::Tuple{T,T}; tol::Real=1e-11, maxit::Int=50, αmin::Real=1e-3,
    lo=(-Inf, -Inf), hi=(Inf, Inf)) where {T<:Real}
    # every iterate is clamped to the box [lo, hi] (the ψ grid domain) so the search can never
    # wander into the extrapolation region outside the grid — the iterate physically cannot escape.
    R, Z = clamp(seed[1], lo[1], hi[1]), clamp(seed[2], lo[2], hi[2])
    F1, F2, J11, J12 = eqs.residual(R, Z)
    J21, J22 = eqs.hessrow(R, Z)
    nrm = hypot(F1, F2)
    for _ in 1:maxit
        nrm <= tol && return ((R, Z), true)
        det = J11 * J22 - J12 * J21
        if isfinite(det) && abs(det) > 1e-13
            dR = (J22 * F1 - J12 * F2) / det
            dZ = (J11 * F2 - J21 * F1) / det
        else  # near-singular J: descend 0.5‖F‖² along Jᵀ F with a small normalized step
            dR = J11 * F1 + J21 * F2
            dZ = J12 * F1 + J22 * F2
            s = hypot(dR, dZ)
            s > 0 && (dR /= s; dZ /= s)
            dR *= T(1e-2); dZ *= T(1e-2)
        end
        α = one(T); stepped = false
        while α >= αmin
            q1, q2 = clamp(R - α * dR, lo[1], hi[1]), clamp(Z - α * dZ, lo[2], hi[2])
            if isfinite(q1) && isfinite(q2)
                g1, g2, j11, j12 = eqs.residual(q1, q2)   # value_gradient only — no Hessian on probes
                if hypot(g1, g2) < nrm
                    R, Z = q1, q2
                    F1, F2, J11, J12 = g1, g2, j11, j12    # reuse the probe's value_gradient
                    J21, J22 = eqs.hessrow(R, Z)           # one Hessian on the accepted step
                    nrm = hypot(F1, F2)
                    stepped = true
                    break
                end
            end
            α /= 2
        end
        stepped || return ((R, Z), nrm <= tol)
    end
    return ((R, Z), nrm <= tol)
end

"""
    _robust_refine_extremum!(H, itp, target_psi, seed, extremum_of, axis; tol=1e-11, maxit=50)

X-point-aware refinement of a flux-surface geometric extremum. Solves the extremum system
`{ψ = target_psi, ∂ψ/∂n = 0}` from `seed` with a globalized [`_damped_newton2d`](@ref)
(no divergence near the separatrix), then validates the result by *physical region* — cheaply,
using the precomputed critical points (`axis` = O-point, plus `xpoints` = X-points) as
references, with NO per-call critical-point search:

  1. **axis-relative direction** — the extremized coordinate must be on the same side of the
     magnetic `axis` as the seed (kills the opposite O-point solution);
  2. **confined side of every X-point** — `(p−xp)·(axis−xp) > 0` for each `xp ∈ xpoints`
     (kills the private-flux-region / across-X-point solution). On-surface already excludes the
     SOL (ψ_N>1), so the only ambiguity left is confined-surface vs private-region, which an
     X-point dot test settles. Generic: `xpoints` empty (limited plasma → no private region)
     reduces to the direction check; one or two X-points (single/double null) just loop.

`xpoints` is a collection of `(R, Z)` (e.g. from `eqt.boundary.x_point`). If the solution is
on-surface and confined it is accepted; otherwise (a private/wrong branch, only reachable from
a bad seed) it re-seeds along the segment toward the axis and re-solves; failing that, falls
back to `seed` (always an on-surface contour vertex).
"""
function _robust_refine_extremum!(H::AbstractMatrix, itp, target_psi::T, seed::Tuple{T,T},
    extremum_of::Symbol, axis::Tuple{T,T}, xpoints=(); tol::Real=1e-11, maxit::Int=50, lo=(-Inf, -Inf), hi=(Inf, Inf)) where {T<:Real}
    d = extremum_of === :R ? 2 : extremum_of === :Z ? 1 :
        throw(ArgumentError("_robust_refine_extremum!: extremum_of must be :R or :Z, got :$extremum_of"))
    e = extremum_of === :R ? 1 : 2
    want_max = seed[e] > axis[e]
    onsurf(q) = abs(itp(q[1], q[2]) - target_psi) < 1e-7
    function confined(q::Tuple{T,T})
        ((q[e] > axis[e]) == want_max) || return false                       # (1) axis-relative direction
        for xp in xpoints                                                    # (2) confined side of every X-point
            (q[1] - xp[1]) * (axis[1] - xp[1]) + (q[2] - xp[2]) * (axis[2] - xp[2]) > 0 || return false
        end
        return true
    end

    eqs = _extremum_eqs(H, itp, target_psi, d)
    p, ok = _damped_newton2d(eqs, seed; tol, maxit, lo, hi)
    (ok && onsurf(p) && confined(p)) && return p

    # wrong branch (private/across-X-point) — only reachable from a bad seed. Re-seed along the
    # segment from `p` toward the magnetic `axis` (always confined) and re-solve: the axis anchor
    # pulls the Newton back onto the confined branch.
    for f in (T(0.3), T(0.5), T(0.7))
        reseed = (p[1] + f * (axis[1] - p[1]), p[2] + f * (axis[2] - p[2]))
        q, okq = _damped_newton2d(eqs, reseed; tol, maxit, lo, hi)
        (okq && onsurf(q) && confined(q)) && return q
    end
    return seed
end

# convenience wrapper for the robust refine (allocates the 2×2 scratch; one-off calls / tests)
_robust_refine_extremum(itp, target_psi::T, seed::Tuple{T,T}, extremum_of::Symbol,
    axis::Tuple{T,T}, xpoints=(); kw...) where {T<:Real} =
    _robust_refine_extremum!(Matrix{T}(undef, 2, 2), itp, target_psi, seed, extremum_of, axis, xpoints; kw...)

"""
    _project_to_level(itp::FI.AbstractInterpolant, c::T, x::Tuple{T,T}; tol::Real=1e-10, maxit::Int=8) where {T<:Real}

Newton corrector that projects `x = (R,Z)` onto the level set ψ=c **along the gradient**:
`x ← x - (ψ(x) - c) · ∇ψ / |∇ψ|²` (the codimension-1 Moore–Penrose Newton step; quadratic
convergence). Returns `((R, Z), converged::Bool)`; `converged=false` if `|∇ψ|` collapses
(a critical point) or iterates go non-finite.
"""
function _project_to_level(itp::FI.AbstractInterpolant, c::T, x::Tuple{T,T}; tol::Real=1e-10, maxit::Int=8) where {T<:Real}
    R, Z = x
    for _ in 1:maxit
        val, g = FI.value_gradient(itp, (R, Z))
        r = val - c
        abs(r) <= tol && return ((R, Z), true)
        n2 = g[1]^2 + g[2]^2
        (isfinite(n2) && n2 > eps(T)) || return ((R, Z), false)   # critical point / degenerate
        R -= r * g[1] / n2
        Z -= r * g[2] / n2
        (isfinite(R) && isfinite(Z)) || return ((R, Z), false)
    end
    val = itp(R, Z)
    return ((R, Z), abs(val - c) <= tol)
end

"""
    _contour_tangent(itp::FI.AbstractInterpolant, x::Tuple{T,T}, sgn::Int) where {T<:Real}

Unit tangent of the ψ iso-contour at `x`: `sgn·(ψ_Z, -ψ_R)/|∇ψ|` (rotate ∇ψ by 90°).
`sgn ∈ {+1,-1}` selects traversal orientation (CW vs CCW).
"""
function _contour_tangent(itp::FI.AbstractInterpolant, x::Tuple{T,T}, sgn::Int) where {T<:Real}
    g = FI.gradient(itp, x)
    nrm = hypot(g[1], g[2])
    return (sgn * g[2] / nrm, -sgn * g[1] / nrm)
end

"""
    _contour_curvature(itp::FI.AbstractInterpolant, x::Tuple{T,T}) where {T<:Real}

Signed curvature of the ψ iso-contour at `x`:
`κ = (ψ_R²ψ_ZZ - 2ψ_Rψ_Zψ_RZ + ψ_Z²ψ_RR) / |∇ψ|³` (needs the Hessian).

`H` is an optional caller-owned 2×2 scratch buffer; when threaded through the tracer's step
loop it lets the per-step `hessian!` avoid allocating a fresh matrix each call.
"""
function _contour_curvature(itp::FI.AbstractInterpolant, x::Tuple{T,T},
    H::AbstractMatrix=Matrix{T}(undef, 2, 2)) where {T<:Real}
    g = FI.gradient(itp, x)
    FI.hessian!(H, itp, x)
    gR, gZ = g[1], g[2]
    n2 = gR^2 + gZ^2
    return (gR^2 * H[2, 2] - 2 * gR * gZ * H[1, 2] + gZ^2 * H[1, 1]) / n2^(T(3) / 2)
end

"""
    _step_pc(itp::FI.AbstractInterpolant, c::T, x::Tuple{T,T}, h::Real, sgn::Int;
             corr_tol::Real=1e-10, corr_max::Int=8) where {T<:Real}

One predictor–corrector step of arclength `h` along the ψ=c contour. Predictor is the
Hessian-osculating 2nd-order Taylor step `x + h·t + ½h²·κ·n_p` (`t` tangent, `n_p` principal
normal `(-t_Z, t_R)`, `κ` curvature); corrector is [`_project_to_level`](@ref) back onto ψ=c.
`H` is an optional caller-owned 2×2 Hessian scratch buffer forwarded to [`_contour_curvature`](@ref).
Returns `(xnew, on_surface::Bool)`.
"""
function _step_pc(itp::FI.AbstractInterpolant, c::T, x::Tuple{T,T}, h::Real, sgn::Int,
    H::AbstractMatrix=Matrix{T}(undef, 2, 2); corr_tol::Real=1e-10, corr_max::Int=8) where {T<:Real}
    t = _contour_tangent(itp, x, sgn)
    np = (-t[2], t[1])                            # principal normal (rotate tangent +90°)
    κ = _contour_curvature(itp, x, H)
    xp = (x[1] + h * t[1] + (h^2 / 2) * κ * np[1],
          x[2] + h * t[2] + (h^2 / 2) * κ * np[2])
    return _project_to_level(itp, c, xp; tol=corr_tol, maxit=corr_max)
end

"""
    _contour_step(itp::FI.AbstractInterpolant, x::Tuple{T,T}; ε::Real=1e-6,
                  h_min::Real, h_max::Real, max_turn::Real=deg2rad(10), κ_floor::Real=1e-8) where {T<:Real}

Arclength step for the contour at `x`: bound the predictor chord error `≈ ½|κ|h² ≤ ε` and
the per-step turning angle `|κ|·h ≤ max_turn`, then clamp to `[h_min, h_max]`. `κ_floor`
prevents division by zero on near-straight pieces. `H` is an optional caller-owned 2×2 Hessian
scratch buffer forwarded to [`_contour_curvature`](@ref).
"""
function _contour_step(itp::FI.AbstractInterpolant, x::Tuple{T,T},
    H::AbstractMatrix=Matrix{T}(undef, 2, 2); ε::Real=1e-6,
    h_min::Real, h_max::Real, max_turn::Real=deg2rad(10), κ_floor::Real=1e-8) where {T<:Real}
    κ = max(abs(_contour_curvature(itp, x, H)), κ_floor)
    h = sqrt(2 * ε / κ)            # chord-error bound
    h = min(h, max_turn / κ)       # turning-angle cap
    return clamp(h, h_min, h_max)
end

# signed angle from vector a to vector b (radians, in (-π, π])
_signed_angle(a::Tuple, b::Tuple) = atan(a[1]*b[2] - a[2]*b[1], a[1]*b[1] + a[2]*b[2])

# does segment p->q cross the ray from x0 along normal m (the Poincaré section ⟂ t0)?
# returns the interpolation fraction in [0,1] if it crosses on the +section side, else nothing
function _section_cross(p::Tuple{T,T}, q::Tuple{T,T}, x0::Tuple{T,T}, t0::Tuple{T,T}) where {T<:Real}
    # Poincaré section through x0, ⟂ t0 (design §5.5). Detect sign change of the tangential
    # coordinate (point - x0)·t0; the section line's normal is t0 itself.
    fp = (p[1]-x0[1])*t0[1] + (p[2]-x0[2])*t0[2]
    fq = (q[1]-x0[1])*t0[1] + (q[2]-x0[2])*t0[2]
    (fp == fq) && return nothing
    (sign(fp) == sign(fq)) && return nothing       # no crossing of the section line
    frac = fp / (fp - fq)                           # in [0,1]
    # require the crossing to be near x0 along the section line (n0 ⟂ t0), not the far side
    cx = p[1] + frac*(q[1]-p[1]); cz = p[2] + frac*(q[2]-p[2])
    n0 = (-t0[2], t0[1])
    along = (cx-x0[1])*n0[1] + (cz-x0[2])*n0[2]
    return abs(along) < 1e-2 ? frac : nothing
end

"""
    _trace_surface_cubic(itp, c, seed; method=:pc, sgn=-1, ε, h_min, h_max, max_turn,
                         κ_floor, s_min, max_steps, domain) -> (Rs, Zs, closed)

Trace the ψ=c contour from `seed` by predictor–corrector continuation. Returns ordered
`(Rs, Zs)` and whether the loop closed. Closure is tested only after the accumulated turning
angle exceeds ~2π (a confined surface winds once), via a Poincaré section ⟂ the start
tangent. `domain=(Rlo,Rhi,Zlo,Zhi)` (or `nothing`) terminates open contours at the boundary.
"""
function _trace_surface_cubic(itp::FI.AbstractInterpolant, c::T, seed::Tuple{T,T};
    method::Symbol=:pc, sgn::Int=-1, ε::Real=1e-6, h_min::Real=1e-4, h_max::Real=0.1,
    max_turn::Real=deg2rad(10), κ_floor::Real=1e-8, s_min::Real=0.0, max_steps::Int=100_000,
    rk4_tol::Real=1e-8, τ_grad::Real=0.0,
    domain=nothing) where {T<:Real}

    x0, ok = _project_to_level(itp, c, seed)
    ok || return (T[seed[1]], T[seed[2]], false)
    t0 = _contour_tangent(itp, x0, sgn)
    Rs = T[x0[1]]; Zs = T[x0[2]]
    x = x0; tprev = t0; Θ = zero(T); s = zero(T)
    hrk = h_max
    H = Matrix{T}(undef, 2, 2)   # 2×2 Hessian scratch reused by every step (in-place hessian!)

    for _ in 1:max_steps
        # near-critical-point guard: when |∇ψ| < τ_grad (close to an X-point saddle), creep
        # at the smallest step via the corrector so the surface traces the confined side
        # instead of leaking through the saddle (τ_grad=0 disables this, default behavior).
        gg = τ_grad > 0 ? FI.gradient(itp, x) : nothing
        if gg !== nothing && hypot(gg[1], gg[2]) < τ_grad
            xnew, sok = _step_pc(itp, c, x, h_min, sgn, H)
        elseif method === :rk4
            xnew, hrk, sok = _step_rk4_adaptive(itp, x, hrk, sgn; tol=rk4_tol, h_min, h_max)
        else
            h = _contour_step(itp, x, H; ε, h_min, h_max, max_turn, κ_floor)
            xnew, sok = _step_pc(itp, c, x, h, sgn, H)
        end
        sok || break
        if domain !== nothing && !(domain[1] <= xnew[1] <= domain[2] && domain[3] <= xnew[2] <= domain[4])
            push!(Rs, xnew[1]); push!(Zs, xnew[2])
            return (Rs, Zs, false)                # open contour: exited domain
        end
        tnew = _contour_tangent(itp, xnew, sgn)
        Θ += _signed_angle(tprev, tnew)
        s += hypot(xnew[1]-x[1], xnew[2]-x[2])
        # closure: only after a full turn, and segment crosses the start section
        if abs(Θ) >= 2π - 0.5 && s > s_min
            frac = _section_cross(x, xnew, x0, t0)
            if frac !== nothing
                cx = x[1] + frac*(xnew[1]-x[1]); cz = x[2] + frac*(xnew[2]-x[2])
                xc, _ = _project_to_level(itp, c, (cx, cz))   # interp crossing is O(h²) off the curve; snap onto ψ=c
                push!(Rs, xc[1]); push!(Zs, xc[2])
                return (Rs, Zs, true)
            end
        end
        push!(Rs, xnew[1]); push!(Zs, xnew[2])
        x = xnew; tprev = tnew
    end
    return (Rs, Zs, false)
end

"""
    _resample_contour(Rs::AbstractVector{T}, Zs::AbstractVector{T}, n::Int) where {T<:Real}

Resample a traced polyline to `n` points equally spaced in cumulative arclength (linear
interpolation between traced vertices). Input is treated as ordered; the last point is the
closure point for closed surfaces.
"""
function _resample_contour(Rs::AbstractVector{T}, Zs::AbstractVector{T}, n::Int) where {T<:Real}
    m = length(Rs)
    if m == 1
        return (fill(Rs[1], n), fill(Zs[1], n))   # degenerate / single-point input
    end
    ll = zeros(T, m)
    for k in 2:m
        ll[k] = ll[k-1] + hypot(Rs[k]-Rs[k-1], Zs[k]-Zs[k-1])
    end
    L = ll[end]
    L > 0 || return (fill(Rs[1], n), fill(Zs[1], n))
    # detect a closed polyline (last vertex coincides with first); if so, sample the
    # half-open interval [0, L) so the n outputs are distinct and the wrap segment equals
    # the interior spacing instead of collapsing to ~0.
    closed = hypot(Rs[end]-Rs[1], Zs[end]-Zs[1]) < 1e-9 * max(L, one(T))
    targets = closed ? range(zero(T), L; length=n+1)[1:n] : range(zero(T), L; length=n)
    Ro = Vector{T}(undef, n); Zo = Vector{T}(undef, n)
    j = 1
    for (i, s) in enumerate(targets)
        while j < m && ll[j+1] < s
            j += 1
        end
        seg = ll[j+1] - ll[j]
        f = seg > 0 ? (s - ll[j]) / seg : zero(T)
        Ro[i] = Rs[j] + f * (Rs[j+1] - Rs[j])
        Zo[i] = Zs[j] + f * (Zs[j+1] - Zs[j])
    end
    return Ro, Zo
end

# one classic RK4 step of size h on dx/ds = tangent(x, sgn)
function _rk4(itp::FI.AbstractInterpolant, x::Tuple{T,T}, h::Real, sgn::Int) where {T<:Real}
    k1 = _contour_tangent(itp, x, sgn)
    x2 = (x[1] + h/2*k1[1], x[2] + h/2*k1[2]);  k2 = _contour_tangent(itp, x2, sgn)
    x3 = (x[1] + h/2*k2[1], x[2] + h/2*k2[2]);  k3 = _contour_tangent(itp, x3, sgn)
    x4 = (x[1] + h*k3[1],   x[2] + h*k3[2]);    k4 = _contour_tangent(itp, x4, sgn)
    return (x[1] + h/6*(k1[1]+2k2[1]+2k3[1]+k4[1]), x[2] + h/6*(k1[2]+2k2[2]+2k3[2]+k4[2]))
end

"""
    _step_rk4_adaptive(itp, x, h, sgn; tol=1e-8, h_min=1e-5, h_max=0.1)

Adaptive RK4 (step-doubling) on `dx/ds = tangent`. Compares one step of `h` against two of
`h/2`; accepts the (extrapolated) half-step if the estimated error ≤ `tol`, and returns the
next step size. **Comparison baseline only — no corrector**, so it measures the intrinsic
off-surface drift of pure integration. Returns `(xnew, h_next, accepted::Bool)`.
"""
function _step_rk4_adaptive(itp::FI.AbstractInterpolant, x::Tuple{T,T}, h::Real, sgn::Int;
    tol::Real=1e-8, h_min::Real=1e-5, h_max::Real=0.1) where {T<:Real}
    hh = clamp(h, h_min, h_max)
    for _ in 1:20
        big = _rk4(itp, x, hh, sgn)
        half = _rk4(itp, _rk4(itp, x, hh/2, sgn), hh/2, sgn)
        err = hypot(big[1]-half[1], big[2]-half[2])
        xnew = (half[1] + (half[1]-big[1])/15, half[2] + (half[2]-big[2])/15)  # local extrapolation
        fac = err > 0 ? 0.9 * (tol/err)^(1/5) : 2.0
        hnext = clamp(hh * clamp(fac, 0.2, 2.0), h_min, h_max)
        (err <= tol || hh <= h_min) && return (xnew, hnext, true)
        hh = hnext
    end
    return (_rk4(itp, x, hh, sgn), hh, false)
end

"""
    _seed_omp(itp::FI.AbstractInterpolant, c::T, RA::T, ZA::T, R_max::T; nscan::Int=200) where {T<:Real}

Outboard-midplane seed for the ψ=c contour: scan `R ∈ [RA, R_max]` at `Z = ZA` for the first
sign change of `ψ(R,ZA) - c`, bisect to a bracket, then project onto ψ=c. Returns
`((R,Z), found::Bool)`.
"""
function _seed_omp(itp::FI.AbstractInterpolant, c::T, RA::T, ZA::T, R_max::T; nscan::Int=200) where {T<:Real}
    Rs = range(RA, R_max; length=nscan)
    f(R) = itp(R, ZA) - c
    prevf = f(Rs[1]); a = Rs[1]
    for i in 2:nscan
        fi = f(Rs[i])
        if sign(fi) != sign(prevf)
            lo, hi = Rs[i-1], Rs[i]
            for _ in 1:60
                mid = (lo + hi) / 2
                (sign(f(mid)) == sign(f(lo))) ? (lo = mid) : (hi = mid)
            end
            return _project_to_level(itp, c, ((lo + hi) / 2, ZA))
        end
        prevf = fi; a = Rs[i]
    end
    return ((a, ZA), false)
end

"""
    trace_surface_cubic(itp::FI.AbstractInterpolant, c::T, RA::T, ZA::T, R_max::T;
                        npoints::Int=361, kw...) where {T<:Real}

Trace a single closed flux surface ψ=c on the cubic interpolant: outboard-midplane seed →
predictor–corrector trace → uniform-arclength resample to `npoints` distinct points →
close and reorder via `reorder_flux_surface!` (default `force_close=true`).

The returned vectors have length `npoints+1`: `npoints` distinct arclength-uniform points
plus a closing duplicate (`r[end] == r[1]`, `z[end] == z[1]`), reordered so the outboard
midplane (OMP) point is first and the polygon is clockwise. This matches the Contour-path
`FluxSurface` representation, so `MXH` and arclength integrals (`int_fluxexpansion_dl`) work
correctly on the output.

Returns `(r, z, closed::Bool)`. Standalone — not wired into `trace_surfaces`.
"""
function trace_surface_cubic(itp::FI.AbstractInterpolant, c::T, RA::T, ZA::T, R_max::T;
    npoints::Int=361, kw...) where {T<:Real}
    seed, ok = _seed_omp(itp, c, RA, ZA, R_max)
    ok || return (T[], T[], false)
    Rs, Zs, closed = _trace_surface_cubic(itp, c, seed; kw...)
    closed || return (Rs, Zs, false)
    R, Z = _resample_contour(Rs, Zs, npoints)
    reorder_flux_surface!(R, Z, RA, ZA)   # close (first==last), reorder OMP-first, clockwise — matches the Contour FluxSurface representation
    return (R, Z, true)
end

"""
    trace_surfaces_cubic(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T},
                         wall_z::AbstractVector{T}; refine_extrema::Bool=true, npoints::Int=361, kw...) where {T<:Real}

Cubic-interpolant drop-in mirror of [`trace_surfaces`](@ref)`(eqt, wall_r, wall_z)`: same
high-level signature, returns the same `Vector{FluxSurface{T}}`, so the two can be compared
1:1 (e.g. `trace_surfaces(eqt, fw.r, fw.z)` vs `trace_surfaces_cubic(eqt, fw.r, fw.z)`).
Each ψ level is traced by predictor–corrector continuation on the bicubic interpolant instead
of marching squares. `kw...` forwards to the tracer (`method=:pc|:rk4`, `ε`, `rk4_tol`, …).

`wall_r`/`wall_z` are accepted for signature parity but currently unused: every level is
traced as a closed surface (open/wall-clipping is not yet implemented in the cubic path), so a
level that does not close raises an error. Standalone — NOT wired into `trace_surfaces`.
"""
function trace_surfaces_cubic(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T},
    wall_z::AbstractVector{T}; refine_extrema::Bool=true, npoints::Int=361, kw...) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, _, itp = ψ_interpolant(eqt2d)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z
    return trace_surfaces_cubic(eqt.profiles_1d.psi, eqt.profiles_1d.f, RA, ZA, maximum(r), itp;
        refine_extrema, npoints, kw...)
end

"""
    trace_surfaces_cubic(psis, f, RA, ZA, R_max, itp; refine_extrema=false, npoints=361, kw...)

Low-level batch flux-surface tracer on the cubic interpolant: for each level in `psis`, trace
with [`trace_surface_cubic`](@ref) and populate a `FluxSurface` by reusing the existing
field-eval helpers (`arc_length`, `Br_Bz`, `trapz`, `fluxsurface_extrema`). The innermost
(k=1) surface is the usual artificial scaled-down copy of the second. When `refine_extrema`,
the four geometric extrema (`max_r`/`min_r`/`max_z`/`min_z`) of each real surface are refined
with [`_robust_refine_extremum!`](@ref) (globalized X-point-aware Newton on the interpolant) —
the cubic analogue of the `refine_extrema` step in [`trace_surfaces`](@ref). Standalone.
"""
function trace_surfaces_cubic(psis::AbstractVector{T}, f::AbstractVector{T}, RA::T, ZA::T,
    R_max::T, itp::FI.AbstractInterpolant; refine_extrema::Bool=false, npoints::Int=361, kw...) where {T<:Real}
    N = length(psis)
    surfaces = Vector{FluxSurface{T}}(undef, N)
    axis = (RA, ZA)
    lo = (first(itp.grids[1]), first(itp.grids[2]))   # ψ grid domain box — clamp the refine search
    hi = (last(itp.grids[1]), last(itp.grids[2]))
    H = Matrix{T}(undef, 2, 2)   # one 2×2 Hessian scratch shared by every _robust_refine_extremum! call
    for k in N:-1:1
        if k == 1
            pr = (surfaces[2].r .- RA) ./ 100 .+ RA
            pz = (surfaces[2].z .- ZA) ./ 100 .+ ZA
        else
            pr, pz, closed = trace_surface_cubic(itp, psis[k], RA, ZA, R_max; npoints, kw...)
            if !closed
                # A level passing through the X-point (separatrix) cannot be closed by the PC
                # tracer — the cubic separatrix/open-surface layer is not implemented. Retry a
                # hair toward the axis as a proxy and warn, so a full-profile call still runs.
                psi_in = psis[k] + T(1e-3) * (psis[1] - psis[k])
                pr, pz, closed = trace_surface_cubic(itp, psi_in, RA, ZA, R_max; npoints, kw...)
                closed && @warn "trace_surfaces_cubic: ψ level $k of $N did not close; traced a slightly inner proxy (separatrix/open handling not implemented in the cubic path)"
            end
            closed || error("trace_surfaces_cubic: failed to close surface $k of $N at ψ=$(psis[k])")
        end
        ll = arc_length(pr, pz)
        Br, Bz = Br_Bz(itp, pr, pz)
        Bp2 = Br .^ 2 .+ Bz .^ 2
        Bp_abs = sqrt.(Bp2)
        Bp = Bp_abs .* sign.((pz .- ZA) .* Br .- (pr .- RA) .* Bz)
        Btot = sqrt.(Bp2 .+ (f[k] ./ pr) .^ 2)
        fluxexpansion = 1.0 ./ Bp_abs
        int_fluxexpansion_dl = trapz(ll, fluxexpansion)
        (_, _, _, _, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) =
            fluxsurface_extrema(pr, pz)
        s = FluxSurface(psis[k], pr, pz, r_at_max_z, max_z, r_at_min_z, min_z,
            z_at_max_r, max_r, z_at_min_r, min_r, Br, Bz, Bp, Btot, ll, fluxexpansion, int_fluxexpansion_dl)
        if refine_extrema && k != 1   # globalized X-point-aware refinement of the 4 geometric extrema (cubic analogue of trace_surfaces' refine_extrema)
            (s.max_r, s.z_at_max_r) = _robust_refine_extremum!(H, itp, psis[k], (s.max_r, s.z_at_max_r), :R, axis; lo, hi)
            (s.min_r, s.z_at_min_r) = _robust_refine_extremum!(H, itp, psis[k], (s.min_r, s.z_at_min_r), :R, axis; lo, hi)
            (s.r_at_max_z, s.max_z) = _robust_refine_extremum!(H, itp, psis[k], (s.r_at_max_z, s.max_z), :Z, axis; lo, hi)
            (s.r_at_min_z, s.min_z) = _robust_refine_extremum!(H, itp, psis[k], (s.r_at_min_z, s.min_z), :Z, axis; lo, hi)
        end
        surfaces[k] = s
    end
    return surfaces
end

"""
    _find_xpoint(itp::FI.AbstractInterpolant, seed::Tuple{T,T}; tol=1e-10, maxit=50) where {T<:Real}

Locate a critical point of ψ (`∇ψ=0`) by globalized 2-D Newton ([`_damped_newton2d`](@ref) +
[`_critical_eqs`](@ref)) from `seed`, and classify it by the Hessian determinant: `:saddle`
(X-point, `det H < 0`) or `:extremum` (O-point/axis, `det H > 0`). Returns
`(point, kind::Symbol, converged::Bool)`. Unclamped (no grid box) so an out-of-domain seed that
converges outside the grid is rejected by the domain guard below.
"""
function _find_xpoint(itp::FI.AbstractInterpolant, seed::Tuple{T,T}; tol::Real=1e-10, maxit::Int=50) where {T<:Real}
    H = Matrix{T}(undef, 2, 2)   # 2×2 Hessian scratch (Newton residual + saddle/extremum classification)
    pt, ok = _damped_newton2d(_critical_eqs(H, itp), seed; tol, maxit)
    ok || return (pt, :none, false)
    g1, g2 = itp.grids[1], itp.grids[2]                       # interpolant grid extents
    (first(g1) <= pt[1] <= last(g1) && first(g2) <= pt[2] <= last(g2)) || return (pt, :none, false)
    FI.hessian!(H, itp, pt)
    detH = H[1, 1] * H[2, 2] - H[1, 2]^2
    return (pt, detH < 0 ? :saddle : :extremum, true)
end

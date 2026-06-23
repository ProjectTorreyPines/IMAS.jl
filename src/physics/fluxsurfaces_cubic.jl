# Cubic-interpolant-based flux-surface geometry.
#
# Flux-surface geometry computed directly from the ψ cubic interpolant (analytic
# gradient/Hessian via FastInterpolations, aliased `FI` in IMAS.jl) with Newton
# iteration — as opposed to the Contour-based tracing in fluxsurfaces.jl.
#
# Currently: X-point-aware refinement of a traced surface's geometric extrema
# (max_r/min_r/max_z/min_z). Intended to also host a general cubic-interpolant-
# based flux-surface tracer. All helpers here are private (no public API yet).

"""
    _newton2d(residual_jacobian::F, seed::Tuple{T,T}; reference::Union{Nothing,Tuple{T,T}}=nothing,
              factor::Real=0.9, tol::Real=1e-10, maxit::Int=30) where {F,T<:Real}

Generic damped 2×2 Newton solver for `F(R,Z) = 0`. `residual_jacobian(R, Z)` returns
`(F1, F2, J11, J12, J21, J22)` — the residual and its 2×2 Jacobian at `(R, Z)`. Passing
the condition in lets the same solver find a flux-surface extremum (via
[`_extremum_residual`](@ref)) or a critical point of ψ (via [`_critical_residual`](@ref)).

If `reference` is given, each step's displacement is capped below `factor ×` the distance
to it, so the iterate cannot cross over that point (used to stay on one side of an
X-point). Returns `((R, Z), converged::Bool)`.
"""
function _newton2d(residual_jacobian::F, seed::Tuple{T,T}; reference::Union{Nothing,Tuple{T,T}}=nothing,
    factor::Real=0.9, tol::Real=1e-10, maxit::Int=30) where {F,T<:Real}
    R, Z = seed
    for _ in 1:maxit
        F1, F2, J11, J12, J21, J22 = residual_jacobian(R, Z)
        (abs(F1) <= tol && abs(F2) <= tol) && return ((R, Z), true)
        det = J11 * J22 - J12 * J21
        (isfinite(det) && abs(det) > eps(T)) || return ((R, Z), false)  # degenerate Jacobian
        dR = (J22 * F1 - J12 * F2) / det
        dZ = (J11 * F2 - J21 * F1) / det
        if reference !== nothing
            maxd = factor * hypot(R - reference[1], Z - reference[2])
            st = hypot(dR, dZ)
            st > maxd && (dR *= maxd / st; dZ *= maxd / st)
        end
        R -= dR
        Z -= dZ
        (isfinite(R) && isfinite(Z)) || return ((R, Z), false)
    end
    return ((R, Z), false)
end

# residual + 2×2 Jacobian for the extremum system {ψ = target_psi, ∂ψ/∂(daxis) = 0}
# (daxis = 2 -> ∂ψ/∂Z = 0 for R-extrema; daxis = 1 -> ∂ψ/∂R = 0 for Z-extrema)
_extremum_residual(itp::FI.AbstractInterpolant, target_psi::Real, daxis::Int) =
    (R, Z) -> begin
        val, g = FI.value_gradient(itp, (R, Z))
        H = FI.hessian(itp, (R, Z))
        return (val - target_psi, g[daxis], g[1], g[2], H[daxis, 1], H[daxis, 2])
    end

# residual + 2×2 Jacobian for the critical-point system ∇ψ = 0 (Jacobian = Hessian)
_critical_residual(itp::FI.AbstractInterpolant) =
    (R, Z) -> begin
        g = FI.gradient(itp, (R, Z))
        H = FI.hessian(itp, (R, Z))
        return (g[1], g[2], H[1, 1], H[1, 2], H[2, 1], H[2, 2])
    end

"""
    _refine_extremum(itp::FI.AbstractInterpolant, target_psi::T, seed::Tuple{T,T},
                     extremum_of::Symbol, axis::Tuple{T,T};
                     factor::Real=0.9, tol::Real=1e-10, maxit::Int=30) where {T<:Real}

Refine a flux-surface geometric extremum from a rough `seed = (R, Z)`, using the analytic
gradient and Hessian of the cubic interpolant `itp`.

`extremum_of = :R` finds an extremum of `R` (`max_r`/`min_r`) by enforcing `∂ψ/∂Z = 0`;
`extremum_of = :Z` finds an extremum of `Z` (`max_z`/`min_z`) by enforcing `∂ψ/∂R = 0`.
The seed selects which root (e.g. outboard vs inboard) is found, since both satisfy the
same 2×2 system `{ψ(R,Z) = target_psi, ∂ψ/∂n(R,Z) = 0}`.

Fast path: a plain Newton ([`_newton2d`](@ref)) on that system. Near an X-point the system
has two solutions and the plain Newton can converge to the wrong one (e.g. above the
X-point, outside the confined region). The candidate is therefore validated by the sign of
the constrained curvature of the extremized coordinate (`< 0` at a maximum, `> 0` at a
minimum); whether a max or min is sought is inferred from the candidate relative to the
magnetic `axis`. When the candidate is the wrong branch, the genuine extremum is recovered:

 1. find the nearby critical point of ψ (`∇ψ = 0`, an X-point saddle or the O-point) with
    [`_newton2d`](@ref) + [`_critical_residual`](@ref), seeded at the candidate;
 2. mirror the candidate across that critical point, back into the confined region;
 3. re-solve the extremum system with [`_newton2d`](@ref) using that critical point as the
    `reference`, so the bounded step cannot cross back over it.

Returns the genuine `(R, Z)`. Degenerate fallbacks: `seed` if the fast-path Newton cannot
converge (e.g. a degenerate Jacobian seeded at the magnetic axis); the candidate if no
critical point is found; the mirror point if the bounded re-solve does not land on a
genuine extremum (already a confined-side estimate).
"""
function _refine_extremum(itp::FI.AbstractInterpolant, target_psi::T, seed::Tuple{T,T},
    extremum_of::Symbol, axis::Tuple{T,T};
    factor::Real=0.9, tol::Real=1e-10, maxit::Int=30) where {T<:Real}
    daxis = extremum_of === :R ? 2 : extremum_of === :Z ? 1 :
            throw(ArgumentError("_refine_extremum: extremum_of must be :R or :Z, got :$extremum_of"))
    eaxis = extremum_of === :R ? 1 : 2

    # fast path: plain Newton on {ψ = target_psi, ∂ψ/∂n = 0}
    candidate, converged = _newton2d(_extremum_residual(itp, target_psi, daxis), seed; tol, maxit)
    converged || return seed   # degenerate Jacobian / no convergence -> give up

    # genuineness test: sign of the constrained curvature of the extremized coordinate
    # (genuine maximum has -ψ_dd/ψ_e < 0, minimum > 0); want_max inferred vs the axis
    want_max = candidate[eaxis] > axis[eaxis]
    function is_genuine(p::Tuple{T,T})
        _, g = FI.value_gradient(itp, p)
        H = FI.hessian(itp, p)
        curv = -H[daxis, daxis] / g[eaxis]
        return want_max ? curv < zero(T) : curv > zero(T)
    end
    is_genuine(candidate) && return candidate   # fast path already genuine

    # wrong branch near an X-point -> recover the genuine extremum
    # 1. nearby critical point of ψ (X-point saddle or O-point) via ∇ψ = 0
    crit, cok = _newton2d(_critical_residual(itp), candidate; tol, maxit)
    cok || return candidate

    # 2. mirror the wrong root across the critical point (back into the confined region)
    mir = (2crit[1] - candidate[1], 2crit[2] - candidate[2])

    # 3. bounded Newton from the mirror seed, capped so it cannot cross back over the X-point
    point, recovered = _newton2d(_extremum_residual(itp, target_psi, daxis), mir;
        reference=crit, factor, tol, maxit)

    # accept only a genuine extremum; else fall back to the mirror point (already a
    # confined-side estimate of the extremum)
    (recovered && is_genuine(point)) && return point
    return mir
end

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
    val, _ = FI.value_gradient(itp, (R, Z))
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
"""
function _contour_curvature(itp::FI.AbstractInterpolant, x::Tuple{T,T}) where {T<:Real}
    g = FI.gradient(itp, x)
    H = FI.hessian(itp, x)
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
Returns `(xnew, on_surface::Bool)`.
"""
function _step_pc(itp::FI.AbstractInterpolant, c::T, x::Tuple{T,T}, h::Real, sgn::Int;
    corr_tol::Real=1e-10, corr_max::Int=8) where {T<:Real}
    t = _contour_tangent(itp, x, sgn)
    np = (-t[2], t[1])                            # principal normal (rotate tangent +90°)
    κ = _contour_curvature(itp, x)
    xp = (x[1] + h * t[1] + (h^2 / 2) * κ * np[1],
          x[2] + h * t[2] + (h^2 / 2) * κ * np[2])
    return _project_to_level(itp, c, xp; tol=corr_tol, maxit=corr_max)
end

"""
    _contour_step(itp::FI.AbstractInterpolant, x::Tuple{T,T}; ε::Real=1e-6,
                  h_min::Real, h_max::Real, max_turn::Real=deg2rad(10), κ_floor::Real=1e-8) where {T<:Real}

Arclength step for the contour at `x`: bound the predictor chord error `≈ ½|κ|h² ≤ ε` and
the per-step turning angle `|κ|·h ≤ max_turn`, then clamp to `[h_min, h_max]`. `κ_floor`
prevents division by zero on near-straight pieces.
"""
function _contour_step(itp::FI.AbstractInterpolant, x::Tuple{T,T}; ε::Real=1e-6,
    h_min::Real, h_max::Real, max_turn::Real=deg2rad(10), κ_floor::Real=1e-8) where {T<:Real}
    κ = max(abs(_contour_curvature(itp, x)), κ_floor)
    h = sqrt(2 * ε / κ)            # chord-error bound
    h = min(h, max_turn / κ)       # turning-angle cap
    return clamp(h, h_min, h_max)
end

# signed angle from vector a to vector b (radians, in (-π, π])
_signed_angle(a::Tuple, b::Tuple) = atan(a[1]*b[2] - a[2]*b[1], a[1]*b[2]*0 + a[1]*b[1] + a[2]*b[2])

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

    for _ in 1:max_steps
        # near-critical-point guard: when |∇ψ| < τ_grad (close to an X-point saddle), creep
        # at the smallest step via the corrector so the surface traces the confined side
        # instead of leaking through the saddle (τ_grad=0 disables this, default behavior).
        gg = τ_grad > 0 ? FI.gradient(itp, x) : nothing
        if gg !== nothing && hypot(gg[1], gg[2]) < τ_grad
            xnew, sok = _step_pc(itp, c, x, h_min, sgn)
        elseif method === :rk4
            xnew, hrk, sok = _step_rk4_adaptive(itp, x, hrk, sgn; tol=rk4_tol, h_min, h_max)
        else
            h = _contour_step(itp, x; ε, h_min, h_max, max_turn, κ_floor)
            xnew, sok = _step_pc(itp, c, x, h, sgn)
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
predictor–corrector trace → uniform-arclength resample to `npoints` → reorder CCW from the
outboard midplane (`reorder_flux_surface!`). Returns `(r, z, closed::Bool)`. Standalone —
not wired into `trace_surfaces`.
"""
function trace_surface_cubic(itp::FI.AbstractInterpolant, c::T, RA::T, ZA::T, R_max::T;
    npoints::Int=361, kw...) where {T<:Real}
    seed, ok = _seed_omp(itp, c, RA, ZA, R_max)
    ok || return (T[], T[], false)
    Rs, Zs, closed = _trace_surface_cubic(itp, c, seed; kw...)
    closed || return (Rs, Zs, false)
    R, Z = _resample_contour(Rs, Zs, npoints)
    reorder_flux_surface!(R, Z, RA, ZA; force_close=false)   # CW from OMP; half-open rep (no duplicate endpoint)
    return (R, Z, true)
end

"""
    trace_surfaces_cubic(psis, f, RA, ZA, R_max, itp; npoints=361, kw...)

Batch flux-surface tracer on the cubic interpolant: for each level in `psis`, trace with
[`trace_surface_cubic`](@ref) and populate a `FluxSurface` by reusing the existing field-eval
helpers (`arc_length`, `Br_Bz`, `trapz`, `fluxsurface_extrema`). The innermost (k=1) surface
is the usual artificial scaled-down copy of the second. Standalone — NOT wired into
`trace_surfaces`; used to validate flux-surface averages against the Contour path.
"""
function trace_surfaces_cubic(psis::AbstractVector{T}, f::AbstractVector{T}, RA::T, ZA::T,
    R_max::T, itp::FI.AbstractInterpolant; npoints::Int=361, kw...) where {T<:Real}
    N = length(psis)
    surfaces = Vector{FluxSurface{T}}(undef, N)
    for k in N:-1:1
        if k == 1
            pr = (surfaces[2].r .- RA) ./ 100 .+ RA
            pz = (surfaces[2].z .- ZA) ./ 100 .+ ZA
        else
            pr, pz, closed = trace_surface_cubic(itp, psis[k], RA, ZA, R_max; npoints, kw...)
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
        surfaces[k] = FluxSurface(psis[k], pr, pz, r_at_max_z, max_z, r_at_min_z, min_z,
            z_at_max_r, max_r, z_at_min_r, min_r, Br, Bz, Bp, Btot, ll, fluxexpansion, int_fluxexpansion_dl)
    end
    return surfaces
end

"""
    _find_xpoint(itp::FI.AbstractInterpolant, seed::Tuple{T,T}; tol=1e-10, maxit=50) where {T<:Real}

Locate a critical point of ψ (`∇ψ=0`) by 2-D Newton ([`_newton2d`](@ref) + [`_critical_residual`](@ref))
from `seed`, and classify it by the Hessian determinant: `:saddle` (X-point, `det H < 0`) or
`:extremum` (O-point/axis, `det H > 0`). Returns `(point, kind::Symbol, converged::Bool)`.
"""
function _find_xpoint(itp::FI.AbstractInterpolant, seed::Tuple{T,T}; tol::Real=1e-10, maxit::Int=50) where {T<:Real}
    pt, ok = _newton2d(_critical_residual(itp), seed; tol, maxit)
    ok || return (pt, :none, false)
    g1, g2 = itp.grids[1], itp.grids[2]                       # interpolant grid extents
    (first(g1) <= pt[1] <= last(g1) && first(g2) <= pt[2] <= last(g2)) || return (pt, :none, false)
    H = FI.hessian(itp, pt)
    detH = H[1, 1] * H[2, 2] - H[1, 2]^2
    return (pt, detH < 0 ? :saddle : :extremum, true)
end

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
    domain=nothing) where {T<:Real}

    x0, ok = _project_to_level(itp, c, seed)
    ok || return (T[seed[1]], T[seed[2]], false)
    t0 = _contour_tangent(itp, x0, sgn)
    Rs = T[x0[1]]; Zs = T[x0[2]]
    x = x0; tprev = t0; Θ = zero(T); s = zero(T)

    for _ in 1:max_steps
        h = _contour_step(itp, x; ε, h_min, h_max, max_turn, κ_floor)
        xnew, sok = _step_pc(itp, c, x, h, sgn)
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

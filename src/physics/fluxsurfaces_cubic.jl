# Cubic-interpolant-based flux-surface geometry.
#
# Flux-surface geometry computed directly from the ψ cubic interpolant (analytic
# gradient/Hessian via FastInterpolations, aliased `FI` in IMAS.jl) with Newton
# iteration — as opposed to the Contour-based tracing in fluxsurfaces.jl.
#
# Currently: fast/slow-path refinement of a traced surface's geometric extrema
# (max_r/min_r/max_z/min_z). Intended to also host a general cubic-interpolant-
# based flux-surface tracer. All helpers here are private (no public API yet).

"""
    _newton2d(residual_jacobian, seed::Tuple{T,T}; reference::Union{Nothing,Tuple{T,T}}=nothing,
              factor::Real=0.9, tol::Real=1e-10, maxit::Int=30) where {T<:Real}

Generic damped 2×2 Newton solver for `F(R,Z) = 0`. `residual_jacobian(R, Z)` returns
`(F1, F2, J11, J12, J21, J22)` — the residual and its 2×2 Jacobian at `(R, Z)`. Passing
the condition in lets the same solver find a flux-surface extremum (via
[`_extremum_residual`](@ref)) or a critical point of ψ (via [`_critical_residual`](@ref)).

If `reference` is given, each step's displacement is capped below `factor ×` the distance
to it, so the iterate cannot cross over that point (used to stay on one side of an
X-point). Returns `((R, Z), converged::Bool)`.
"""
function _newton2d(residual_jacobian, seed::Tuple{T,T}; reference::Union{Nothing,Tuple{T,T}}=nothing,
    factor::Real=0.9, tol::Real=1e-10, maxit::Int=30) where {T<:Real}
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
_extremum_residual(itp, target_psi, daxis::Int) =
    (R, Z) -> begin
        val, g = FI.value_gradient(itp, (R, Z))
        H = FI.hessian(itp, (R, Z))
        return (val - target_psi, g[daxis], g[1], g[2], H[daxis, 1], H[daxis, 2])
    end

# residual + 2×2 Jacobian for the critical-point system ∇ψ = 0 (Jacobian = Hessian)
_critical_residual(itp) =
    (R, Z) -> begin
        _, g = FI.value_gradient(itp, (R, Z))
        H = FI.hessian(itp, (R, Z))
        return (g[1], g[2], H[1, 1], H[1, 2], H[2, 1], H[2, 2])
    end

"""
    _refine_extremum(itp, target_psi::T, seed::Tuple{T,T}, extremum_of::Symbol; tol::Real=1e-10, maxit::Int=30) where {T<:Real}

Fast-path refinement of a flux-surface geometric extremum: solve the 2×2 system
`{ψ(R,Z) = target_psi, ∂ψ/∂n(R,Z) = 0}` with a plain Newton ([`_newton2d`](@ref)) seeded
at `seed = (R, Z)`, using the analytic gradient and Hessian of the cubic interpolant `itp`.

`extremum_of = :R` finds an extremum of `R` (`max_r`/`min_r`) by enforcing `∂ψ/∂Z = 0`;
`extremum_of = :Z` finds an extremum of `Z` (`max_z`/`min_z`) by enforcing `∂ψ/∂R = 0`.
The seed selects which root (e.g. outboard vs inboard) is found, since both satisfy
the same system.

Returns the refined `(R, Z)`, or `seed` if Newton fails to converge or hits a
degenerate Jacobian. Near an X-point the system has two solutions and this plain Newton
can converge to the wrong one — use [`_refine_extremum_bounded`](@ref) to verify and
recover in that case.
"""
function _refine_extremum(itp, target_psi::T, seed::Tuple{T,T}, extremum_of::Symbol; tol::Real=1e-10, maxit::Int=30) where {T<:Real}
    daxis = extremum_of === :R ? 2 : extremum_of === :Z ? 1 :
            throw(ArgumentError("_refine_extremum: extremum_of must be :R or :Z, got :$extremum_of"))
    point, converged = _newton2d(_extremum_residual(itp, target_psi, daxis), seed; tol, maxit)
    return converged ? point : seed
end

"""
    _refine_extremum_bounded(itp, target_psi::T, candidate::Tuple{T,T}, extremum_of::Symbol, axis::Tuple{T,T};
                             factor::Real=0.9, tol::Real=1e-10, maxit::Int=30) where {T<:Real}

Slow-path recovery for [`_refine_extremum`](@ref): verify that `candidate` (a root of
`{ψ=target_psi, ∂ψ/∂n=0}` returned by the fast path) is the *genuine* extremum and, if
not, recover the correct one.

The genuineness test is the sign of the constrained curvature of the extremized
coordinate along the surface (`< 0` at a maximum, `> 0` at a minimum); whether a maximum
or minimum is sought is inferred from `candidate` relative to the magnetic `axis`. When
the test fails (the fast path landed on the wrong branch near an X-point), the routine:

 1. finds the nearby critical point of ψ (`∇ψ = 0`, an X-point or the O-point) with
    [`_newton2d`](@ref) + [`_critical_residual`](@ref), seeded at `candidate`;
 2. mirrors `candidate` across that critical point, back into the confined region;
 3. re-solves the extremum system with [`_newton2d`](@ref) using that critical point as
    the `reference`, so the bounded step cannot cross back over it.

Returns the genuine `(R, Z)`; if the bounded solve does not land on a genuine extremum
it returns the mirror point (already a confined-side estimate); if no critical point is
found it returns `candidate` unchanged.
"""
function _refine_extremum_bounded(itp, target_psi::T, candidate::Tuple{T,T}, extremum_of::Symbol, axis::Tuple{T,T};
    factor::Real=0.9, tol::Real=1e-10, maxit::Int=30) where {T<:Real}
    daxis = extremum_of === :R ? 2 : extremum_of === :Z ? 1 :
            throw(ArgumentError("_refine_extremum_bounded: extremum_of must be :R or :Z, got :$extremum_of"))
    eaxis = extremum_of === :R ? 1 : 2
    want_max = candidate[eaxis] > axis[eaxis]

    # constrained-curvature sign test: genuine maximum has -ψ_dd/ψ_e < 0 (minimum > 0)
    function is_genuine(p::Tuple{T,T})
        _, g = FI.value_gradient(itp, p)
        H = FI.hessian(itp, p)
        curv = -H[daxis, daxis] / g[eaxis]
        return want_max ? curv < zero(T) : curv > zero(T)
    end

    is_genuine(candidate) && return candidate   # already the genuine extremum

    # 1. nearby critical point of ψ (X-point saddle or O-point) via ∇ψ = 0
    crit, cok = _newton2d(_critical_residual(itp), candidate; tol, maxit)
    cok || return candidate

    # 2. mirror the wrong root across the critical point (back into the confined region)
    mir = (2crit[1] - candidate[1], 2crit[2] - candidate[2])

    # 3. bounded Newton from the mirror seed, capped so it cannot cross back over the X-point
    point, converged = _newton2d(_extremum_residual(itp, target_psi, daxis), mir;
        reference=crit, factor, tol, maxit)

    # accept only a genuine extremum; else fall back to the mirror point (already a
    # confined-side estimate of the extremum)
    (converged && is_genuine(point)) && return point
    return mir
end

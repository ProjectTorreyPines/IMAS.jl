using LinearAlgebra
import SimpleNonlinearSolve

document[Symbol("Physics fields")] = Symbol[]

# gradient of a ψ interpolant: FI gets the native fast path,
# the open version is for users still passing an Interpolations.jl interpolant
@inline _psi_gradient(itp::FI.AbstractInterpolant, r, z) = FI.gradient(itp, (r, z))
@inline _psi_gradient(itp, r, z) = parentmodule(typeof(itp)).gradient(itp, r, z)

# value+gradient and Hessian with the same FI-fast / open-fallback split, so the Newton extremum
# refine stays backend-agnostic. Open value_gradient emulates (value, gradient); open hessian!
# copies the backend's returned Hessian into the caller's buffer.
@inline _psi_value_gradient(itp::FI.AbstractInterpolant, r, z) = FI.value_gradient(itp, (r, z))
@inline _psi_value_gradient(itp, r, z) = (itp(r, z), parentmodule(typeof(itp)).gradient(itp, r, z))
@inline _psi_hessian!(H, itp::FI.AbstractInterpolant, r, z) = FI.hessian!(H, itp, (r, z))
@inline _psi_hessian!(H, itp, r, z) = (H .= parentmodule(typeof(itp)).hessian(itp, r, z); H)

"""
    Br_Bz(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)

Returns Br and Bz named tuple evaluated at r and z starting from ψ interpolant
"""
function Br_Bz(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    return Br_Bz_meshgrid(PSI_interpolant, r, z)
end

"""
    Br_Bz(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
"""
function Br_Bz(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return Br_Bz(eqt2d)
end

"""
    Br_Bz(PSI_interpolant, r::T, z::T) where {T<:Real}
"""
function Br_Bz(PSI_interpolant, r::T, z::T) where {T<:Real}
    grad = _psi_gradient(PSI_interpolant, r, z)
    inv_twopi_r = 1.0 / (2π * r)
    Br = grad[2] * inv_twopi_r
    Bz = -grad[1] * inv_twopi_r
    return (Br=Br, Bz=Bz)
end

"""
    Br_Bz(PSI_interpolant, r::AbstractArray{T}, z::AbstractArray{T}) where {T<:Real}
"""
function Br_Bz(PSI_interpolant, r::AbstractArray{T}, z::AbstractArray{T}) where {T<:Real}
    @assert size(r) == size(z)
    Br, Bz = similar(r), similar(r)
    for k in eachindex(r)
        Br[k], Bz[k] = Br_Bz(PSI_interpolant, r[k], z[k])
    end
    return (Br=Br, Bz=Bz)
end

@compat public Br_Bz
push!(document[Symbol("Physics fields")], :Br_Bz)

"""
    Br_Bphi_Bz(PSI_interpolant, B0::Real, R0::Real, r::Real, phi::Real, z::Real)

Returns Br, Bphi, Bz named tuple evaluated at (r, z) starting from ψ interpolant and B0 and R0.
`phi` is accepted for signature symmetry with the (x, y, z) callers but is unused: the
equilibrium is axisymmetric, so all components are independent of the toroidal angle.
"""
function Br_Bphi_Bz(PSI_interpolant,
    B0::Real, R0::Real, r::Real, phi::Real, z::Real)

    Br, Bz = IMAS.Br_Bz(PSI_interpolant, r, z)
    Bphi = B0 * R0 / r

    return (Br=Br, Bphi=Bphi, Bz=Bz) #B-field components in this order in accordance with cocos 11
end

@compat public Br_Bphi_Bz
push!(document[Symbol("Physics fields")], :Br_Bphi_Bz)

"""
    Bx_By_Bz(PSI_interpolant, B0::Real, R0::Real, x::Real, y::Real, z::Real)

Returns Bx, By, Bz named tuple evaluated at x,y,z starting from ψ interpolant and B0 and R0
"""
function Bx_By_Bz(PSI_interpolant,
    B0::Real, R0::Real, x::Real, y::Real, z::Real)

    r = hypot(x, y)
    phi = atan(y, x)

    Br, Bphi, Bz = IMAS.Br_Bphi_Bz(PSI_interpolant, B0, R0, r, phi, z)

    sp, cp = sincos(phi)
    Bx = Br * cp - Bphi * sp
    By = Br * sp + Bphi * cp

    return (Bx=Bx, By=By, Bz=Bz)
end

@compat public Bx_By_Bz
push!(document[Symbol("Physics fields")], :Bx_By_Bz)

"""
    Br_Bz_meshgrid(PSI_interpolant, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
"""
function Br_Bz_meshgrid(PSI_interpolant, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
    Br = Matrix{T}(undef, length(r), length(z))
    Bz = Matrix{T}(undef, length(r), length(z))
    for kr in eachindex(r)
        for kz in eachindex(z)
            Br[kr, kz], Bz[kr, kz] = Br_Bz(PSI_interpolant, r[kr], z[kz])
        end
    end
    return (Br=Br, Bz=Bz)
end

"""
    Bp(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)

Returns Bp evaluated at r and z starting from ψ interpolant
"""
function Bp(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    return Bp_meshgrid(PSI_interpolant, r, z)
end

"""
    Bp(PSI_interpolant, r::T, z::T) where {T<:Real}
"""
function Bp(PSI_interpolant, r::T, z::T) where {T<:Real}
    Br, Bz = Br_Bz(PSI_interpolant, r, z)
    return sqrt(Br^2 + Bz^2)
end

"""
    Bp_meshgrid(PSI_interpolant, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
"""
function Bp_meshgrid(PSI_interpolant, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
    Bp = Matrix{T}(undef, length(r), length(z))
    for kr in eachindex(r)
        for kz in eachindex(z)
            Bp[kr, kz] = Bp(PSI_interpolant, r[kr], z[kz])
        end
    end
    return Bp
end

"""
    Bp(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
"""
function Bp(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return Bp(eqt2d)
end

@compat public Bp
push!(document[Symbol("Physics fields")], :Bp)

"""
    _toroidal_angle(p, q)

Unsigned angle swept in the (x, y) plane between two consecutive points `p` and
`q`, via `atan(cross, dot)`. Shared by both field-line integrators so they
measure `Δϕ` identically. Unlike `acos(dot / (‖p‖‖q‖))` this needs no
normalization and no domain clamp: it is exact for nearly-collinear points and
returns 0 when a point lies on the axis (r = 0), instead of `NaN`.
"""
function _toroidal_angle(p, q)
    cross = p[1] * q[2] - p[2] * q[1]
    dot2d = p[1] * q[1] + p[2] * q[2]
    return abs(atan(cross, dot2d))
end

"""
    ImplicitMidpointState{T,F<:Function}(vector_field, current_point, step_size, count, Δϕ)

State for the implicit-midpoint field-line integrator. Each step solves the
implicit midpoint equation `next = current + h·vector_field((current+next)/2)`
for the next point. (The RK4 path needs no state object — see [`_trace_rk4`](@ref).)

# Fields

  - `vector_field`: Function representing the vector field to be integrated
  - `current_point`: Current position in the vector field
  - `step_size`: Step size for numerical integration
  - `count`: Step counter
  - `Δϕ`: Toroidal angle traversed
"""
mutable struct ImplicitMidpointState{T,F<:Function}
    vector_field::F
    current_point::Vector{T}
    step_size::T
    count::Int
    Δϕ::T
end

"""
    (obj::ImplicitMidpointState)(next_point)

Residual of the implicit midpoint equation, driven to zero each step to find
`next_point`.
"""
function (obj::ImplicitMidpointState)(next_point)
    midpoint = 0.5 * (obj.current_point .+ next_point)
    return next_point .- obj.current_point .- obj.step_size * obj.vector_field(midpoint)
end

"""
    _next!(obj::ImplicitMidpointState)

Advance one step with the implicit midpoint rule (2nd-order, symplectic): solve
`obj(next) = 0`, then record the new point and the toroidal angle swept.
"""
function _next!(obj::ImplicitMidpointState)
    # Solve the implicit-midpoint equation obj(next)=0 (copy guess to avoid aliasing)
    prob = SimpleNonlinearSolve.NonlinearProblem((u, _p) -> obj(u), copy(obj.current_point))
    sol = SimpleNonlinearSolve.solve(prob, SimpleNonlinearSolve.SimpleTrustRegion(); abstol=1e-10, reltol=1e-10)
    next_point = sol.u

    dphi = _toroidal_angle(obj.current_point, next_point)
    obj.current_point .= next_point
    obj.count += 1
    obj.Δϕ += dphi

    return next_point
end

"""
    _trace_rk4(vector_field, start_point, step_size, stop_condition)

Trace a field line with the explicit classical Runge–Kutta (RK4) method.
Self-contained: a single loop, no per-step nonlinear solve. `stop_condition`
receives a lightweight object exposing `current_point`, `count`, and `Δϕ`.
"""
function _trace_rk4(vector_field::VF, start_point, step_size, stop_condition::SC) where {VF,SC}
    h = step_size
    p = copy(start_point)
    trajectory = [p]
    count = 0
    Δϕ = zero(h)
    while !stop_condition((current_point=p, count=count, Δϕ=Δϕ))
        k1 = vector_field(p)
        k2 = vector_field(p .+ (h / 2) .* k1)
        k3 = vector_field(p .+ (h / 2) .* k2)
        k4 = vector_field(p .+ h .* k3)
        pnext = p .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        dphi = _toroidal_angle(p, pnext)
        Δϕ += dphi
        count += 1
        p = pnext
        push!(trajectory, pnext)
    end
    return trajectory
end

"""
    trace_field_line(vector_field, start_point, step_size, stop_condition; method=:implicit_midpoint)

Compute a trajectory along a vector field.

# Arguments

  - `vector_field`: Vector field function
  - `start_point`: Initial point of the trajectory
  - `step_size`: Size of each integration step
  - `stop_condition`: Function determining when to stop integration. Receives an object exposing `current_point`, `count`, and `Δϕ`, and returns a `Bool`.

# Keywords

  - `method`: integration method, `:implicit_midpoint` (default; 2nd-order, symplectic) or `:rk4` (explicit 4th-order, faster)

# Returns

Trajectory represented as a vector of points
"""
function trace_field_line(vector_field::VF, start_point, step_size, stop_condition::SC; method::Symbol=:implicit_midpoint) where {VF,SC}
    if method === :implicit_midpoint
        return _trace_implicit_midpoint(vector_field, start_point, step_size, stop_condition)
    elseif method === :rk4
        return _trace_rk4(vector_field, start_point, step_size, stop_condition)
    else
        error("unknown field-line integration method :$(method); use :implicit_midpoint or :rk4")
    end
end

function _trace_implicit_midpoint(vector_field::VF, start_point, step_size, stop_condition::SC) where {VF,SC}
    # copy start_point so the caller's array is not mutated (matches _trace_rk4)
    integrator = ImplicitMidpointState(vector_field, copy(start_point), step_size, 0, zero(step_size))

    next_point = copy(start_point)
    trajectory = [next_point]

    while !stop_condition(integrator)
        next_point = _next!(integrator)
        push!(trajectory, next_point)
    end

    return trajectory
end

"""
    trace_field_line(eqt::IMAS.equilibrium__time_slice, r, z; phi=zero(r), step_size=0.01, max_turns=1, stop_condition=nothing, method=:implicit_midpoint)

Compute a field line trajectory for a specific equilibrium time slice.

# Arguments

  - `eqt`: Equilibrium time slice
  - `r`: Radial coordinate
  - `z`: Vertical coordinate
  - `phi`: Toroidal angle (default: 0)
  - `step_size`: Size of each integration step (default: 0.01)
  - `max_turns`: Maximum number of toroidal turns used in default stop condition (default: 1)
  - `stop_condition`: User defined stop condition that takes an object exposing `current_point`, `count`, and `Δϕ` and returns a Boolean (default: nothing)
  - `method`: integration method, `:implicit_midpoint` (default) or `:rk4` (explicit, faster)

# Returns

Named tuple with x, y, z, r, and phi coordinates of the trajectory
"""
function trace_field_line(eqt::IMAS.equilibrium__time_slice, r, z;
    phi=zero(r), step_size=0.01, max_turns=1,
    stop_condition::Union{Nothing,Function}=nothing, method::Symbol=:implicit_midpoint)

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    rg, zg, PSI_interpolant = ψ_interpolant(eqt2d)
    R0 = eqt.global_quantities.vacuum_toroidal_field.r0
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

    # Define vector field function (unit B)
    vector_field = (xyz) -> begin
        x = xyz[1]
        y = xyz[2]
        z = xyz[3]
        Bx, By, Bz = IMAS.Bx_By_Bz(PSI_interpolant, B0, R0, x, y, z)
        B = [Bx, By, Bz]
        B / norm(B)
    end

    # Define start point
    x = r * cos(phi)
    y = r * sin(phi)
    start_point = [x, y, z]

    # Define stop condition
    if stop_condition === nothing
        stop = (obj) -> begin
            turns = floor(Int, obj.Δϕ / (2pi))
            x1, y1, z1 = obj.current_point
            r1 = hypot(x1, y1)
            r1 < rg[1] || r1 > rg[end] || z1 < zg[1] || z1 > zg[end] || turns >= max_turns
        end
    else
        stop = stop_condition
    end

    # Calculate field line trajectory
    trajectory = trace_field_line(vector_field, start_point, step_size, stop; method)

    # Extract positions and return named tuple
    xt = getindex.(trajectory, 1)
    yt = getindex.(trajectory, 2)
    zt = getindex.(trajectory, 3)
    rt = hypot.(xt, yt)
    phit = atan.(yt, xt)

    return (x=xt, y=yt, z=zt, r=rt, phi=phit)
end

@compat public trace_field_line
push!(document[Symbol("Physics fields")], :trace_field_line)

using LinearAlgebra
import NLsolve
import Interpolations

document[Symbol("Physics fields")] = Symbol[]

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
    Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::T, z::T) where {T<:Real}
"""
function Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::T, z::T) where {T<:Real}
    grad = Interpolations.gradient(PSI_interpolant, r, z)
    inv_twopi_r = 1.0 / (2π * r)
    Br = grad[2] * inv_twopi_r
    Bz = -grad[1] * inv_twopi_r
    return (Br=Br, Bz=Bz)
end

"""
    Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::AbstractArray{T}, z::AbstractArray{T}) where {T<:Real}
"""
function Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::AbstractArray{T}, z::AbstractArray{T}) where {T<:Real}
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
    Br_Bphi_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, B0::T, R0::T, r::T, z::T) where {T<:Real}

Returns Br, Bphi, Bz named tuple evaluated at r and z starting from ψ interpolant and B0 and R0
"""
function Br_Bphi_Bz(PSI_interpolant::Interpolations.AbstractInterpolation,
    B0::T, R0::T, r::T, phi::T, z::T) where {T<:Real}

    Br, Bz = IMAS.Br_Bz(PSI_interpolant, r, z)
    Bphi = B0 * R0 / r

    return (Br=Br, Bphi=Bphi, Bz=Bz) #B-field components in this order in accordance with cocos 11
end

@compat public Br_Bphi_Bz
push!(document[Symbol("Physics fields")], :Br_Bphi_Bz)

"""
    Bx_By_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, B0::T, R0::T, r::T, z::T) where {T<:Real}

Returns Bx, By, Bz named tuple evaluated at x,y,z starting from ψ interpolant and B0 and R0
"""
function Bx_By_Bz(PSI_interpolant::Interpolations.AbstractInterpolation,
    B0::T, R0::T, x::T, y::T, z::T) where {T<:Real}

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
    Br_Bz_meshgrid(PSI_interpolant::Interpolations.AbstractInterpolation, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
"""
function Br_Bz_meshgrid(PSI_interpolant::Interpolations.AbstractInterpolation, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
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
    Bp(PSI_interpolant::Interpolations.AbstractInterpolation, r::T, z::T) where {T<:Real}
"""
function Bp(PSI_interpolant::Interpolations.AbstractInterpolation, r::T, z::T) where {T<:Real}
    Br, Bz = Br_Bz(PSI_interpolant, r, z)
    return sqrt(Br^2 + Bz^2)
end

"""
    Bp_meshgrid(PSI_interpolant::Interpolations.AbstractInterpolation, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
"""
function Bp_meshgrid(PSI_interpolant::Interpolations.AbstractInterpolation, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
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
    ImplicitMidpointState{T,F<:Function}(vector_field, current_point, step_size, count, Δϕ)

Represents an implicit midpoint update equation for numerical integration.

# Fields

  - `vector_field`: Function representing the vector field to be integrated
  - `current_point`: Current position in the vector field
  - `step_size`: Step size for numerical integration
  - `count`: Counter
  - `Δϕ`: Phi angle traversed
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

Implements the implicit midpoint method update equation.

# Arguments

  - `next_point`: Proposed next point in the integration

# Returns

Residual of the implicit midpoint equation
"""
function (obj::ImplicitMidpointState)(next_point)
    midpoint = 0.5 * (obj.current_point .+ next_point)
    return next_point .- obj.current_point .- obj.step_size * obj.vector_field(midpoint)
end

"""
    _next!(obj::ImplicitMidpointState)

Advances the trajectory by one step using the implicit midpoint method.

# Returns

Next point in the trajectory
"""
function _next!(obj::ImplicitMidpointState)
    # Use root finding to solve the implicit equation
    sol = NLsolve.nlsolve(obj, obj.current_point)
    next_point = sol.zero

    # Update implicit equation
    dphi = abs(acos(dot(next_point[1:2], obj.current_point[1:2]) / (norm(next_point[1:2]) * norm(obj.current_point[1:2]))))
    obj.current_point .= next_point
    obj.count += 1
    obj.Δϕ += dphi

    return next_point
end

"""
    trace_field_line(vector_field, start_point, step_size, stop_condition)

Compute a trajectory along a vector field using the implicit midpoint method.

# Arguments

  - `vector_field`: Vector field function
  - `start_point`: Initial point of the trajectory
  - `step_size`: Size of each integration step
  - `stop_condition`: Function determining when to stop integration. Takes a ImplicitMidpointState struct as an argument.

# Returns

Trajectory represented as a vector of points
"""
function trace_field_line(vector_field, start_point, step_size, stop_condition)

    # Initialize implicit equation struct
    implicit_equation = ImplicitMidpointState(vector_field, start_point, step_size, 0, 0.0)

    # Initialize trajectory
    next_point = copy(start_point)
    trajectory = [next_point]

    while !stop_condition(implicit_equation)
        # Update position
        next_point = _next!(implicit_equation)

        # Add to trajectory
        push!(trajectory, next_point)
    end

    return trajectory
end

"""
    trace_field_line(eqt::IMAS.equilibrium__time_slice, r, z; phi=zero(r), step_size=0.01, max_turns=1, stop_condition=nothing)

Compute a field line trajectory for a specific equilibrium time slice.

# Arguments

  - `eqt`: Equilibrium time slice
  - `r`: Radial coordinate
  - `z`: Vertical coordinate
  - `phi`: Toroidal angle (default: 0)
  - `step_size`: Size of each integration step (default: 0.01)
  - `max_turns`: Maximum number of toroidal turns used in default stop condition (default: 1)
  - `stop_condition`: User defined stop condition that takes a ImplicitMidpointState struct and returns a Boolean (default: nothing)

# Returns

Named tuple with x, y, z, r, and phi coordinates of the trajectory
"""
function trace_field_line(eqt::IMAS.equilibrium__time_slice, r, z;
    phi=zero(r), step_size=0.01, max_turns=1,
    stop_condition::Union{Nothing,Function}=nothing)

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    rg, zg, PSI_interpolant = ψ_interpolant(eqt2d)
    R0 = eqt.global_quantities.vacuum_toroidal_field.r0
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

    # Define vector field function
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
    trajectory = trace_field_line(vector_field, start_point, step_size, stop)

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

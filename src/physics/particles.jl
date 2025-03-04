document[Symbol("Physics particles")] = Symbol[]

#= ======================================= =#
#  Particle tracker in a toroidal geometry  #
#= ======================================= =#

"""
    x::T
    y::T
    z::T
    δvx::T
    δvy::T
    δvz::T

Cartesian coordinate system centered in (R=0, Z=0); R = sqrt(X^2 + Y^2) Z = Z
"""
mutable struct Particle{T<:Real}
    x::T
    y::T
    z::T
    δvx::T
    δvy::T
    δvz::T
end

@compat public Particle
push!(document[Symbol("Physics particles")], :Particle)

@recipe function plot_particles(particles::Vector{Particle{T}}, eqt::IMAS.equilibrium__time_slice) where {T<:Real}
    N = length(particles)
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r = eqt2d.grid.dim1
    z = eqt2d.grid.dim2

    # normalization weight to bring bin value in the unitary range
    plot_norm = 0.5 / (N * (r[2] - r[1]) * (z[2] - z[1]) / eqt.profiles_1d.area[end])

    @series begin
        seriestype --> :histogram2d
        nbins := (range(minimum(r), maximum(r), length(r) - 1), range(minimum(z), maximum(z), length(z) - 1))
        aspect_ratio := :equal
        grid := false
        weights := zeros(N) .+ plot_norm
        Rcoord.(particles), Zcoord.(particles)
    end
end

"""
    Rcoord(p::Particle)

Return R coordinate of a Particle
"""
function Rcoord(p::Particle)
    return sqrt(p.x^2 + p.y^2)
end

@compat public Rcoord
push!(document[Symbol("Physics particles")], :Rcoord)

"""
    Zcoord(p::Particle)

Return Z coordinate of a Particle
"""
function Zcoord(p::Particle)
    return p.z
end

@compat public Particle
push!(document[Symbol("Physics particles")], :Particle)

"""
    define_particles(eqt::IMAS.equilibrium__time_slice, psi::Vector{T}, source_1d::Vector{T}, N::Int) where {T<:Real}

Creates a vector of particles from a 1D source (psi, source_1d) launching N particles.

Returns also a scalar (I_per_trace) which is the intensity per trace.
"""
function define_particles(eqt::IMAS.equilibrium__time_slice, psi::Vector{T}, source_1d::Vector{T}, N::Int) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)

    # in-plasma mask
    r = eqt2d.grid.dim1
    z = eqt2d.grid.dim2
    dr = (maximum(diff(r))) # save grid dimension and carry the info to find_flux
    dz = (maximum(diff(z)))

    R = [rr for rr in r, zz in z] # 2D
    Z = [zz for rr in r, zz in z] # 2D
    rz_lcfs = collect(zip(eqt.boundary.outline.r, eqt.boundary.outline.z))
    mask = [PolygonOps.inpolygon((rr, zz), rz_lcfs) for rr in r, zz in z]

    # 2D source
    tmp = eqt2d.psi .* (mask .== 1) .+ eqt2d.psi[end] .* .*(mask .== 0) # to avoid extrapolation
    source_2d = interp1d(psi, source_1d).(tmp) # intensity/m^3
    source_2d .*= mask .== 1 # set things to zero outsize of lcfs
    source_2d .*= R .* (2π * (r[2] - r[1]) * (z[2] - z[1])) # unit of intensity = volume integral of source_2d

    # cumulative distribution function
    CDF = cumsum(source_2d[:])
    I_per_trace = CDF[end] / N
    CDF .= (CDF .- CDF[1]) ./ (CDF[end] - CDF[1])
    ICDF = interp1d(CDF, Float64.(1:length(CDF)), :linear)

    # particles structures
    particles = Vector{Particle{Float64}}(undef, N)
    for k in eachindex(particles)
        dk = Int(ceil(ICDF(rand())))
        ϕ = rand() * 2π # toroidal angle of position - put to 0 for 2D

        θv = rand() * 2π            # toroidal angle of velocity - put to 0 for 2D
        ϕv = acos(rand() * 2.0 - 1.0) # poloidal angle of velocity 0 <= ϕv <= π; comment for 2D and use ϕv = rand() * 2π

        Rk = R[dk]
        Zk = Z[dk]
        yk, xk = Rk .* sincos(ϕ)
        zk = Zk
        δvyk, δvxk = (sin(ϕv)) .* sincos(θv)
        δvzk = cos(ϕv)
        particles[k] = Particle(xk, yk, zk, δvxk, δvyk, δvzk)
    end

    return (particles=particles, I_per_trace=I_per_trace, dr=dr, dz=dz)
end

@compat public define_particles
push!(document[Symbol("Physics particles")], :define_particles)

"""
    find_flux(particles::Vector{Particle{T}}, I_per_trace::T, rwall::Vector{T}, zwall::Vector{T}, dr::T, dz::T; ns::Int=10, debug::Bool=false) where {T<:Real}

Returns the flux at the wall of the quantity brought by particles, together with a vector of the surface elements on the wall (wall_s)

`I_per_trace` is the intensity per trace

`(rwall, zwall)` are the coordinates of the wall

`dr`, `dz` are the grid sizes used for the 2d source generation

`ns` is the window size
"""
function find_flux(particles::Vector{Particle{T}}, I_per_trace::T, rwall::Vector{T}, zwall::Vector{T}, dr::T, dz::T; ns::Int=10, debug::Bool=false) where {T<:Real}
    @assert is_closed_polygon(rwall, zwall) "Wall polygong must be closed"

    # advance particles until they hit the wall
    for p in particles
        ti = toroidal_intersection(rwall, zwall, p.x, p.y, p.z, p.δvx, p.δvy, p.δvz)

        if ti == Inf
            @show p
            error("Infinity while moving particles")
        end

        p.x += p.δvx * ti
        p.y += p.δvy * ti
        p.z += p.δvz * ti
    end

    # find flux [quantity/s/m²]
    # smooth the load of each particle within a window
    # note: flux is defined at the cells, not at the nodes

    d = sqrt.(diff(rwall) .^ 2.0 .+ diff(zwall) .^ 2.0) # length of each cell
    l = cumsum(d)
    wall_r = (rwall[1:end-1] .+ rwall[2:end]) ./ 2.0
    wall_z = (zwall[1:end-1] .+ zwall[2:end]) ./ 2.0
    wall_s = d .* wall_r .* 2π

    if debug
        @show minimum(d)
        @show length(d)
        @show length(l)
        @show ns
        @show length(particles)
        # @show σw = (l[end] / length(l)) * (2ns + 1.0) / 5.0 / sqrt(2) # old standard deviation of the distribution
        # @show σw = 2*(l[end] / length(l))
        # @show σw = l[1]
        # @show σw = 2*maximum(d)
        # @show σw = (l[end] / length(l))
        # @show σw = minimum(d)
        # @show σw = 0.001
        # @show σw = 2*sqrt(dr^2 + dz^2)
        # @show σw = l[end]/200
        # @show σw = 2*maximum([dr,dz])
        @show 2 * maximum([dr, dz])
        @show 1250 * l[end] / length(particles)
        @show σw = maximum([1250 * l[end] / length(particles), 2 * maximum([dr, dz])]) # standard deviation of the distribution
        # @show length(particles)*σw/l[end]
        # @show length(particles)*σw/minimum(d)
        @show I_per_trace / sqrt(2π) / σw # - W/m
        @show I_per_trace / sqrt(2π) / σw / 2 / π / minimum(wall_r) # - W/m2
    else
        # when not debugging, just define std dev
        # 2*maximum([dr,dz]) ensures gaussian is wider than the resolution of the grid of the 2d source (smooth solution)
        # 1250*l[end]/N ensures enough particle density to have decent statistics

        σw = maximum([1250 * l[end] / length(particles), 2 * maximum([dr, dz])])  # standard deviation of the distribution
    end
    @assert minimum(d) !== 0.0 "Error in mesher_heat_flux: mesh has repeating points, instead of unique ones"
    ns = maximum([ns, Int(ceil(5 * σw / minimum(d)))])

    # check that the window size is smaller than the whole wall
    if ns > length(wall_r) / 2
        # what is the minimum set of ns elements that covers at least 5*σw?
        ns = 10
        counter_max = Int(floor(length(wall_r) / 2) - 1)
        dist_min = 5σw
        for counter in 1:counter_max
            # find minimum distance for a set of succesive ns element in d
            for k in 1:length(d)-ns
                dist = sum(d[k:k+ns])
                if dist < dist_min
                    dist_min = dist
                end
            end
            if dist_min < 5σw
                # need to increase ns
                dist_min = 5σw
                ns = ns + 1
            else
                # ns is enough to cover 5σw
                break
            end
        end
    end

    if debug
        norm_th = sqrt(2π) * σw
        @show ns
        @show minimum(d)
        @show argmin(d)
        @show d[end]
        @show norm_th
        max_norm = 0
    end

    stencil = collect(-ns:ns)

    flux_r = zero(wall_r)
    flux_z = zero(wall_z)
    dwall = zero(wall_r)
    index = zero(stencil)
    ll = zeros(2ns + 1)
    window = zeros(2ns + 1)

    for p in particles
        old_r = Rcoord(p)
        old_z = Zcoord(p)
        @. dwall = (wall_r - old_r)^2 + (wall_z - old_z)^2
        index0 = argmin(dwall)
        index .= mod.(stencil .+ index0 .- 1, length(wall_r)) .+ 1

        p.x += p.δvx
        p.y += p.δvy
        p.z += p.δvz
        new_r = Rcoord(p)
        new_z = Zcoord(p)

        @views cumsum!(ll, d[index])
        ll .-= ll[ns+1]
        window .= exp.(-(1 / 2) .* (ll ./ σw) .^ 2)
        norm = trapz(ll, window) # - m

        if debug && (norm > max_norm)
            # this is to check that the numerical norm is very close to the theoretical norm of the gaussian
            @show norm
            @show abs(norm - norm_th) / norm_th
            max_norm = norm
        end

        window ./= norm #  - 1/m
        unit_vector = sqrt((new_r - old_r)^2 + (new_z - old_z)^2)

        @inbounds for (k, i) in enumerate(index)
            flux_r[i] += @. (new_r - old_r) / unit_vector * window[k] * I_per_trace / 2 / π / wall_r[i] # - W/m2
            flux_z[i] += @. (new_z - old_z) / unit_vector * window[k] * I_per_trace / 2 / π / wall_r[i] # - W/m2
        end
    end

    return (flux_r=flux_r, flux_z=flux_z, wall_s=wall_s)
end

@compat public find_flux
push!(document[Symbol("Physics particles")], :find_flux)

"""
    toroidal_intersection(r1::Real, z1::Real, r2::Real, z2::Real, px::Real, py::Real, pz::Real, vx::Real, vy::Real, vz::Real, v2::Real, vz2::Real)

Compute the time of intersection between a moving particle and a toroidal surface defined by two points `(r1, z1)` and `(r2, z2)`.

Returns the smallest positive time at which the particle intersects the surface segment. Returns `Inf` if no valid intersection occurs.

  - `r1`, `z1`: Coordinates of the first endpoint of the toroidal surface segment.
  - `r2`, `z2`: Coordinates of the second endpoint of the toroidal surface segment.
  - `x`, `py`, `pz`: Current position of the particle in Cartesian coordinates.
  - `vx`, `vy`, `vz`: Velocity components of the particle.
  - `v2`: Squared radial velocity component, `vx^2 + vy^2`.
  - `vz2`: Squared z-velocity component, `vz^2`.
"""
function toroidal_intersection(r1::Real, z1::Real, r2::Real, z2::Real, px::Real, py::Real, pz::Real, vx::Real, vy::Real, vz::Real, v2::Real, vz2::Real)
    t = Inf

    if isapprox(z1, z2)
        # a horizontal segment: time for how long it takes to get to this z, then check if r is between segment bounds
        ti = (vz == 0.0) ? -Inf : (z1 - pz) / vz
        if ti >= 0
            ri = sqrt((px + vx * ti)^2 + (py + vy * ti)^2)
            if ri >= r1 && ri <= r2
                t = ti
            end
        end

    else
        # parameterize parabola intersection with a line
        m = (z2 - z1) / (r2 - r1)
        z0 = z1 - m * r1
        m2 = m^2

        a = 0.0
        b = 0.0
        c = 0.0
        if isfinite(m)
            z_z0 = pz - z0
            a = m2 * v2 - vz2
            b = 2 * (m2 * (vx * px + vy * py) - z_z0 * vz)
            c = m2 * (px^2 + py^2) - z_z0^2
        else
            a = v2
            b = 2 * (vx * px + vy * py)
            c = px^2 + py^2 - r1^2
        end

        t1 = -Inf
        t2 = -Inf
        if a == 0
            b != 0 && (t2 = -c / b)
        else
            nb_2a = -b / (2a)
            b2a_ca = (b^2 - 4 * a * c) / (4 * a^2)
            if b2a_ca > 0
                t_sq = sqrt(b2a_ca)
                t1 = nb_2a - t_sq
                t2 = nb_2a + t_sq
            elseif b2a_ca ≈ 0.0
                t2 = nb_2a
            end
        end

        # check if intersection points are at positive time and fall between segment
        t = Inf
        if t1 >= 0
            zi = vz * t1 + pz
            if zi >= z1 && zi <= z2
                t = t1
            end
        end
        if !isfinite(t) && t2 >= 0
            zi = vz * t2 + pz
            if zi >= z1 && zi <= z2
                t = t2
            end
        end
    end

    return t
end

"""
    toroidal_intersection(wallr::Vector{T}, wallz::vector{T}, px::Real, py::Real, pz::Real, vx::Real, vy::Real, vz::Real) where {T<:Real}

Compute the time of intersection between a moving particle and a wall

Returns the smallest positive time at which the particle intersects the surface segment. Returns `Inf` if no valid intersection occurs.

  - `wallr`: Vector with r coordinates of the wall (must be closed)
  - `wallz`: Vector with z coordinates of the wall (must be closed)
  - `px`, `py`, `pz`: Current position of the particle in Cartesian coordinates.
  - `vx`, `vy`, `vz`: Velocity components of the particle
"""
function toroidal_intersection(wallr::Vector{T}, wallz::Vector{T}, px::Real, py::Real, pz::Real, vx::Real, vy::Real, vz::Real) where {T<:Real}
    Nw = length(wallr)

    v2 = vx^2 + vy^2
    vz2 = vz^2

    a = v2 # V_R^2
    b = 2 * (vx * px + vy * py)
    c = px^2 + py^2 # R^2 coordinate of the particle
    rmin2 = c - b^2 / (4a)
    if 1 + rmin2 ≈ 1.0
        # this is to avoid issues when c - b^2 / (4a) is essentially zero but negative
        rmin = 0.0
    else
        rmin = sqrt(rmin2)
    end

    ti = Inf
    @inbounds for k in eachindex(wallr)
        rw1 = wallr[k]
        zw1 = wallz[k]
        rw2 = k == Nw ? wallr[1] : wallr[k+1]
        zw2 = k == Nw ? wallz[1] : wallz[k+1]

        rw1 == rw2 && zw1 == zw2 && continue

        if isapprox(zw1, zw2)
            if (rw1 > rw2)
                rw1, rw2 = rw2, rw1
                zw1, zw2 = zw2, zw1
            end
        elseif (zw1 > zw2)
            rw1, rw2 = rw2, rw1
            zw1, zw2 = zw2, zw1
        end

        (vz > 0 && pz > zw2) && continue
        (vz < 0 && pz < zw1) && continue
        (rw1 < rmin && rw2 < rmin) && continue

        t = toroidal_intersection(rw1, zw1, rw2, zw2, px, py, pz, vx, vy, vz, v2, vz2)
        if t < ti
            ti = t
        end
    end

    return ti
end

"""
    toroidal_intersection(wallr::Vector{T}, wallz::vector{T}, px::Real, py::Real, pz::Real, vx::Real, vy::Real, vz::Real) where {T<:Real}

Compute the time of intersection between a moving particle and a wall

Returns the smallest positive time at which the particle intersects the surface segment. Returns `Inf` if no valid intersection occurs.

  - `wallr`: Vector with r coordinates of the wall (must be closed)
  - `wallz`: Vector with z coordinates of the wall (must be closed)
  - `px`, `py`, `pz`: Current position of the particle in Cartesian coordinates.
  - `pol_angle`, `tor_angle`: poloidal and toroidal angles of injection
"""
function toroidal_intersection(wallr::Vector{T}, wallz::Vector{T}, px::Real, py::Real, pz::Real, pol_angle::Real, tor_angle::Real) where {T<:Real}
    vx, vy, vz = pol_tor_angles_2_vector(pol_angle, tor_angle)
    return toroidal_intersection(wallr, wallz, px, py, pz, vx, vy, vz)
end

@compat public toroidal_intersection
push!(document[Symbol("Physics particles")], :toroidal_intersection)

"""
    pol_tor_angles_2_vector(pol_angle::T, tor_angle::T) where {T<:Real}

Conversion from toroidal and poloidal angles expressed in radians to unit vector
"""
function pol_tor_angles_2_vector(pol_angle::T, tor_angle::T) where {T<:Real}
    vx = cos(pol_angle - pi) * cos(tor_angle)
    vy = cos(pol_angle - pi) * sin(tor_angle)
    vz = sin(pol_angle - pi)
    return (vx=vx, vy=vy, vz=vz)
end

@compat public pol_tor_angles_2_vector
push!(document[Symbol("Physics particles")], :pol_tor_angles_2_vector)

"""
    pencil_beam(starting_position::Vector{T}, velocity_vector::Vector{T}, time::AbstractVector{Float64})

returns named tuple with (x, y, z, r) of a beam injected at given (px, py, pz) starting_position with velocity (vx, vy, vz)
"""
function pencil_beam(starting_position::Vector{T}, velocity_vector::Vector{T}, time::AbstractVector{Float64}) where {T<:Real}
    x = starting_position[1] .+ velocity_vector[1] .* time
    y = starting_position[2] .+ velocity_vector[2] .* time
    z = starting_position[3] .+ velocity_vector[3] .* time
    r = sqrt.(x .^ 2.0 .+ y .^ 2.0)
    return (x=x, y=y, z=z, r=r)
end

"""
    pencil_beam(starting_position::Vector{T}, pol_angle::T, tor_angle::T, time::AbstractVector{Float64}) where {T<:Real}

returns named tuple with (x, y, z, r) of a beam injected at given (px, py, pz) starting_position with given poloidal and toroidal angles
"""
function pencil_beam(starting_position::Vector{T}, pol_angle::T, tor_angle::T, time::AbstractVector{Float64}) where {T<:Real}
    vx, vy, vz = pol_tor_angles_2_vector(pol_angle, tor_angle)
    return pencil_beam(starting_position, [vx, vy, vz], time)
end

@compat public pencil_beam
push!(document[Symbol("Physics particles")], :pencil_beam)
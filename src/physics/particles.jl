#= ======================================= =#
#  Particle tracker in a toroidal geometry  #
#= ======================================= =#

mutable struct particle{T<:Real}
    x::T
    y::T
    z::T
    δvx::T
    δvy::T
    δvz::T
end

function Rcoord(p::particle)
    return sqrt(p.x^2 + p.y^2)
end

function Zcoord(p::particle)
    return p.z
end


"""
define_particles(dd::IMAS.dd, rho::Vector{T}, source_1d::rho::Vector{T} , N::Int) where {T<:Real}

Creates a vector of particles from a 1D source (rho, source_1d) launching N particles.
Returns also a vector (I_per_trace) of the intensity per trace.

"""

function define_particles(dd::IMAS.dd, psi::Vector{T}, source_1d::Vector{T} , N::Int) where {T<:Real}
        
    eqt = dd.equilibrium.time_slice[]
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)

    # in-plasma mask
    r = eqt2d.grid.dim1
    z = eqt2d.grid.dim2
    R = [rr for rr in r, zz in z] # 2D
    Z = [zz for rr in r, zz in z] # 2D
    rz_lcfs = collect(zip(eqt.boundary.outline.r, eqt.boundary.outline.z))
    mask = [IMAS.PolygonOps.inpolygon((rr, zz), rz_lcfs) for rr in r, zz in z]

    # 2D source
    tmp = eqt2d.psi .* (mask .== 1) .+ eqt2d.psi[end] .* .*(mask .== 0) # to avoid extrapolation
    source_2d = IMAS.interp1d(psi, source_1d).(tmp) # intensity/m^3
    source_2d .*= mask .== 1 # set things to zero outsize of lcfs
    source_2d .*= R .* (2π * (r[2] - r[1]) * (z[2] - z[1])) # unit of intensity = volume integral of source_2d


    # cumulative distribution function
    CDF = cumsum(source_2d[:])
    I_per_trace = CDF[end] / N
    CDF .= (CDF .- CDF[1]) ./ (CDF[end] - CDF[1])
    ICDF = IMAS.interp1d(CDF, Float64.(1:length(CDF)), :linear)

    # particles structures
    particles = Vector{particle{Float64}}(undef, N)
    for k in eachindex(particles)
        dk = Int(ceil(ICDF(rand())))
        ϕ = rand() * 2π
        θv = rand() * 2π
        ϕv = acos(rand() * 2.0 - 1.0)
        Rk = R[dk]
        Zk = Z[dk]
        yk, xk = Rk .* sincos(ϕ)
        zk = Zk
        δvyk, δvxk = (sin(ϕv)) .* sincos(θv)
        δvzk = cos(ϕv)
        particles[k] = particle(xk, yk, zk, δvxk, δvyk, δvzk)
    end

    return particles, I_per_trace

end

"""

find_flux(particles::Vector{particle{T}}, I_per_trace:T, rwall::Vector{T}, zwall::Vector{T}; ns::Int = 10)

Returns the flux at the wall of the quantity brought by particles, together with a vector of the surface elements on the wall (wall_s)

I_per_trace is the intensity per trace
(rwall, zwall) is the coordinate of the wall
ns is the window size

"""
function find_flux(particles::Vector{particle{T}}, I_per_trace::T, rwall::Vector{T}, zwall::Vector{T}; ns::Int = 10) where {T<:Real}
    rz_wall = collect(zip(rwall, zwall))
    Nw = length(rz_wall)

    # advance neutrons until they hit the wall
    for p in particles
        #println(n)
        ti = Inf

        xn = p.x
        yn = p.y
        zn = p.z
        vx = p.δvx
        vy = p.δvy
        vz = p.δvz

        v2 = vx^2 + vy^2
        vz2 = vz^2

        a = v2
        b = 2 * (vx * xn + vy * yn)
        c = xn^2 + yn^2
        rmin = sqrt(c - b^2 / (4a))

        @inbounds for k in eachindex(rz_wall)
            rw1, zw1 = rz_wall[k]
            rw2, zw2 = k == Nw ? rz_wall[1] : rz_wall[k+1]

            rw1 == rw2 && zw1 == zw2 && continue

            if zw1 > zw2
                rw1, rw2 = rw2, rw1
                zw1, zw2 = zw2, zw1
            end

            (vz > 0 && zn > zw2) && continue
            (vz < 0 && zn < zw1) && continue
            (rw1 < rmin && rw2 < rmin) && continue

            t = toroidal_intersection(rw1, zw1, rw2, zw2, xn, yn, zn, vx, vy, vz, v2, vz2)
            if t < ti
                ti = t
            end
        end
        p.x += vx * ti
        p.y += vy * ti
        p.z += vz * ti
    end

    # find flux [quantity/s/m²]
    # smooth the load of each particle within a window
    # note: flux is defined at the cells, not at the nodes
    wall_r = (rwall[1:end-1] .+ rwall[2:end]) ./ 2.0
    wall_z = (zwall[1:end-1] .+ zwall[2:end]) ./ 2.0
    d = sqrt.(IMAS.gradient(rwall) .^ 2.0 .+ IMAS.gradient(zwall) .^ 2.0)
    d = @views (d[1:end-1] .+ d[2:end]) ./ 2.0
    l = cumsum(d)
    wall_s = d .* wall_r .* 2π

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
        window .= exp.(.-(ll ./ (l[end] / length(l)) ./ (2ns + 1.0) .* 5.0) .^ 2)
        window ./= sum(window)
        unit_vector = sqrt((new_r - old_r)^2 + (new_z - old_z)^2)

        @inbounds for (k, i) in enumerate(index)
            flux_r[i] += @. (new_r - old_r) / unit_vector * window[k] * I_per_trace / wall_s[i]
            flux_z[i] += @. (new_z - old_z) / unit_vector * window[k] * I_per_trace / wall_s[i]
        end
    end

    return flux_r, flux_z, wall_s
end


@recipe function plot_particles(particles::Vector{particle{T}}, eqt::IMAS.equilibrium__time_slice) where {T<:Real}
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


function toroidal_intersection(r1::Real, z1::Real, r2::Real, z2::Real, x::Real, y::Real, z::Real, vx::Real, vy::Real, vz::Real, v2::Real, vz2::Real)
    m = (z2 - z1) / (r2 - r1)
    z0 = z1 - m * r1
    m2 = m^2

    a = 0.0
    b = 0.0
    c = 0.0
    if isfinite(m)
        z_z0 = z - z0
        a = m2 * v2 - vz2
        b = 2 * (m2 * (vx * x + vy * y) - z_z0 * vz)
        c = m2 * (x^2 + y^2) - z_z0^2
    else
        a = v2
        b = 2 * (vx * x + vy * y)
        c = x^2 + y^2 - r1^2
    end

    t1 = -Inf
    t2 = -Inf
    if a == 0
        b != 0 && (t2 = -c / b)
    else
        nb_2a = -b / (2a)
        b2a_ca = nb_2a^2 - c / a
        if b2a_ca >= 0
            t_sq = sqrt(b2a_ca)
            t1 = nb_2a - t_sq
            t2 = nb_2a + t_sq
        end
    end

    t = Inf
    if t1 > 0
        zi = vz * t1 + z
        if zi > z1 && zi < z2
            t = t1
        end
    end
    if !isfinite(t) && t2 > 0
        zi = vz * t2 + z
        if zi > z1 && zi < z2
            t = t2
        end
    end

    return t
end
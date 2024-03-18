#= ======================================= =#
#  Particle tracker in a toroidal geometry  #
#= ======================================= =#

mutable struct particle{T<:Real}
    x::T        # cartesian coordinate system centered in (R=0, Z=0); R = sqrt(X^2 + Y^2) Z = Z
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
    define_particles(eqt::IMAS.equilibrium__time_slice, psi::Vector{T}, source_1d::Vector{T} , N::Int) where {T<:Real}

Creates a vector of particles from a 1D source (psi, source_1d) launching N particles.
Returns also a scalar (I_per_trace) which is the intensity per trace.

"""
function define_particles(eqt::IMAS.equilibrium__time_slice, psi::Vector{T}, source_1d::Vector{T} , N::Int) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)

    # in-plasma mask
    r = eqt2d.grid.dim1
    z = eqt2d.grid.dim2
    dr = (maximum(diff(r))) # save grid dimension and carry the info to find_flux
    dz = (maximum(diff(z)))

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
        ϕ = rand() * 2π # toroidal angle of position - put to 0 for 2D 

        θv = rand() * 2π            # toroidal angle of velocity - put to 0 for 2D
        ϕv  = acos(rand() * 2.0 - 1.0) # poloidal angle of velocity 0 <= ϕv <= π; comment for 2D and use ϕv = rand() * 2π

        Rk = R[dk]
        Zk = Z[dk]
        yk, xk = Rk .* sincos(ϕ)
        zk = Zk
        δvyk, δvxk = (sin(ϕv)) .* sincos(θv)
        δvzk = cos(ϕv)
        particles[k] = particle(xk, yk, zk, δvxk, δvyk, δvzk)
    end

    return (particles = particles, I_per_trace = I_per_trace, dr = dr, dz = dz)

end

"""

find_flux(particles::Vector{particle{T}}, I_per_trace:T, rwall::Vector{T}, zwall::Vector{T}, dr::T, dz::T; ns::Int = 10)

Returns the flux at the wall of the quantity brought by particles, together with a vector of the surface elements on the wall (wall_s)

I_per_trace is the intensity per trace
(rwall, zwall) is the coordinate of the wall
dr,dz is the grid size used for the 2d source generation
ns is the window size
if debug = true activates code for debugging

"""
function find_flux(particles::Vector{particle{T}}, I_per_trace::T, rwall::Vector{T}, zwall::Vector{T}, dr::T, dz::T; ns::Int = 10, debug::Bool = false) where {T<:Real}
    rz_wall = collect(zip(rwall, zwall))
    Nw = length(rz_wall)

    # advance particles until they hit the wall
    for p in particles
        
        ti = Inf

        xn = p.x
        yn = p.y
        zn = p.z
        vx = p.δvx
        vy = p.δvy
        vz = p.δvz

        v2 = vx^2 + vy^2
        vz2 = vz^2

        a = v2 # V_R^2
        b = 2 * (vx * xn + vy * yn)
        c = xn^2 + yn^2 # R^2 coordinate of the particle
        rmin = sqrt(c - b^2 / (4a))

        @inbounds for k in eachindex(rz_wall)
            rw1, zw1 = rz_wall[k]
            rw2, zw2 = k == Nw ? rz_wall[1] : rz_wall[k+1]

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

            (vz > 0 && zn > zw2) && continue
            (vz < 0 && zn < zw1) && continue
            (rw1 < rmin && rw2 < rmin) && continue

            t = toroidal_intersection(rw1, zw1, rw2, zw2, xn, yn, zn, vx, vy, vz, v2, vz2)
            if t < ti
                ti = t
            end
        end

        if ti == Inf
            @show p
            error("Infinity while moving particles")
        end
        
        p.x += vx * ti
        p.y += vy * ti
        p.z += vz * ti

    end

    # find flux [quantity/s/m²]
    # smooth the load of each particle within a window
    # note: flux is defined at the cells, not at the nodes

      # Check if (rwall, zwall) is closed, if not add a point at the end equal to the first
    if (rwall[1], zwall[1]) == (rwall[end], zwall[end])
        # wall is closed - length(wall_s) = length(rwall) - 1 
        d = sqrt.(diff(rwall) .^ 2.0 .+ diff(zwall) .^ 2.0) # length of each cell
        l = cumsum(d)

        wall_r = (rwall[1:end-1] .+ rwall[2:end]) ./ 2.0
        wall_z = (zwall[1:end-1] .+ zwall[2:end]) ./ 2.0
        wall_s = d .* wall_r .* 2π

    else
        # wall NOT closed - length(wall_s) = length(rwall)
        error(" Wall mesh does not close on itself ")
    end
    
    if debug
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
        @show 2*maximum([dr,dz])
        @show 1250*l[end]/length(particles)
        @show σw = maximum([1250*l[end]/length(particles), 2*maximum([dr,dz])]) # standard deviation of the distribution
        # @show length(particles)*σw/l[end]
        # @show length(particles)*σw/minimum(d)
        @show I_per_trace/sqrt(2π)/σw # - W/m
        @show I_per_trace/sqrt(2π)/σw/2/π/minimum(wall_r) # - W/m2
    else
        # when not debugging, just define std dev
        # 2*maximum([dr,dz]) ensures gaussian is wider than the resolution of the grid of the 2d source (smooth solution)
        # 1250*l[end]/N ensures enough particle density to have decent statistics

        σw = maximum([1250*l[end]/length(particles), 2*maximum([dr,dz])])  # standard deviation of the distribution
    end

    ns = maximum([ns, Int(ceil(5*σw/minimum(d)))])

    # check that the window size is smaller than the whole wall
    if ns > length(wall_r) / 2
        # what is the minimum set of ns elements that covers at least 5*σw?
        ns = 10
        counter_max = floor(length(wall_r)/2)-1
        counter = 0
        dist_min = 5σw
        while dist_min <= 5*σw && counter < counter_max
            # find minimum distance for a set of succesive NN element in d
            for k in 1:length(d) - ns
                dist = sum(d[k:k + ns])
                if dist < dist_min
                    dist_min = dist
                end
            end
            if dist_min < 5σw
                # need to increase
                dist_min = 5σw
                ns = ns + 1
            else
                # ns is enough to cover 5σw
                dist_min = Inf
            end
            counter = counter + 1
        end        
    end
  
    if debug
        norm_th = sqrt(2π)*σw
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
        window .= exp.(-(1/2) .* (ll ./ σw) .^ 2)
        norm = integrate(ll,window) # - m
        
        if debug && (norm > max_norm)
            # this is to check that the numerical norm is very close to the theoretical norm of the gaussian
            @show norm
            @show abs(norm-norm_th)/norm_th
            max_norm = norm
        end
        
        window ./= norm #  - 1/m 
        unit_vector = sqrt((new_r - old_r)^2 + (new_z - old_z)^2)

        @inbounds for (k, i) in enumerate(index)
            flux_r[i] += @. (new_r - old_r) / unit_vector * window[k] * I_per_trace / 2 / π / wall_r[i] # - W/m2
            flux_z[i] += @. (new_z - old_z) / unit_vector * window[k] * I_per_trace / 2 / π / wall_r[i] # - W/m2
        end
    end
    
    if debug
        power_linear_density = zero(wall_r)
        flux = zero(wall_r)
        power_linear_density[index] .= window .* I_per_trace 
        flux[index] .= window .* I_per_trace ./ 2 ./ π ./ wall_r[index]
    end

    return (flux_r = flux_r, flux_z = flux_z, wall_s = wall_s)
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

    t = Inf
    if isapprox(z1, z2)
        # a horizontal segment: time for how long it takes to get to this z, then check if r is between segment bounds
        ti = (vz == 0.0) ? -Inf : (z1 - z) / vz
        if ti >= 0
            ri = sqrt((x + vx * ti) ^ 2 + (y + vy * ti) ^ 2)
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
            b2a_ca = (b^2 - 4 * a * c) / (4 * a ^ 2)
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
            zi = vz * t1 + z
            if zi >= z1 && zi <= z2
                t = t1
            end
        end
        if !isfinite(t) && t2 >= 0
            zi = vz * t2 + z
            if zi >= z1 && zi <= z2
                t = t2
            end
        end
    end

    return t
end
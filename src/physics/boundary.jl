document[Symbol("Physics boundary")] = Symbol[]

"""
    contour_no_edges(X, Y, F, l, n::Int=1)

Returns the contour of F on the domain X, Y,
but removes all points that are on the edges of the domain.
"""
function contour_no_edges(X, Y, F, l, n::Int=1)
    x, y = Contour.coordinates(Contour.lines(Contour.contour(X, Y, F, l))[n])
    index = (x .!= X[1]) .&& (x .!= X[end]) .&& (y .!= Y[1]) .&& (y .!= Y[end])
    return x[index], y[index]
end

"""
    boundary_shape(;
        a::T,
        eps::T,
        kapu::T,
        kapl::T,
        delu::T,
        dell::T,
        zetaou::T,
        zetaiu::T,
        zetail::T,
        zetaol::T,
        zoffset::T,
        upnull::Bool=false,
        lonull::Bool=false,
        npts::Int=90
    ) where {T<:Real}

Function used to generate boundary shapes based on `T. C. Luce, PPCF, 55 9 (2013)`

  - `a`: minor radius
  - `eps`: aspect ratio
  - `kapu`: upper elongation
  - `lkap`: lower elongation
  - `delu`: upper triangularity
  - `dell`: lower triangularity
  - `zetaou`: upper outer squareness
  - `zetaiu`: upper inner squareness
  - `zetail`: lower inner squareness
  - `zetaol`: lower outer squareness
  - `zoffset`: z-offset
  - `upnull`: toggle upper x-point
  - `lonull`: toggle lower x-point
  - `npts`: number of points (per quadrant)

returns tuple with arrays of (r, z, zref)

    >> boundary_shape(;a=0.608,eps=0.374,kapu=1.920,kapl=1.719,delu=0.769,dell=0.463,zetaou=-0.155,zetaiu=-0.255,zetail=-0.174,zetaol=-0.227,zoffset=0.000,upnull=true,lonull=false)
"""
function boundary_shape(;
    a::T,
    eps::T,
    kapu::T,
    kapl::T,
    delu::T,
    dell::T,
    zetaou::T,
    zetaiu::T,
    zetail::T,
    zetaol::T,
    zoffset::T,
    upnull::Bool=false,
    lonull::Bool=false,
    npts::Int=90
) where {T<:Real}

    amin = a

    i5 = 7
    i95 = 100 - i5

    r1 = zeros(npts)
    r2 = zeros(npts)
    r3 = zeros(npts - 1)
    r4 = zeros(npts - 1)
    z1 = zeros(npts)
    z2 = zeros(npts)
    z3 = zeros(npts - 1)
    z4 = zeros(npts - 1)
    z1ref = zeros(npts)
    z2ref = zeros(npts)
    z3ref = zeros(npts - 1)
    z4ref = zeros(npts - 1)

    for is_upper in (false, true)

        if is_upper
            ukap = kapu
            utri = delu
            uosq = zetaou
            uisq = zetaiu
        else
            ukap = kapl
            utri = dell
            uosq = zetaol
            uisq = zetail
        end

        ang = range(0, 2 * π, (npts * 4 + 1))
        i1 = findall((ang .>= 0 * π / 2.0) .&& (ang .< 1 * π / 2.0))
        i2 = findall((ang .>= 1 * π / 2.0) .&& (ang .< 2 * π / 2.0))
        ang1 = ang[i1]
        ang2 = ang[i2]
        rsr2 = 1.0 / sqrt(2.0)
        cc = 1.0 - rsr2

        if uosq < -1.0 * rsr2
            uosq = -1.0 * rsr2
        end
        if uisq < -1.0 * rsr2
            uisq = -1.0 * rsr2
        end

        n1 = @. -log(2.0) / log(uosq * cc + rsr2)
        r1 .= @. amin * (1.0 / eps - utri) + amin * (1.0 + utri) * cos(ang1)^(2.0 / n1)
        z1 .= @. zoffset + amin * ukap * sin(ang1)^(2.0 / n1)
        z1ref .= @. zoffset + amin * ukap * sin(ang1)
        n2 = @. -log(2.0) / log(uisq * cc + rsr2)
        r2 .= @. amin * (1.0 / eps - utri) - amin * (1.0 - utri) * abs(cos(ang2))^(2.0 / n2)
        z2 .= @. zoffset + amin * ukap * sin(ang2)^(2.0 / n2)
        z2ref .= @. zoffset + amin * ukap * sin(ang2)

        if (upnull && is_upper) || (lonull && !is_upper)
            f = range(0.0, 1.0, 100)[2:end]
            n = range(1.0, 5.0, 100)[2:end]
            h1 = @. 1.0 - (1.0 - uosq) * cc
            h2 = @. 1.0 - (1.0 - uisq) * cc
            a1 = @. amin * (1.0 + utri)
            a2 = @. amin * (1.0 - utri)
            b = @. amin * ukap

            if utri >= 0
                c1 = @. utri - 1.0
            else
                c1 = @. -1.0 / (1.0 + utri)
            end

            if utri >= 0.0
                c2 = @. -1.0 / (1.0 - utri)
            else
                c2 = @. -1.0 * (1.0 + utri)
            end

            y1q1 = zeros((99, 99))
            y2q1 = zeros((99, 99))
            y1q2 = zeros((99, 99))
            y2q2 = zeros((99, 99))

            for i in 1:size(y1q1)[1]
                for j in 1:size(y1q1)[2]
                    y1q1[j, i] = (f[i] + h1 * (1.0 - f[i]))^n[j] + (1.0 - f[i]^n[j]) * h1^n[j] - 1.0
                    y2q1[j, i] = f[i]^(n[j] - 1.0) * (f[i] * (c1 + b / a1) - b / a1) - c1
                    y1q2[j, i] = (f[i] + h2 * (1.0 - f[i]))^n[j] + (1.0 - f[i]^n[j]) * h2^n[j] - 1.0
                    y2q2[j, i] = f[i]^(n[j] - 1.0) * (f[i] * (c2 + b / a2) - b / a2) - c2
                end
            end

            xy1q1 = contour_no_edges(f, n, y1q1, 0.0)
            xy2q1 = contour_no_edges(f, n, y2q1, 0.0)
            xy1q2 = contour_no_edges(f, n, y1q2, 0.0)
            xy2q2 = contour_no_edges(f, n, y2q2, 0.0)

            y1q1sol = interp1d(xy1q1[1], xy1q1[2]).(f)
            y2q1sol = interp1d(xy2q1[1], xy2q1[2]).(f)
            y1q2sol = interp1d(xy1q2[1], xy1q2[2]).(f)
            y2q2sol = interp1d(xy2q2[1], xy2q2[2]).(f)

            maxdiffq1 = maximum(y1q1sol .- y2q1sol)
            mindiffq1 = minimum(y1q1sol .- y2q1sol)
            maxdiffq2 = maximum(y1q2sol .- y2q2sol)
            mindiffq2 = minimum(y1q2sol .- y2q2sol)

            if maxdiffq1 / mindiffq1 < 0.0
                imin = argmin_abs(y1q1sol, y2q1sol)
                fsolq1 = f[imin]
                nsolq1 = y1q1sol[imin]
                gsolq1 = @. (1.0 - fsolq1^nsolq1)^(1.0 / nsolq1)
            else
                if maxdiffq1 > 0
                    y1new = @. (f[i95] + h1 * (1.0 - f[i95]))^y2q1sol[i95] + (1.0 - f[i95]^y2q1sol[i95]) * h1^y2q1sol[i95] - 1.0
                    y2new = @. f[i95]^(y2q1sol[i95] - 1.0) * (f[i95] * (c1 + b / a1) - b / a1) - c1
                    while y1new > y2new
                        h1 = @. h1 - 0.01
                        y1new = @. (f[i95] + h1 * (1.0 - f[i95]))^y2q1sol[i95] + (1.0 - f[i95]^y2q1sol[i95]) * h1^y2q1sol[i95] - 1.0
                    end
                    fsolq1 = f[i95]
                    nsolq1 = y2q1sol[i95]
                    gsolq1 = @. (1.0 - fsolq1^nsolq1)^(1.0 / nsolq1)
                else
                    y1new = @. (f[i5] + h1 * (1.0 - f[i5]))^y2q1sol[i5] + (1.0 - f[i5]^y2q1sol[i5]) * h1^y2q1sol[i5] - 1.0
                    y2new = @. f[i5]^(y2q1sol[i5] - 1.0) * (f[i5] * (c1 + b / a1) - b / a1) - c1
                    while y1new < y2new
                        h1 = @. h1 + 0.01
                        y1new = @. (f[i5] + h1 * (1.0 - f[i5]))^y2q1sol[i5] + (1.0 - f[i5]^y2q1sol[i5]) * h1^y2q1sol[i5] - 1.0
                    end
                    fsolq1 = f[i5]
                    nsolq1 = y2q1sol[i5]
                    gsolq1 = @. (1.0 - fsolq1^nsolq1)^(1.0 / nsolq1)
                end
            end

            alpha1 = @. a1 / (1.0 - fsolq1)
            beta1 = @. b / gsolq1
            y1 = @. beta1 * (1.0 - ((r1 - amin * (1.0 / eps + 1.0)) / alpha1 + 1.0)^nsolq1)^(1.0 / nsolq1)
            z1 .= @. y1 + zoffset

            if maxdiffq2 / mindiffq2 < 0.0
                imin = argmin_abs(y1q2sol, y2q2sol)
                fsolq2 = f[imin]
                nsolq2 = y1q2sol[imin]
                gsolq2 = @. (1.0 - fsolq2^nsolq2)^(1.0 / nsolq2)
            else
                if maxdiffq2 > 0.0
                    y1new = @. (f[i95] + h2 * (1.0 - f[i95]))^y2q2sol[i95] + (1.0 - f[i95]^y2q2sol[i95]) * h2^y2q2sol[i95] - 1.0
                    y2new = @. f[i95]^(y2q2sol[i95] - 1.0) * (f[i95] * (c2 + b / a2) - b / a2) - c2
                    while y1new > y2new
                        h2 = @. h2 - 0.01
                        y1new = @. (f[i95] + h2 * (1.0 - f[i95]))^y2q2sol[i95] + (1.0 - f[i95]^y2q2sol[i95]) * h2^y2q2sol[i95] - 1.0
                    end
                    fsolq2 = f[i95]
                    nsolq2 = y2q2sol[i95]
                    gsolq2 = @. (1.0 - fsolq2^nsolq2)^(1.0 / nsolq2)
                else
                    y1new = @. (f[i5] + h2 * (1.0 - f[i5]))^y2q2sol[i5] + (1.0 - f[i5]^y2q2sol[i5]) * h2^y2q2sol[i5] - 1.0
                    y2new = @. f[i5]^(y2q2sol[i5] - 1.0) * (f[i5] * (c2 + b / a2) - b / a2) - c2
                    while y1new < y2new
                        h2 = @. h2 + 0.01
                        y1new = @. (f[i5] + h2 * (1.0 - f[i5]))^y2q2sol[i5] + (1.0 - f[i5]^y2q2sol[i5]) * h2^y2q2sol[i5] - 1.0
                    end
                    fsolq2 = f[i5]
                    nsolq2 = y2q2sol[i5]
                    gsolq2 = @. (1.0 - fsolq2^nsolq2)^(1.0 / nsolq2)
                end
            end

            alpha2 = @. a2 / (1.0 - fsolq2)
            beta2 = @. b / gsolq2
            y2 = @. beta2 * (1.0 - (1.0 + (amin * (1.0 / eps - 1.0) - r2) / alpha2)^nsolq2)^(1.0 / nsolq2)
            z2 .= @. y2 + zoffset
        end

        if !is_upper
            r4 .= reverse!(r1)[2:end]
            z4 .= -reverse!(z1)[2:end]
            z4ref .= -reverse!(z1ref)[2:end]
            r3 .= reverse!(r2)[2:end]
            z3 .= -reverse!(z2)[2:end]
            z3ref .= -reverse!(z2ref)[2:end]
        end

    end

    r = vcat(r1, r2, r3, r4)
    z = vcat(z1, z2, z3, z4)
    zref = vcat(z1ref, z2ref, z3ref, z4ref)

    return (r=r, z=z, zref=zref)
end

@compat public boundary_shape
push!(document[Symbol("Physics boundary")], :boundary_shape)

"""
    boundary(pc::IMAS.pulse_schedule__position_control{T}, time0::Float64)

return boundary from pulse_schedule.position_control at a given time0
"""
function boundary(pc::IMAS.pulse_schedule__position_control{T}, time0::Float64) where {T<:Real}
    r = T[extrap1d(interp1d_itp(pc.time, pcb.r.reference); first=:constant, last=:constant).(time0) for pcb in pc.boundary_outline]
    z = T[extrap1d(interp1d_itp(pc.time, pcb.z.reference); first=:constant, last=:constant).(time0) for pcb in pc.boundary_outline]
    reorder_flux_surface!(r, z)
    return (r=r, z=z)
end

"""
    boundary(pc::IMAS.pulse_schedule__position_control{T}, time_index::Int)

returns boundary from pulse_schedule.position_control at a given time_index
"""
function boundary(pc::IMAS.pulse_schedule__position_control{T}, time_index::Int) where {T<:Real}
    r = T[pcb.r.reference[time_index] for pcb in pc.boundary_outline]
    z = T[pcb.z.reference[time_index] for pcb in pc.boundary_outline]
    reorder_flux_surface!(r, z)
    return (r=r, z=z)
end

"""
    boundary(pc::IMAS.pulse_schedule__position_control; time0::Float64=global_time(pc))

Beturns r,z vectors from pulse_schedule.position_control.equilibrium__time_slice___boundary__outline
"""
function boundary(pc::IMAS.pulse_schedule__position_control; time0::Float64=global_time(pc))
    return boundary(pc, time0)
end

@compat public boundary
push!(document[Symbol("Physics boundary")], :boundary)

"""
    x_points(x_points::IMAS.IDSvector{<:IMAS.pulse_schedule__position_control__x_point{T}}; time0::Float64=global_time(x_points)) where {T<:Real}

Beturns vector with tuples of R,Z coordinates of x-points in pulse_schedule at time0
"""
function x_points(x_points::IMAS.IDSvector{<:IMAS.pulse_schedule__position_control__x_point{T}}; time0::Float64=global_time(x_points)) where {T<:Real}
    x_points0 = Tuple{T,T}[]
    for x_point in x_points
        Rx = get_time_array(x_point.r, :reference, time0)
        if Rx > 0.0 # discard NaN and points with Rx==0
            Zx = get_time_array(x_point.z, :reference, time0)
            push!(x_points0, (Rx, Zx))
        end
    end
    return x_points0
end

@compat public x_points
push!(document[Symbol("Physics boundary")], :x_points)

"""
    strike_points(strike_points::IMAS.IDSvector{<:IMAS.pulse_schedule__position_control__strike_point{T}}; time0::Float64=global_time(strike_points)) where {T<:Real}

Beturns vector with tuples of R,Z coordinates of x-points in pulse_schedule at time0
"""
function strike_points(strike_points::IMAS.IDSvector{<:IMAS.pulse_schedule__position_control__strike_point{T}}; time0::Float64=global_time(strike_points)) where {T<:Real}
    strike_points0 = Tuple{T,T}[]
    for strike_point in strike_points
        Rs = get_time_array(strike_point.r, :reference, time0)
        if Rs > 0.0 # discard NaN and points with Rx==0
            Zs = get_time_array(strike_point.z, :reference, time0)
            push!(strike_points0, (Rs, Zs))
        end
    end
    return strike_points0
end

@compat public strike_points
push!(document[Symbol("Physics boundary")], :strike_points)

"""
    arc_length(pr::AbstractVector{<:Real}, pz::AbstractVector{<:Real}; include_zero::Bool=true)

Compute the cumulative arc length of a 2D curve defined from the coordinate vectors `pr` and `pz`
"""
function arc_length(pr::AbstractVector{<:Real}, pz::AbstractVector{<:Real}; include_zero::Bool=true)
    n = length(pr)
    if include_zero
        ll = Vector{Float64}(undef, n)
        ll[1] = 0.0
        @inbounds for i in 2:n
            dx = pr[i] - pr[i-1]
            dz = pz[i] - pz[i-1]
            ll[i] = ll[i-1] + sqrt(dx * dx + dz * dz)
        end
    else
        ll = Vector{Float64}(undef, n - 1)
        @inbounds for i in 2:n
            dx = pr[i] - pr[i-1]
            dz = pz[i] - pz[i-1]
            if i == 2
                ll[i-1] = sqrt(dx * dx + dz * dz)
            else
                ll[i-1] = ll[i-2] + sqrt(dx * dx + dz * dz)
            end
        end
    end
    return ll
end

"""
    arc_length(px::AbstractVector{<:Real}, py::AbstractVector{<:Real}, pz::AbstractVector{<:Real}; include_zero::Bool=true)

Compute the cumulative arc length of a 3D curve from the Cartesian coordinate vectors `px`, `py`, and `pz`.
"""
function arc_length(px::AbstractVector{<:Real}, py::AbstractVector{<:Real}, pz::AbstractVector{<:Real}; include_zero::Bool=true)
    n = length(px)
    if include_zero
        ll = Vector{Float64}(undef, n)
        ll[1] = 0.0
        @inbounds for i in 2:n
            dx = px[i] - px[i-1]
            dy = py[i] - py[i-1]
            dz = pz[i] - pz[i-1]
            ll[i] = ll[i-1] + sqrt(dx * dx + dy * dy + dz * dz)
        end
    else
        ll = Vector{Float64}(undef, n - 1)
        @inbounds for i in 2:n
            dx = px[i] - px[i-1]
            dy = py[i] - py[i-1]
            dz = pz[i] - pz[i-1]
            if i == 2
                ll[i-1] = sqrt(dx * dx + dy * dy + dz * dz)
            else
                ll[i-1] = ll[i-2] + sqrt(dx * dx + dy * dy + dz * dz)
            end
        end
    end
    return ll
end

"""
    arc_length_cylindrical(r::AbstractVector{<:Real}, phi::AbstractVector{<:Real}, z::AbstractVector{<:Real}; include_zero::Bool=true)

Compute the cumulative arc length of a 3D curve in cylindrical coordinates using  radial (`r`), angular (`phi`, in radians), and axial (`z`) vectors.

The arc length is calculated as:

    Δs = sqrt((Δr)^2 + (r_avg * Δphi)^2 + (Δz)^2)
"""
function arc_length_cylindrical(r::AbstractVector{<:Real}, phi::AbstractVector{<:Real}, z::AbstractVector{<:Real}; include_zero::Bool=true)
    n = length(r)
    if include_zero
        ll = Vector{Float64}(undef, n)
        ll[1] = 0.0
        @inbounds for i in 2:n
            dr = r[i] - r[i-1]
            dphi = phi[i] - phi[i-1]
            dz = z[i] - z[i-1]
            r_avg = 0.5 * (r[i] + r[i-1])
            arc = sqrt(dr^2 + (r_avg^2 * dphi^2) + dz^2)
            ll[i] = ll[i-1] + arc
        end
    else
        ll = Vector{Float64}(undef, n - 1)
        @inbounds for i in 2:n
            dr = r[i] - r[i-1]
            dphi = phi[i] - phi[i-1]
            dz = z[i] - z[i-1]
            r_avg = 0.5 * (r[i] + r[i-1])
            arc = sqrt(dr^2 + (r_avg^2 * dphi^2) + dz^2)
            if i == 2
                ll[i-1] = arc
            else
                ll[i-1] = ll[i-2] + arc
            end
        end
    end
    return ll
end

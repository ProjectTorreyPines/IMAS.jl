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

  - a: minor radius
  - eps: aspect ratio
  - kapu: upper elongation
  - lkap: lower elongation
  - delu: upper triangularity
  - dell: lower triangularity
  - zetaou: upper outer squareness
  - zetaiu: upper inner squareness
  - zetail: lower inner squareness
  - zetaol: lower outer squareness
  - zoffset: z-offset
  - upnull: toggle upper x-point
  - lonull: toggle lower x-point
  - npts: number of points (per quadrant)

returns tuple with arrays of (r, z, zref)

> > boundary_shape(;a=0.608,eps=0.374,kapu=1.920,kapl=1.719,delu=0.769,dell=0.463,zetaou=-0.155,zetaiu=-0.255,zetail=-0.174,zetaol=-0.227,zoffset=0.000,upnull=true,lonull=false)
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
                imin = argmin(abs.(y1q1sol .- y2q1sol))
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
                imin = argmin(abs.(y1q2sol .- y2q2sol))
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

function boundary(pc::IMAS.pulse_schedule__position_control{T}, time0::Float64) where {T<:Real}
    return (r=T[extrap1d(interp1d_itp(pc.time, pcb.r.reference); first=:flat, last=:flat).(time0) for pcb in pc.boundary_outline],
        z=T[extrap1d(interp1d_itp(pc.time, pcb.z.reference); first=:flat, last=:flat).(time0) for pcb in pc.boundary_outline])
end

function boundary(pc::IMAS.pulse_schedule__position_control{T}, time_index::Int) where {T<:Real}
    return (r=T[pcb.r.reference[time_index] for pcb in pc.boundary_outline],
        z=T[pcb.z.reference[time_index] for pcb in pc.boundary_outline])
end

"""
    boundary(pc::IMAS.pulse_schedule__position_control; time0::Float64=global_time(pc))

Beturns r,z vectors from pulse_schedule.position_control.equilibrium__time_slice___boundary__outline
"""
function boundary(pc::IMAS.pulse_schedule__position_control; time0::Float64=global_time(pc))
    return boundary(pc, time0)
end

"""
    x_points(x_points::IMAS.IDSvector{<:IMAS.pulse_schedule__position_control__x_point{T}}; time0::Float64=global_time(x_points)) where {T<:Real}

Beturns vector with tuples of R,Z coordinates of x-points in pulse_schedule at time0
"""
function x_points(x_points::IMAS.IDSvector{<:IMAS.pulse_schedule__position_control__x_point{T}}; time0::Float64=global_time(x_points)) where {T<:Real}
    x_points0 = Tuple{T,T}[]
    for x_point in x_points
        Rx = get_time_array(x_point.r, :reference, time0)
        Zx = get_time_array(x_point.z, :reference, time0)
        push!(x_points0, (Rx, Zx))
    end
    return x_points0
end

"""
    strike_points(strike_points::IMAS.IDSvector{<:IMAS.pulse_schedule__position_control__strike_point{T}}; time0::Float64=global_time(strike_points)) where {T<:Real}

Beturns vector with tuples of R,Z coordinates of x-points in pulse_schedule at time0
"""
function strike_points(strike_points::IMAS.IDSvector{<:IMAS.pulse_schedule__position_control__strike_point{T}}; time0::Float64=global_time(strike_points)) where {T<:Real}
    strike_points0 = Tuple{T,T}[]
    for x_point in strike_points
        Rxx = get_time_array(x_point.r, :reference, time0)
        Zxx = get_time_array(x_point.z, :reference, time0)
        push!(strike_points0, (Rxx, Zxx))
    end
    return strike_points0
end

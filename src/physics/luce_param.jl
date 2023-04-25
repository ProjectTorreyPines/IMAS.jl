import Contour

function contour_noedges(X, Y, F, l, n=1)
    x, y = Contour.coordinates(Contour.lines(Contour.contour(X, Y, F, l))[n])
    index = (x .!= X[1]) .&& (x .!= X[end]) .&& (y .!= Y[1]) .&& (y .!= Y[end])
    return x[index], y[index]
end

"""
Function used to generate boundary shapes based on `T. C. Luce, PPCF, 55 9 (2013)`
Direct Python translation of the IDL program /u/luce/idl/shapemaker3.pro

:param a: minor radius

:param eps: aspect ratio

:param kapu: upper elongation

:param lkap: lower elongation

:param delu: upper triangularity

:param dell: lower triangularity

:param zetaou: upper outer squareness

:param zetaiu: upper inner squareness

:param zetail: lower inner squareness

:param zetaol: lower outer squareness

:param zoffset: z-offset

:param upnull: toggle upper x-point

:param lonull: toggle lower x-point

:param npts: number of points (per quadrant)

:param newsq: A 4 element array, into which the new squareness values are stored

:return: tuple with arrays of r,z,zref

>> boundaryShape(a=0.608,eps=0.374,kapu=1.920,kapl=1.719,delu=0.769,dell=0.463,zetaou=-0.155,zetaiu=-0.255,zetail=-0.174,zetaol=-0.227,zoffset=0.000,upnull=False,lonull=False,doPlot=True)
"""
function boundaryShape(;
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
    npts::Int=90,
    newsq=zeros(4)
) where {T<:Real}

    ukap = kapu
    lkap = kapl
    utri = delu
    ltri = dell
    uosq = zetaou
    uisq = zetaiu
    lisq = zetail
    losq = zetaol
    amin = a
    newsq[1:4] = [zetaou, uisq, lisq, losq]
    ang = LinRange(0, 2 * π, (npts * 4 + 1))

    i1 = findall((ang .>= 0 * π / 2.0) .&& (ang .< 1 * π / 2.0))
    i2 = findall((ang .>= 1 * π / 2.0) .&& (ang .< 2 * π / 2.0))
    i3 = findall((ang .>= 2 * π / 2.0) .&& (ang .< 3 * π / 2.0))
    i4 = findall((ang .>= 3 * π / 2.0) .&& (ang .< 4 * π / 2.0))

    ang1 = ang[i1]
    ang2 = ang[i2]
    ang3 = ang[i3]
    ang4 = ang[i4]

    rsr2 = 1.0 / sqrt(2.0)
    cc = 1.0 - rsr2

    if uosq < -1.0 * rsr2
        uosq = -1.0 * rsr2
    end

    if uisq < -1.0 * rsr2
        uisq = -1.0 * rsr2
    end

    if lisq < -1.0 * rsr2
        lisq = -1.0 * rsr2
    end

    if losq < -1.0 * rsr2
        losq = -1.0 * rsr2
    end
    n1 = @. -log(2.0) / log(uosq * cc + rsr2)
    r1 = @. amin * (1.0 / eps - utri) + amin * (1.0 + utri) * cos(ang1)^(2.0 / n1)
    z1 = @. zoffset + amin * ukap * sin(ang1)^(2.0 / n1)
    z1ref = @. zoffset + amin * ukap * sin(ang1)
    n2 = @. -log(2.0) / log(uisq * cc + rsr2)
    r2 = @. amin * (1.0 / eps - utri) - amin * (1.0 - utri) * abs(cos(ang2))^(2.0 / n2)
    z2 = @. zoffset + amin * ukap * sin(ang2)^(2.0 / n2)
    z2ref = @. zoffset + amin * ukap * sin(ang2)

    if upnull
        f = LinRange(0.01, 1, 100)
        n = LinRange(1.04, 5, 100)
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
        y1q1 = zeros((100, 100))
        y2q1 = zeros((100, 100))
        y1q2 = zeros((100, 100))
        y2q2 = zeros((100, 100))
        for i in 1:size(y1q1)[1]
            for j in 1:size(y1q1)[2]
                y1q1[j, i] = (f[i] + h1 * (1.0 - f[i]))^n[j] + (1.0 - f[i]^n[j]) * h1^n[j] - 1.0
                y2q1[j, i] = f[i]^(n[j] - 1.0) * (f[i] * (c1 + b / a1) - b / a1) - c1
                y1q2[j, i] = (f[i] + h2 * (1.0 - f[i]))^n[j] + (1.0 - f[i]^n[j]) * h2^n[j] - 1.0
                y2q2[j, i] = f[i]^(n[j] - 1.0) * (f[i] * (c2 + b / a2) - b / a2) - c2
            end
        end

        xy1q1 = contour_noedges(f, n, y1q1, 0.0)
        xy2q1 = contour_noedges(f, n, y2q1, 0.0)
        xy1q2 = contour_noedges(f, n, y1q2, 0.0)
        xy2q2 = contour_noedges(f, n, y2q2, 0.0)

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
                y1new = @. (f[95] + h1 * (1.0 - f[95]))^y2q1sol[95] + (1.0 - f[95]^y2q1sol[95]) * h1^y2q1sol[95] - 1.0
                y2new = @. f[95]^(y2q1sol[95] - 1.0) * (f[95] * (c1 + b / a1) - b / a1) - c1
                while y1new > y2new
                    h1 = @. h1 - 0.01
                    y1new = @. (f[95] + h1 * (1.0 - f[95]))^y2q1sol[95] + (1.0 - f[95]^y2q1sol[95]) * h1^y2q1sol[95] - 1.0
                end
                fsolq1 = f[95]
                nsolq1 = y2q1sol[95]
                gsolq1 = @. (1.0 - fsolq1^nsolq1)^(1.0 / nsolq1)
            else
                y1new = @. (f[5] + h1 * (1.0 - f[5]))^y2q1sol[5] + (1.0 - f[5]^y2q1sol[5]) * h1^y2q1sol[5] - 1.0
                y2new = @. f[5]^(y2q1sol[5] - 1.0) * (f[5] * (c1 + b / a1) - b / a1) - c1
                while y1new < y2new
                    h1 = @. h1 + 0.01
                    y1new = @. (f[5] + h1 * (1.0 - f[5]))^y2q1sol[5] + (1.0 - f[5]^y2q1sol[5]) * h1^y2q1sol[5] - 1.0
                end
                fsolq1 = f[5]
                nsolq1 = y2q1sol[5]
                gsolq1 = @. (1.0 - fsolq1^nsolq1)^(1.0 / nsolq1)
            end
            sqnew1 = @. 1.0 - (1.0 - h1) / cc
            newsq[1] = sqnew1
        end

        alpha1 = @. a1 / (1.0 - fsolq1)
        beta1 = @. b / gsolq1

        y1 = @. beta1 * (1.0 - ((r1 - amin * (1.0 / eps + 1.0)) / alpha1 + 1.0)^nsolq1)^(1.0 / nsolq1)
        z1 = @. y1 + zoffset

        if maxdiffq2 / mindiffq2 < 0.0
            imin = argmin(abs.(y1q2sol .- y2q2sol))
            fsolq2 = f[imin]
            nsolq2 = y1q2sol[imin]
            gsolq2 = @. (1.0 - fsolq2^nsolq2)^(1.0 / nsolq2)

        else
            if maxdiffq2 > 0.0
                y1new = @. (f[95] + h2 * (1.0 - f[95]))^y2q2sol[95] + (1.0 - f[95]^y2q2sol[95]) * h2^y2q2sol[95] - 1.0
                y2new = @. f[95]^(y2q2sol[95] - 1.0) * (f[95] * (c2 + b / a2) - b / a2) - c2
                while y1new > y2new
                    h2 = @. h2 - 0.01
                    y1new = @. (f[95] + h2 * (1.0 - f[95]))^y2q2sol[95] + (1.0 - f[95]^y2q2sol[95]) * h2^y2q2sol[95] - 1.0
                end
                fsolq2 = f[95]
                nsolq2 = y2q2sol[95]
                gsolq2 = @. (1.0 - fsolq2^nsolq2)^(1.0 / nsolq2)

            else
                y1new = @. (f[5] + h2 * (1.0 - f[5]))^y2q2sol[5] + (1.0 - f[5]^y2q2sol[5]) * h2^y2q2sol[5] - 1.0
                y2new = @. f[5]^(y2q2sol[5] - 1.0) * (f[5] * (c2 + b / a2) - b / a2) - c2
                while y1new < y2new
                    h2 = @. h2 + 0.01
                    y1new = @. (f[5] + h2 * (1.0 - f[5]))^y2q2sol[5] + (1.0 - f[5]^y2q2sol[5]) * h2^y2q2sol[5] - 1.0
                end
                fsolq2 = f[5]
                nsolq2 = y2q2sol[5]
                gsolq2 = @. (1.0 - fsolq2^nsolq2)^(1.0 / nsolq2)
            end

            sqnew2 = @. 1.0 - (1.0 - h2) / cc
            newsq[2] = sqnew2
        end

        alpha2 = @. a2 / (1.0 - fsolq2)
        beta2 = @. b / gsolq2
        y2 = @. beta2 * (1.0 - (1.0 + (amin * (1.0 / eps - 1.0) - r2) / alpha2)^nsolq2)^(1.0 / nsolq2)
        z2 = @. y2 + zoffset
    end

    n3 = @. -log(2.0) / log(lisq * cc + rsr2)
    r3 = @. amin * (1.0 / eps - ltri) - amin * (1.0 - ltri) * abs(cos(ang3))^(2.0 / n3)
    z3 = @. zoffset - amin * lkap * abs(sin(ang3))^(2.0 / n3)
    z3ref = @. zoffset + amin * lkap * sin(ang3)
    n4 = @. -log(2.0) / log(losq * cc + rsr2)
    r4 = @. amin * (1.0 / eps - ltri) + amin * (1.0 + ltri) * abs(cos(ang4))^(2.0 / n4)
    z4 = @. zoffset - amin * lkap * abs(sin(ang4))^(2.0 / n4)
    z4ref = @. zoffset + amin * lkap * sin(ang4)

    if lonull
        f = LinRange(0.01, 1, 99)
        n = LinRange(1.04, 5, 99)
        h4 = @. 1.0 - (1.0 - losq) * cc
        h3 = @. 1.0 - (1.0 - lisq) * cc
        a4 = @. amin * (1.0 + ltri)
        a3 = @. amin * (1.0 - ltri)
        b = @. amin * lkap

        if (ltri >= 0.0)
            c4 = @. ltri - 1.0
        else
            c4 = @. -1.0 / (1.0 + ltri)
        end

        if (ltri >= 0.0)
            c3 = @. -1.0 / (1.0 - ltri)
        else
            c3 = @. -1.0 * (1.0 + ltri)
        end

        y1q4 = zeros((99, 99))
        y2q4 = zeros((99, 99))
        y1q3 = zeros((99, 99))
        y2q3 = zeros((99, 99))

        for i in 1:size(y1q4)[1]
            for j in 1:size(y1q4)[2]
                y1q4[j, i] = (f[i] + h4 * (1.0 - f[i]))^n[j] + (1.0 - f[i]^n[j]) * h4^n[j] - 1.0
                y2q4[j, i] = f[i]^(n[j] - 1.0) * (f[i] * (c4 + b / a4) - b / a4) - c4
                y1q3[j, i] = (f[i] + h3 * (1.0 - f[i]))^n[j] + (1.0 - f[i]^n[j]) * h3^n[j] - 1.0
                y2q3[j, i] = f[i]^(n[j] - 1.0) * (f[i] * (c3 + b / a3) - b / a3) - c3
            end
        end

        xy1q4 = contour_noedges(f, n, y1q4, 0.0)
        xy2q4 = contour_noedges(f, n, y2q4, 0.0)
        xy1q3 = contour_noedges(f, n, y1q3, 0.0)
        xy2q3 = contour_noedges(f, n, y2q3, 0.0)

        y1q4sol = interp1d(xy1q4[1], xy1q4[2]).(f)
        y2q4sol = interp1d(xy2q4[1], xy2q4[2]).(f)
        y1q3sol = interp1d(xy1q3[1], xy1q3[2]).(f)
        y2q3sol = interp1d(xy2q3[1], xy2q3[2]).(f)

        maxdiffq4 = maximum(y1q4sol .- y2q4sol)
        mindiffq4 = minimum(y1q4sol .- y2q4sol)
        maxdiffq3 = maximum(y1q3sol .- y2q3sol)
        mindiffq3 = minimum(y1q3sol .- y2q3sol)

        if maxdiffq4 / mindiffq4 < 0.0
            imin = argmin(abs.(y1q4sol .- y2q4sol))
            fsolq4 = f[imin]
            nsolq4 = @. y1q4sol[imin]
            gsolq4 = @. (1.0 - fsolq4^nsolq4)^(1.0 / nsolq4)

        else
            if maxdiffq4 > 0.0
                y1new = @. (f[95] + h4 * (1.0 - f[95]))^y2q4sol[95] + (1.0 - f[95]^y2q4sol[95]) * h4^y2q4sol[95] - 1.0
                y2new = @. f[95]^(y2q4sol[95] - 1.0) * (f[95] * (c4 + b / a4) - b / a4) - c4
                while y1new > y2new
                    h4 = @. h4 - 0.01
                    y1new = @. (f[95] + h4 * (1.0 - f[95]))^y2q4sol[95] + (1.0 - f[95]^y2q4sol[95]) * h4^y2q4sol[95] - 1.0
                end
                fsolq4 = f[95]
                nsolq4 = y2q4sol[95]
                gsolq4 = @. (1.0 - fsolq4^nsolq4)^(1.0 / nsolq4)

            else
                y1new = @. (f[5] + h4 * (1.0 - f[5]))^y2q4sol[5] + (1.0 - f[5]^y2q4sol[5]) * h4^y2q4sol[5] - 1.0
                y2new = @. f[5]^(y2q4sol[5] - 1.0) * (f[5] * (c4 + b / a4) - b / a4) - c4
                while y1new < y2new
                    h4 = @. h4 + 0.01
                    y1new = @. (f[5] + h4 * (1.0 - f[5]))^y2q4sol[5] + (1.0 - f[5]^y2q4sol[5]) * h4^y2q4sol[5] - 1.0
                end
                fsolq4 = f[5]
                nsolq4 = @. y2q4sol[5]
                gsolq4 = @. (1.0 - fsolq4^nsolq4)^(1.0 / nsolq4)
            end

            sqnew4 = @. 1.0 - (1.0 - h4) / cc
            newsq[4] = sqnew4
        end

        alpha4 = @. a4 / (1.0 - fsolq4)
        beta4 = @. b / gsolq4
        y4 = @. -1.0 * (beta4 * (1.0 - ((r4 - amin * (1.0 / eps + 1.0)) / alpha4 + 1.0)^nsolq4)^(1.0 / nsolq4))
        z4 = @. y4 + zoffset
        if maxdiffq3 / mindiffq3 < 0.0
            imin = argmin(abs.(y1q3sol .- y2q3sol))
            fsolq3 = f[imin]
            nsolq3 = y1q3sol[imin]
            gsolq3 = @. (1.0 - fsolq3^nsolq3)^(1.0 / nsolq3)

        else
            if maxdiffq3 > 0.0
                y1new = @. (f[95] + h3 * (1.0 - f[95]))^y2q3sol[95] + (1.0 - f[95]^y2q3sol[95]) * h3^y2q3sol[95] - 1.0
                y2new = @. f[95]^(y2q3sol[95] - 1.0) * (f[95] * (c3 + b / a3) - b / a3) - c3
                while y1new > y2new
                    h3 = @. h3 - 0.01
                    y1new = @. (f[95] + h3 * (1.0 - f[95]))^y2q3sol[95] + (1.0 - f[95]^y2q3sol[95]) * h3^y2q3sol[95] - 1.0
                end
                fsolq3 = f[95]
                nsolq3 = y2q3sol[95]
                gsolq3 = @. (1.0 - fsolq3^nsolq3)^(1.0 / nsolq3)

            else
                y1new = @. (f[5] + h3 * (1.0 - f[5]))^y2q3sol[5] + (1.0 - f[5]^y2q3sol[5]) * h3^y2q3sol[5] - 1.0
                y2new = @. f[5]^(y2q3sol[5] - 1.0) * (f[5] * (c3 + b / a3) - b / a3) - c3
                while y1new < y2new
                    h3 = @. h3 + 0.01
                    y1new = @. (f[5] + h3 * (1.0 - f[5]))^y2q3sol[5] + (1.0 - f[5]^y2q3sol[5]) * h3^y2q3sol[5] - 1.0
                end
                fsolq3 = f[5]
                nsolq3 = y2q3sol[5]
                gsolq3 = @. (1.0 - fsolq3^nsolq3)^(1.0 / nsolq3)

            end
            sqnew3 = @. 1.0 - (1.0 - h3) / cc
            newsq[3] = sqnew3
        end

        alpha3 = @. a3 / (1.0 - fsolq3)
        beta3 = @. b / gsolq3
        y3 = @. -1.0 * (beta3 * (1.0 - (1.0 + (amin * (1.0 / eps - 1.0) - r3) / alpha3)^nsolq3)^(1.0 / nsolq3))
        z3 = @. y3 + zoffset
    end

    r = vcat(r1, r2, r3, r4)
    z = vcat(z1, z2, z3, z4)
    zref = vcat(z1ref, z2ref, z3ref, z4ref)

    return r, z, zref
end


import CoordinateConventions: cocos
import Interpolations
import Contour
import StaticArrays
import PolygonOps
import Optim
import NumericalIntegration: integrate, cumul_integrate

function Br_Bz_interpolant(r::AbstractRange, z::AbstractRange, psi::AbstractMatrix; cocos_number::Int = 11)
    cc = cocos(cocos_number)
    PSI_interpolant = Interpolations.CubicSplineInterpolation((r, z), psi)
    Br_vector_interpolant = (x, y) -> cc.sigma_RpZ * Interpolations.gradient(PSI_interpolant, x, y)[2] / x / (2 * pi)^cc.exp_Bp
    Bz_vector_interpolant = (x, y) -> -cc.sigma_RpZ * Interpolations.gradient(PSI_interpolant, x, y)[1] / x / (2 * pi)^cc.exp_Bp
    Br_Bz_vector_interpolant = (x, y) -> (Br_vector_interpolant(x, y), Bz_vector_interpolant(x, y))
    return Br_Bz_vector_interpolant
end


"""
    flux_surfaces(eq::equilibrium; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities in the equilibrium IDS
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eq::equilibrium; upsample_factor::Int = 1)
    for time_index = 1:length(eq.time_slice)
        flux_surfaces(eq.time_slice[time_index]; upsample_factor)
    end
    return eq
end

"""
    flux_surfaces(eqt::equilibrium__time_slice; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eqt::equilibrium__time_slice; upsample_factor::Int = 1)
    R0 = eqt.boundary.geometric_axis.r
    B0 = eqt.profiles_1d.f[end] / R0
    return flux_surfaces(eqt, B0, R0; upsample_factor)
end

"""
    flux_surfaces(eqt::equilibrium__time_slice, B0::Real, R0::Real; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice, B0 and R0
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eqt::equilibrium__time_slice, B0::Real, R0::Real; upsample_factor::Int = 1)
    cc = cocos(11)

    r_upsampled = r = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length = length(eqt.profiles_2d[1].grid.dim1))
    z_upsampled = z = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length = length(eqt.profiles_2d[1].grid.dim2))
    PSI_interpolant = Interpolations.CubicSplineInterpolation((r, z), eqt.profiles_2d[1].psi)
    PSI_upsampled = eqt.profiles_2d[1].psi

    # upsampling for high-resolution r,z flux surface coordinates
    if upsample_factor > 1
        r_upsampled = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length = length(eqt.profiles_2d[1].grid.dim1) * upsample_factor)
        z_upsampled = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length = length(eqt.profiles_2d[1].grid.dim2) * upsample_factor)
        PSI_upsampled = PSI_interpolant(r_upsampled, z_upsampled)
    end

    # Br and Bz evaluated through spline gradient
    Br_vector_interpolant = (x, y) -> [cc.sigma_RpZ * Interpolations.gradient(PSI_interpolant, x[k], y[k])[2] / x[k] / (2 * pi)^cc.exp_Bp for k = 1:length(x)]
    Bz_vector_interpolant = (x, y) -> [-cc.sigma_RpZ * Interpolations.gradient(PSI_interpolant, x[k], y[k])[1] / x[k] / (2 * pi)^cc.exp_Bp for k = 1:length(x)]

    # find magnetic axis
    res = Optim.optimize(x -> PSI_interpolant(x[1], x[2]), [r[Int(round(length(r) / 2))], z[Int(round(length(z) / 2))]], Optim.Newton(), Optim.Options(g_tol = 1E-8); autodiff = :forward)
    eqt.global_quantities.magnetic_axis.r = res.minimizer[1]
    eqt.global_quantities.magnetic_axis.z = res.minimizer[2]

    for item in [
        :elongation,
        :triangularity_lower,
        :triangularity_upper,
        :r_inboard,
        :r_outboard,
        :q,
        :dvolume_dpsi,
        :j_tor,
        :area,
        :volume,
        :gm1,
        :gm2,
        :gm4,
        :gm5,
        :gm8,
        :gm9,
        :phi,
        :trapped_fraction
    ]
        setproperty!(eqt.profiles_1d, item, zeros(eltype(eqt.profiles_1d.psi), size(eqt.profiles_1d.psi)))
    end

    PR = []
    PZ = []
    LL = []
    FLUXEXPANSION = []
    INT_FLUXEXPANSION_DL = zeros(length(eqt.profiles_1d.psi))
    BPL = zeros(length(eqt.profiles_1d.psi))
    for (k, psi_level0) in reverse(collect(enumerate(eqt.profiles_1d.psi)))

        if k == 1 # on axis flux surface is a synthetic one
            eqt.profiles_1d.elongation[1] = eqt.profiles_1d.elongation[2] - (eqt.profiles_1d.elongation[3] - eqt.profiles_1d.elongation[2])
            eqt.profiles_1d.triangularity_upper[1] = 0.0
            eqt.profiles_1d.triangularity_lower[1] = 0.0

            a = (eqt.profiles_1d.r_outboard[2] - eqt.profiles_1d.r_inboard[2]) / 100.0
            b = eqt.profiles_1d.elongation[1] * a

            t = range(0, 2 * pi, length = 17)
            pr = cos.(t) .* a .+ eqt.global_quantities.magnetic_axis.r
            pz = sin.(t) .* b .+ eqt.global_quantities.magnetic_axis.z

            # Extrema on array indices
            _, imaxr = findmax(pr)
            _, iminr = findmin(pr)
            _, imaxz = findmax(pz)
            _, iminz = findmin(pz)
            r_at_max_z, max_z = pr[imaxz], pz[imaxz]
            r_at_min_z, min_z = pr[iminz], pz[iminz]
            z_at_max_r, max_r = pz[imaxr], pr[imaxr]
            z_at_min_r, min_r = pz[iminr], pr[iminr]

        else  # other flux surfaces
            # trace flux surface
            pr, pz, psi_level = flux_surface(
                r_upsampled,
                z_upsampled,
                PSI_upsampled,
                eqt.profiles_1d.psi,
                eqt.global_quantities.magnetic_axis.r,
                eqt.global_quantities.magnetic_axis.z,
                psi_level0,
                true,
            )
            if length(pr) == 0
                error("Could not trace closed flux surface $k out of $(length(eqt.profiles_1d.psi)) at ψ = $(psi_level)")
            end

            # Extrema on array indices
            _, imaxr = findmax(pr)
            _, iminr = findmin(pr)
            _, imaxz = findmax(pz)
            _, iminz = findmin(pz)
            r_at_max_z, max_z = pr[imaxz], pz[imaxz]
            r_at_min_z, min_z = pr[iminz], pz[iminz]
            z_at_max_r, max_r = pz[imaxr], pr[imaxr]
            z_at_min_r, min_r = pz[iminr], pr[iminr]

            # accurate geometric quantities by finding geometric extrema as optimization problem
            w = 1E-4 # push away from magnetic axis
            function fx(x, psi_level)
                try
                    (PSI_interpolant(x[1], x[2]) - psi_level)^2 - (x[1] - eqt.global_quantities.magnetic_axis.r)^2 * w
                catch
                    return 100
                end
            end
            function fz(x, psi_level)
                try
                    (PSI_interpolant(x[1], x[2]) - psi_level)^2 - (x[2] - eqt.global_quantities.magnetic_axis.z)^2 * w
                catch
                    return 100
                end
            end
            res = Optim.optimize(x -> fx(x, psi_level), [max_r, z_at_max_r], Optim.Newton(), Optim.Options(g_tol = 1E-8); autodiff = :forward)
            (max_r, z_at_max_r) = (res.minimizer[1], res.minimizer[2])
            res = Optim.optimize(x -> fx(x, psi_level), [min_r, z_at_min_r], Optim.Newton(), Optim.Options(g_tol = 1E-8); autodiff = :forward)
            (min_r, z_at_min_r) = (res.minimizer[1], res.minimizer[2])
            if psi_level0 != eqt.profiles_1d.psi[end]
                res = Optim.optimize(x -> fz(x, psi_level), [r_at_max_z, max_z], Optim.Newton(), Optim.Options(g_tol = 1E-8); autodiff = :forward)
                (r_at_max_z, max_z) = (res.minimizer[1], res.minimizer[2])
                res = Optim.optimize(x -> fz(x, psi_level), [r_at_min_z, min_z], Optim.Newton(), Optim.Options(g_tol = 1E-8); autodiff = :forward)
                (r_at_min_z, min_z) = (res.minimizer[1], res.minimizer[2])
            end
            # p = plot(pr, pz, label = "")
            # plot!([max_r], [z_at_max_r], marker = :cicle)
            # plot!([min_r], [z_at_min_r], marker = :cicle)
            # plot!([r_at_max_z], [max_z], marker = :cicle)
            # plot!([r_at_min_z], [min_z], marker = :cicle)
            # display(p)
        end

        # geometric
        a = 0.5 * (max_r - min_r)
        b = 0.5 * (max_z - min_z)
        R = 0.5 * (max_r + min_r)
        eqt.profiles_1d.r_outboard[k] = max_r
        eqt.profiles_1d.r_inboard[k] = min_r
        eqt.profiles_1d.elongation[k] = b / a
        eqt.profiles_1d.triangularity_upper[k] = (R - r_at_max_z) / a
        eqt.profiles_1d.triangularity_lower[k] = (R - r_at_min_z) / a

        # poloidal magnetic field (with sign)
        Br = Br_vector_interpolant(pr, pz)
        Bz = Bz_vector_interpolant(pr, pz)
        Bp2 = Br .^ 2.0 .+ Bz .^ 2.0
        Bp_abs = sqrt.(Bp2)
        Bp = (
            Bp_abs .* cc.sigma_rhotp * cc.sigma_RpZ .*
            sign.((pz .- eqt.global_quantities.magnetic_axis.z) .* Br .- (pr .- eqt.global_quantities.magnetic_axis.r) .* Bz)
        )

        # flux expansion
        ll = cumsum(vcat(0.0, sqrt.(diff(pr) .^ 2 + diff(pz) .^ 2)))
        fluxexpansion = 1.0 ./ Bp_abs
        int_fluxexpansion_dl = integrate(ll, fluxexpansion)
        Bpl = integrate(ll, Bp)

        # save flux surface coordinates for later use
        pushfirst!(PR, pr)
        pushfirst!(PZ, pz)
        pushfirst!(LL, ll)
        pushfirst!(FLUXEXPANSION, fluxexpansion)
        INT_FLUXEXPANSION_DL[k] = int_fluxexpansion_dl
        BPL[k] = Bpl

        # flux-surface averaging function
        function flxAvg(input)
            return integrate(ll, input .* fluxexpansion) / int_fluxexpansion_dl
        end

        # trapped fraction
        Bt = eqt.profiles_1d.f[k] ./ pr
        Btot = sqrt.(Bp2 .+ Bt .^ 2)
        Bmax = maximum(Btot)
        Bratio = Btot ./ Bmax
        avg_Btot = flxAvg(Btot)
        avg_Btot2 = flxAvg(Btot .^ 2)
        hf = flxAvg((1.0 .- sqrt.(1.0 .- Bratio) .* (1.0 .+ Bratio ./ 2.0)) ./ Bratio .^ 2)
        h = avg_Btot / Bmax
        h2 = avg_Btot2 / Bmax^2
        ftu = 1.0 - h2 / (h^2) * (1.0 - sqrt(1.0 - h) * (1.0 + 0.5 * h))
        ftl = 1.0 - h2 * hf
        eqt.profiles_1d.trapped_fraction[k] = 0.75 * ftu + 0.25 * ftl

        # gm1 = <1/R^2>
        eqt.profiles_1d.gm1[k] = flxAvg(1.0 ./ pr .^ 2)

        # gm4 = <1/B^2>
        eqt.profiles_1d.gm4[k] = flxAvg(1.0 ./ Btot .^ 2)

        # gm5 = <B^2>
        eqt.profiles_1d.gm5[k] = avg_Btot2

        # gm8 = <R>
        eqt.profiles_1d.gm8[k] = flxAvg(pr)

        # gm9 = <1/R>
        eqt.profiles_1d.gm9[k] = flxAvg(1.0 ./ pr)

        # j_tor = <j_tor/R> / <1/R>
        eqt.profiles_1d.j_tor[k] =
            (
                -cc.sigma_Bp .* (eqt.profiles_1d.dpressure_dpsi[k] + eqt.profiles_1d.f_df_dpsi[k] * eqt.profiles_1d.gm1[k] / (4 * pi * 1e-7)) *
                (2.0 * pi)^cc.exp_Bp
            ) / eqt.profiles_1d.gm9[k]

        # dvolume_dpsi
        eqt.profiles_1d.dvolume_dpsi[k] = (cc.sigma_rhotp * cc.sigma_Bp * sign(flxAvg(Bp)) * int_fluxexpansion_dl * (2.0 * pi)^(1.0 - cc.exp_Bp))

        # q
        eqt.profiles_1d.q[k] = (
            cc.sigma_rhotp .* cc.sigma_Bp .* eqt.profiles_1d.dvolume_dpsi[k] .* eqt.profiles_1d.f[k] .* eqt.profiles_1d.gm1[k] ./
            ((2.0 * pi)^(2.0 - cc.exp_Bp))
        )

        # quantities calculated on the last closed flux surface
        if k == length(eqt.profiles_1d.psi)
            # ip
            eqt.global_quantities.ip = cc.sigma_rhotp * Bpl / (4e-7 * pi)

            # perimeter
            eqt.global_quantities.length_pol = ll[end]
        end
    end

    function flxAvg(input, ll, fluxexpansion, int_fluxexpansion_dl)
        return integrate(ll, input .* fluxexpansion) / int_fluxexpansion_dl
    end

    function volume_integrate(what)
        return integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what)
    end

    function surface_integrate(what)
        return integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what .* eqt.profiles_1d.gm9) ./ 2pi
    end

    function cumlul_volume_integrate(what)
        return cumul_integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what)
    end

    function cumlul_surface_integrate(what)
        return cumul_integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what .* eqt.profiles_1d.gm9) ./ 2pi
    end

    # integral quantities
    for k = 2:length(eqt.profiles_1d.psi)
        # area
        eqt.profiles_1d.area[k] = integrate(eqt.profiles_1d.psi[1:k], eqt.profiles_1d.dvolume_dpsi[1:k] .* eqt.profiles_1d.gm9[1:k]) ./ 2pi

        # volume
        eqt.profiles_1d.volume[k] = integrate(eqt.profiles_1d.psi[1:k], eqt.profiles_1d.dvolume_dpsi[1:k])

        # phi
        eqt.profiles_1d.phi[k] = (cc.sigma_Bp * cc.sigma_rhotp * integrate(eqt.profiles_1d.psi[1:k], eqt.profiles_1d.q[1:k]) * (2.0 * pi)^(1.0 - cc.exp_Bp))
    end

    R = (eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]) / 2.0
    a = (eqt.profiles_1d.r_outboard[end] - eqt.profiles_1d.r_inboard[end]) / 2.0

    # vacuum magnetic field at the geometric center
    Btvac = B0 * R0 / R

    # average poloidal magnetic field
    Bpave = eqt.global_quantities.ip * (4.0 * pi * 1e-7) / eqt.global_quantities.length_pol

    # li
    Bp2v = integrate(eqt.profiles_1d.psi, BPL * (2.0 * pi)^(1.0 - cc.exp_Bp))
    eqt.global_quantities.li_3 = 2 * Bp2v / R0 / (eqt.global_quantities.ip * (4.0 * pi * 1e-7))^2

    # beta_tor
    eqt.global_quantities.beta_tor = abs(volume_integrate(eqt.profiles_1d.pressure) / (Btvac^2 / 2.0 / 4.0 / pi / 1e-7) / eqt.profiles_1d.volume[end])

    # beta_pol
    eqt.global_quantities.beta_pol = abs(volume_integrate(eqt.profiles_1d.pressure) / eqt.profiles_1d.volume[end] / (Bpave^2 / 2.0 / 4.0 / pi / 1e-7))

    # beta_normal
    ip = eqt.global_quantities.ip / 1e6
    eqt.global_quantities.beta_normal = eqt.global_quantities.beta_tor / abs(ip / a / Btvac) * 100

    # rho_tor_norm
    rho = sqrt.(abs.(eqt.profiles_1d.phi ./ (pi * B0)))
    rho_meters = rho[end]
    eqt.profiles_1d.rho_tor_norm = rho ./ rho_meters

    # phi 2D
    eqt.profiles_2d[1].phi =
        Interpolations.CubicSplineInterpolation(eqt.profiles_1d.psi, eqt.profiles_1d.phi, extrapolation_bc = Interpolations.Line()).(eqt.profiles_2d[1].psi)

    # rho 2D in meters
    RHO = sqrt.(abs.(eqt.profiles_2d[1].phi ./ (pi * B0)))

    # gm2: <∇ρ²/R²>
    if false
        RHO_interpolant = Interpolations.CubicSplineInterpolation((r, z), RHO)
        for k = 1:length(eqt.profiles_1d.psi)
            tmp = [Interpolations.gradient(RHO_interpolant, PR[k][j], PZ[k][j]) for j = 1:length(PR[k])]
            dPHI2 = [j[1] .^ 2.0 .+ j[2] .^ 2.0 for j in tmp]
            eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k] .^ 2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
        end
    else
        dRHOdR, dRHOdZ = gradient(RHO, collect(r), collect(z))
        dPHI2_interpolant = Interpolations.CubicSplineInterpolation((r, z), dRHOdR .^ 2.0 .+ dRHOdZ .^ 2.0)
        for k = 1:length(eqt.profiles_1d.psi)
            dPHI2 = dPHI2_interpolant.(PR[k], PZ[k])
            eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k] .^ 2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
        end
    end

    # fix quantities on axis
    for quantity in [:gm2]
        eqt.profiles_1d.gm2[1] =
            Interpolations.CubicSplineInterpolation(
                eqt.profiles_1d.psi[2:end],
                getproperty(eqt.profiles_1d, quantity)[2:end],
                extrapolation_bc = Interpolations.Line(),
            ).(eqt.profiles_1d.psi[1])
    end

    return eqt
end

"""
    flux_surface(eqt::equilibrium__time_slice, psi_level::Real)

returns r,z coordiates of closed flux surface at given psi_level
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real)
    return flux_surface(eqt, psi_level, true)
end

"""
    flux_surface(eqt::equilibrium__time_slice, psi_level::Real, closed::Union{Nothing,Bool})

returns r,z coordiates of open or closed flux surface at given psi_level
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real, closed::Union{Nothing,Bool})
    dim1 = eqt.profiles_2d[1].grid.dim1
    dim2 = eqt.profiles_2d[1].grid.dim2
    PSI = eqt.profiles_2d[1].psi
    psi = eqt.profiles_1d.psi
    r0 = eqt.global_quantities.magnetic_axis.r
    z0 = eqt.global_quantities.magnetic_axis.z
    flux_surface(dim1, dim2, PSI, psi, r0, z0, psi_level, closed)
end

function flux_surface(
    dim1::Union{AbstractVector,AbstractRange},
    dim2::Union{AbstractVector,AbstractRange},
    PSI::AbstractArray,
    psi::Union{AbstractVector,AbstractRange},
    r0::Real,
    z0::Real,
    psi_level::Real,
    closed::Union{Nothing,Bool},
)
    # handle on axis value as the first flux surface
    if psi_level == psi[1]
        psi_level = psi[2]
        # handle boundary by finding accurate lcfs psi
    elseif psi_level == psi[end]
        psi__boundary_level = find_psi_boundary(dim1, dim2, PSI, psi, r0, z0; raise_error_on_not_open = false)
        if psi__boundary_level !== nothing
            if abs(psi__boundary_level - psi_level) < abs(psi[end] - psi[end-1])
                psi_level = psi__boundary_level
            end
        end
    end

    # contouring routine
    cl = Contour.contour(dim1, dim2, PSI, psi_level)

    prpz = []
    # if no open/closed check, then return all contours
    if closed === nothing
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            push!(prpz, (pr, pz))
        end
        return prpz
        # look for closed flux-surface
    elseif closed
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surface that close and contain magnetic axis
            if (pr[1] == pr[end]) && (pz[1] == pz[end]) && (PolygonOps.inpolygon((r0, z0), collect(zip(pr, pz))) == 1)
                return pr, pz, psi_level
            end
        end
        return [], [], psi_level
        # look for open flux-surfaces
    elseif !closed
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surfaces that close or that do not contain magnetic axis
            if (pr[1] != pr[end]) || (pz[1] != pz[end]) || (PolygonOps.inpolygon((r0, z0), collect(zip(pr, pz))) != 1)
                push!(prpz, (pr, pz))
            end
        end
        return prpz
    end
end

"""
    find_psi_boundary(eqt; precision=1e-6, raise_error_on_not_open=true)

Find psi value of the last closed flux surface
"""
function find_psi_boundary(eqt; precision = 1e-6, raise_error_on_not_open = true)
    dim1 = eqt.profiles_2d[1].grid.dim1
    dim2 = eqt.profiles_2d[1].grid.dim2
    PSI = eqt.profiles_2d[1].psi
    psi = eqt.profiles_1d.psi
    r0 = eqt.global_quantities.magnetic_axis.r
    z0 = eqt.global_quantities.magnetic_axis.z
    find_psi_boundary(dim1, dim2, PSI, psi, r0, z0; precision, raise_error_on_not_open)
end

function find_psi_boundary(dim1, dim2, PSI, psi, r0, z0; precision = 1e-6, raise_error_on_not_open)
    psirange_init = [psi[1] * 0.9 + psi[end] * 0.1, psi[end] + 0.5 * (psi[end] - psi[1])]

    pr, pz = flux_surface(dim1, dim2, PSI, psi, r0, z0, psirange_init[1], true)
    if length(pr) == 0
        error("Flux surface at ψ=$(psirange_init[1]) is not closed")
    end

    pr, pz = flux_surface(dim1, dim2, PSI, psi, r0, z0, psirange_init[end], true)
    if length(pr) > 0
        if raise_error_on_not_open
            error("Flux surface at ψ=$(psirange_init[end]) is not open")
        else
            return nothing
        end
    end

    psirange = deepcopy(psirange_init)
    for k = 1:100
        psimid = (psirange[1] + psirange[end]) / 2.0
        pr, pz = flux_surface(dim1, dim2, PSI, psi, r0, z0, psimid, true)
        # closed flux surface
        if length(pr) > 0
            psirange[1] = psimid
            if (abs(psirange[end] - psirange[1]) / abs(psirange[end] + psirange[1]) / 2.0) < precision
                return psimid
            end
            # open flux surface
        else
            psirange[end] = psimid
        end
    end

    error("Could not find closed boundary between ψ=$(psirange_init[1]) and ψ=$(psirange_init[end])")
end

"""
    build_radii(bd::IMAS.build)

Return list of radii in the build
"""
function build_radii(bd::IMAS.build)
    layers_radii = Real[]
    layer_start = 0
    for l in bd.layer
        push!(layers_radii, layer_start)
        layer_start = layer_start + l.thickness
    end
    push!(layers_radii, layer_start)
    return layers_radii
end

"""
    function get_build(
        bd::IMAS.build;
        type::Union{Nothing,Int} = nothing,
        name::Union{Nothing,String} = nothing,
        identifier::Union{Nothing,UInt,Int} = nothing,
        hfs::Union{Nothing,Int,Array} = nothing,
        return_only_one = true,
        return_index = false,
        raise_error_on_missing = true
    )

Select layer(s) in build based on a series of selection criteria
"""
function get_build(
    bd::IMAS.build;
    type::Union{Nothing,Int} = nothing,
    name::Union{Nothing,String} = nothing,
    identifier::Union{Nothing,UInt,Int} = nothing,
    hfs::Union{Nothing,Int,Array} = nothing,
    return_only_one = true,
    return_index = false,
    raise_error_on_missing = true
)

    if isa(hfs, Int)
        hfs = [hfs]
    end

    valid_layers = []
    for (k, l) in enumerate(bd.layer)
        if (name === nothing || l.name == name) &&
           (type === nothing || l.type == type) &&
           (identifier === nothing || l.identifier == identifier) &&
           (hfs === nothing || l.hfs in hfs)
            if return_index
                push!(valid_layers, k)
            else
                push!(valid_layers, l)
            end
        end
    end
    if length(valid_layers) == 0
        if raise_error_on_missing
            error("Did not find build.layer: name=$name type=$type identifier=$identifier hfs=$hfs")
        else
            return nothing
        end
    end
    if return_only_one
        if length(valid_layers) == 1
            return valid_layers[1]
        else
            error("Found multiple layers that satisfy name:$name type:$type identifier:$identifier hfs:$hfs")
        end
    else
        return valid_layers
    end
end

"""
    structures_mask(bd::IMAS.build; ngrid::Int = 257, border_fraction::Real = 0.1, one_is_for_vacuum::Bool = false)

return rmask, zmask, mask of structures that are not vacuum
"""
function structures_mask(bd::IMAS.build; ngrid::Int = 257, border_fraction::Real = 0.1, one_is_for_vacuum::Bool = false)
    border = maximum(bd.layer[end].outline.r) * border_fraction
    xlim = [0.0, maximum(bd.layer[end].outline.r) + border]
    ylim = [minimum(bd.layer[end].outline.z) - border, maximum(bd.layer[end].outline.z) + border]
    rmask = range(xlim[1], xlim[2], length = ngrid)
    zmask = range(ylim[1], ylim[2], length = ngrid * Int(round((ylim[2] - ylim[1]) / (xlim[2] - xlim[1]))))
    mask = ones(length(rmask), length(zmask))

    valid = true
    for layer in vcat(bd.layer[end], bd.layer)
        if layer.type == -1
            valid = false
        end
        if valid && !ismissing(layer.outline, :r)
            outline = StaticArrays.SVector.(layer.outline.r, layer.outline.z)
            if !ismissing(layer, :material) && layer.material == "vacuum"
                for (kr, rr) in enumerate(rmask)
                    for (kz, zz) in enumerate(zmask)
                        if PolygonOps.inpolygon((rr, zz), outline) == 1
                            mask[kr, kz] = 0.0
                        end
                    end
                end
            else
                for (kr, rr) in enumerate(rmask)
                    for (kz, zz) in enumerate(zmask)
                        if PolygonOps.inpolygon((rr, zz), outline) == 1
                            mask[kr, kz] = 1.0
                        end
                    end
                end
            end
        end
    end
    rlim_oh = IMAS.get_build(bd, type = 1).start_radius
    for (kr, rr) in enumerate(rmask)
        for (kz, zz) in enumerate(zmask)
            if rr < rlim_oh
                mask[kr, kz] = 1.0
            end
        end
    end
    if one_is_for_vacuum
        return rmask, zmask, 1.0 .- mask
    else
        return rmask, zmask, mask
    end
end

function total_pressure_thermal!(core_profiles)
    prof1d = core_profiles.profiles_1d[]
    pressure = prof1d.electrons.density .* prof1d.electrons.temperature
    for ion in prof1d.ion
        pressure += ion.density .* ion.temperature
    end
    return pressure * constants.e
end

function calc_beta_thermal_norm!(summary::IMAS.summary, equilibrium::IMAS.equilibrium, core_profiles::IMAS.core_profiles)
    eqt = equilibrium.time_slice[]
    eq1d = eqt.profiles_1d
    cp1d = core_profiles.profiles_1d[]
    pressure_thermal = cp1d.pressure_thermal
    rho = cp1d.grid.rho_tor_norm
    Bt = @ddtime(equilibrium.vacuum_toroidal_field.b0)
    Ip = eqt.global_quantities.ip
    volume_cp = IMAS.interp(eq1d.rho_tor_norm, eq1d.volume)[rho]

    pressure_thermal_avg = integrate(volume_cp, pressure_thermal) / volume_cp[end]
    beta_tor = 2 * constants.μ_0 * pressure_thermal_avg / Bt^2
    @ddtime (summary.global_quantities.beta_tor.value = beta_tor)
    @ddtime (summary.global_quantities.beta_tor_thermal_norm.value = beta_tor * eqt.boundary.minor_radius * abs(Bt) / abs(Ip / 1e6) * 1.0e2)
end

"""
    ion_element(species::Symbol)

returns a `core_profiles__profiles_1d___ion` structure populated with the element information
"""
function ion_element(species::Symbol)
    ion = IMAS.core_profiles__profiles_1d___ion()
    element = resize!(ion.element, 1)
    if species == :H
        element.z_n = 1
        element.a = 1
    elseif species == :D
        element.z_n = 1
        element.a = 2
    elseif species == :DT
        element.z_n = 1
        element.a = 2.5
    elseif species == :T
        element.z_n = 1
        element.a = 3
    elseif species == :He
        element.z_n = 2
        element.a = 4
    elseif species == :C
        element.z_n = 6
        element.a = 12
    elseif species == :Ne
        element.z_n = 10
        element.a = 20
    else
        error("Element $species is not recognized. Add it to the `IMAS.ion_element()` function.")
    end
    ion.label = String(species)
    return ion
end

"""
    new_source(
        source::IMAS.core_sources__source,
        index::Int,
        name::String,
        rho::Union{AbstractVector,AbstractRange},
        volume::Union{AbstractVector,AbstractRange};
        qe::Union{AbstractVector,Missing} = missing,
        qi::Union{AbstractVector,Missing} = missing,
        se::Union{AbstractVector,Missing} = missing,
        jpar::Union{AbstractVector,Missing} = missing,
        momentum::Union{AbstractVector,Missing} = missing)

Populates the IMAS.core_sources__source with given heating, particle, current, momentun profiles


"""
function new_source(
    source::IMAS.core_sources__source,
    index::Int,
    name::String,
    rho::Union{AbstractVector,AbstractRange},
    volume::Union{AbstractVector,AbstractRange};
    electrons_energy::Union{AbstractVector,Missing} = missing,
    electrons_power_inside::Union{AbstractVector,Missing} = missing,
    total_ion_energy::Union{AbstractVector,Missing} = missing,
    total_ion_power_inside::Union{AbstractVector,Missing} = missing,
    electrons_particles::Union{AbstractVector,Missing} = missing,
    electrons_particles_inside::Union{AbstractVector,Missing} = missing,
    j_parallel::Union{AbstractVector,Missing} = missing,
    current_parallel_inside::Union{AbstractVector,Missing} = missing,
    momentum_tor::Union{AbstractVector,Missing} = missing,
    torque_tor_inside::Union{AbstractVector,Missing} = missing)

    source.identifier.name = name
    source.identifier.index = index
    resize!(source.profiles_1d)
    cs1d = source.profiles_1d[]
    cs1d.grid.rho_tor_norm = rho
    cs1d.grid.volume = volume

    if electrons_energy !== missing
        cs1d.electrons.energy = IMAS.interp(LinRange(0, 1, length(electrons_energy)), electrons_energy)(cs1d.grid.rho_tor_norm)
    end
    if electrons_power_inside !== missing
        cs1d.electrons.power_inside = IMAS.interp(LinRange(0, 1, length(electrons_power_inside)), electrons_power_inside)(cs1d.grid.rho_tor_norm)
    end

    if total_ion_energy !== missing
        cs1d.total_ion_energy = IMAS.interp(LinRange(0, 1, length(total_ion_energy)), total_ion_energy)(cs1d.grid.rho_tor_norm)
    end
    if total_ion_power_inside !== missing
        cs1d.total_ion_power_inside = IMAS.interp(LinRange(0, 1, length(total_ion_power_inside)), total_ion_power_inside)(cs1d.grid.rho_tor_norm)
    end

    if electrons_particles !== missing
        cs1d.electrons.particles = IMAS.interp(LinRange(0, 1, length(electrons_particles)), electrons_particles)(cs1d.grid.rho_tor_norm)
    end
    if electrons_particles_inside !== missing
        cs1d.electrons.particles_inside = IMAS.interp(LinRange(0, 1, length(electrons_particles_inside)), electrons_particles_inside)(cs1d.grid.rho_tor_norm)
    end

    if j_parallel !== missing
        cs1d.j_parallel = IMAS.interp(LinRange(0, 1, length(j_parallel)), j_parallel)(cs1d.grid.rho_tor_norm)
    end
    if current_parallel_inside !== missing
        cs1d.current_parallel_inside = IMAS.interp(LinRange(0, 1, length(current_parallel_inside)), current_parallel_inside)(cs1d.grid.rho_tor_norm)
    end

    if momentum_tor !== missing
        cs1d.momentum_tor = IMAS.interp(LinRange(0, 1, length(momentum_tor)), momentum_tor)(cs1d.grid.rho_tor_norm)
    end
    if torque_tor_inside !== missing
        cs1d.torque_tor_inside = IMAS.interp(LinRange(0, 1, length(torque_tor_inside)), torque_tor_inside)(cs1d.grid.rho_tor_norm)
    end

    return source
end
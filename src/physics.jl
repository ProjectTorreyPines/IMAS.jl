import Equilibrium:cocos
import Interpolations
import Contour
import StaticArrays
import PolygonOps
import Optim
using ForwardDiff

function flux_surfaces(eqt::equilibrium__time_slice, B0::Real, r0::Real)
    cc = cocos(3) # for now hardcoded to 3 because testing for 

    r = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length=length(eqt.profiles_2d[1].grid.dim1))
    z = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length=length(eqt.profiles_2d[1].grid.dim2))
    Br_interpolant = Interpolations.CubicSplineInterpolation((r, z), eqt.profiles_2d[1].b_field_r)
    Bz_interpolant = Interpolations.CubicSplineInterpolation((r, z), eqt.profiles_2d[1].b_field_z)
    PSI_interpolant = Interpolations.CubicSplineInterpolation((r, z), eqt.profiles_2d[1].psi)

    eqt.profiles_1d.elongation = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.triangularity_lower = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.triangularity_upper = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.r_inboard = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.r_outboard = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.q = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.dvolume_dpsi = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.j_tor = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.j_parallel = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.volume = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.gm1 = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.gm9 = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.phi = zero(eqt.profiles_1d.psi)

    for (k, psi_level0) in reverse(collect(enumerate(eqt.profiles_1d.psi)))

        # on axis flux surface is a synthetic one
        if k == 1
            eqt.profiles_1d.elongation[1] = eqt.profiles_1d.elongation[2] - (eqt.profiles_1d.elongation[3] - eqt.profiles_1d.elongation[2])
            eqt.profiles_1d.triangularity_upper[1] = 0.0
            eqt.profiles_1d.triangularity_lower[1] = 0.0
            
            t = range(0, 2 * pi, length=17)
            a = 1E-3
            b = eqt.profiles_1d.elongation[1] * a
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
        
        # other flux surfaces
        else
            # trace flux surface
            pr, pz, psi_level = flux_surface(eqt, psi_level0)
            if length(pr) == 0
                error("Could not trace a closed flux surface at ψₙ of $(k / length(eqt.profiles_1d.psi))")
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

            w = 1E-6

            # find extrema as optimization problem
            fx(x, psi_level) = (PSI_interpolant(x[1], x[2]) - psi_level)^2 - (x[1] - eqt.global_quantities.magnetic_axis.r)^2 * w
            res = Optim.optimize(x -> fx(x, psi_level), x -> ForwardDiff.gradient(xx -> fx(xx, psi_level), x), [max_r, z_at_max_r], inplace=false)
            (max_r, z_at_max_r) = (res.minimizer[1], res.minimizer[2])
            res = Optim.optimize(x -> fx(x, psi_level), x -> ForwardDiff.gradient(xx -> fx(xx, psi_level), x), [min_r, z_at_min_r], inplace=false)
            (min_r, z_at_min_r) = (res.minimizer[1], res.minimizer[2])
            fz(x, psi_level) = (PSI_interpolant(x[1], x[2]) - psi_level)^2 - (x[2] - eqt.global_quantities.magnetic_axis.z)^2 * w
            res = Optim.optimize(x -> fz(x, psi_level), x -> ForwardDiff.gradient(xx -> fz(xx, psi_level), x), [r_at_max_z, max_z], inplace=false)
            (r_at_max_z, max_z) = (res.minimizer[1], res.minimizer[2])
            res = Optim.optimize(x -> fz(x, psi_level), x -> ForwardDiff.gradient(xx -> fz(xx, psi_level), x), [r_at_min_z, min_z], inplace=false)
            (r_at_min_z, min_z) = (res.minimizer[1], res.minimizer[2])
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
        Br = Br_interpolant.(pr, pz)
        Bz = Bz_interpolant.(pr, pz)
        Bp_abs = sqrt.(Br.^2.0 .+ Bz.^2.0)
        Bp = (Bp_abs
        .* cc.sigma_rhotp * cc.sigma_RpZ
        .* sign.((pz .- eqt.global_quantities.magnetic_axis.z) .* Br
              .- (pr .- eqt.global_quantities.magnetic_axis.r) .* Bz))

        # flux expansion
        ll = cumsum(vcat([0], sqrt.(diff(pr).^2 + diff(pz).^2)))
        fluxexpansion = 1.0 ./ Bp_abs
        int_fluxexpansion_dl = trapz(ll, fluxexpansion)

        # flux-surface averaging function
        function flxAvg(input)
            return trapz(ll, input .* fluxexpansion) / int_fluxexpansion_dl
        end

        # gm1 = <1/R^2>
        eqt.profiles_1d.gm1[k] = flxAvg(1.0 ./ pr.^2)

        # gm9 = <1/R>
        eqt.profiles_1d.gm9[k] = flxAvg(1.0 ./ pr)

        # j_tor = <j_tor/R> / <1/R>
        eqt.profiles_1d.j_tor[k] = (
                -cc.sigma_Bp
                .* (eqt.profiles_1d.dpressure_dpsi[k] + eqt.profiles_1d.f_df_dpsi[k] * eqt.profiles_1d.gm9[k] / (4 * pi * 1e-7))
                * (2.0 * pi)^cc.exp_Bp
            )

        # dvolume_dpsi
        eqt.profiles_1d.dvolume_dpsi[k] = (
            cc.sigma_rhotp
            * cc.sigma_Bp
            * sign(flxAvg(Bp))
            * int_fluxexpansion_dl
            * (2.0 * pi)^(1.0 - cc.exp_Bp)
        )

        # q
        eqt.profiles_1d.q[k] = (
            cc.sigma_rhotp
            .* cc.sigma_Bp
            .* eqt.profiles_1d.dvolume_dpsi[k]
            .* eqt.profiles_1d.f[k]
            .* eqt.profiles_1d.gm1[k]
            ./ ((2 * pi)^(2.0 - cc.exp_Bp)))

        # quantities calculated on the last closed flux surface
        if k == length(eqt.profiles_1d.psi)
            # ip
            eqt.global_quantities.ip = cc.sigma_rhotp * trapz(ll, Bp) / (4e-7 * pi)
        end
    end

    # integral quantities are defined later
    for (k, psi_level) in enumerate(eqt.profiles_1d.psi)
        # volume
        eqt.profiles_1d.volume[k] = trapz(eqt.profiles_1d.psi[1:k], eqt.profiles_1d.dvolume_dpsi[1:k])

        # phi
        eqt.profiles_1d.phi[k] = (
            cc.sigma_Bp
            * cc.sigma_rhotp
            * trapz(eqt.profiles_1d.psi[1:k], eqt.profiles_1d.q[1:k])
            * (2.0 * pi)^(1.0 - cc.exp_Bp)
        )
    end

    # rho_tor_norm
    rho = sqrt.(abs.(eqt.profiles_1d.phi ./ (pi * B0)))
    eqt.profiles_1d.rho_tor_norm = rho / rho[end]
end

"""
    flux_surface(eqt, psi_level)

returns r,z coordiates of flux surface at given psi_level
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real)
    if psi_level == eqt.profiles_1d.psi[1]
        psi_level = eqt.profiles_1d.psi[2]
    end
    cl = Contour.contour(eqt.profiles_2d[1].grid.dim1,
                         eqt.profiles_2d[1].grid.dim2,
                         eqt.profiles_2d[1].psi,
                         psi_level)
    for line in Contour.lines(cl)
        pr, pz = Contour.coordinates(line)
        if !((pr[1] == pr[end]) & (pz[1] == pz[end]))
            continue
        end
        polygon = StaticArrays.SVector.(pr, pz)
        if PolygonOps.inpolygon((eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z), polygon) == 1
            return pr, pz, psi_level
        end
    end
    return [], [], psi_level
end

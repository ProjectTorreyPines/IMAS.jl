import Equilibrium:cocos
import Interpolations
import Contour
import StaticArrays
import PolygonOps
import Optim

"""
    flux_surfaces(eq::equilibrium)

update flux surface averaged and geometric quantities in the equilibrium IDS
"""
function flux_surfaces(eq::equilibrium)
    for time_index in 1:length(eq.time_slice)
        flux_surfaces(eq, time_index)
    end
    return eq
end

"""
    flux_surfaces(eq::equilibrium, time_index::Int)

Update flux surface averaged and geometric quantities for a given equilibrium IDS time slice
"""
function flux_surfaces(eq::equilibrium, time_index::Int)
    flux_surfaces(eq.time_slice[time_index],
                  eq.vacuum_toroidal_field.b0[time_index],
                  eq.vacuum_toroidal_field.r0)
    return eq.time_slice[time_index]
end

"""
    flux_surfaces(eqt::equilibrium__time_slice, B0::Real, R0::Real)

update flux surface averaged and geometric quantities for a given equilibrum IDS time slice
"""
function flux_surfaces(eqt::equilibrium__time_slice, B0::Real, R0::Real)
    cc = cocos(3) # for now hardcoded to 3 because testing with Solovev

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
    eqt.profiles_1d.gm2 = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.gm9 = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.phi = zero(eqt.profiles_1d.psi)

    PR = []
    PZ = []
    LL = []
    FLUXEXPANSION = []
    INT_FLUXEXPANSION_DL = []
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
            pr, pz, psi_level = flux_surface(eqt, psi_level0, true)
            if length(pr) == 0
                error("Could not trace $(k) closed flux surface at ψ = $(psi_level)  which is  ψₙ = $(k / length(eqt.profiles_1d.psi))")
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
            w = 1E-6
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
            res = Optim.optimize(x -> fx(x, psi_level), [max_r, z_at_max_r], Optim.Newton(), Optim.Options(g_tol=1E-8); autodiff=:forward)
            (max_r, z_at_max_r) = (res.minimizer[1], res.minimizer[2])
            res = Optim.optimize(x -> fx(x, psi_level), [min_r, z_at_min_r], Optim.Newton(), Optim.Options(g_tol=1E-8); autodiff=:forward)
            (min_r, z_at_min_r) = (res.minimizer[1], res.minimizer[2])
            res = Optim.optimize(x -> fz(x, psi_level), [r_at_max_z, max_z], Optim.Newton(), Optim.Options(g_tol=1E-8); autodiff=:forward)
            (r_at_max_z, max_z) = (res.minimizer[1], res.minimizer[2])
            res = Optim.optimize(x -> fz(x, psi_level), [r_at_min_z, min_z], Optim.Newton(), Optim.Options(g_tol=1E-8); autodiff=:forward)
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

        # save flux surface coordinates for later use
        pushfirst!(PR, pr)
        pushfirst!(PZ, pz)
        pushfirst!(LL, ll)
        pushfirst!(FLUXEXPANSION, fluxexpansion)
        pushfirst!(INT_FLUXEXPANSION_DL, int_fluxexpansion_dl)

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

            # perimeter
            eqt.global_quantities.length_pol = ll[end]
        end

    end

    # integral quantities are defined later
    for k in 1:length(eqt.profiles_1d.psi)
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

    R = (eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]) / 2.0
    a = (eqt.profiles_1d.r_outboard[end] - eqt.profiles_1d.r_inboard[end]) / 2.0

    # vacuum magnetic field at the geometric center
    Btvac = B0 * R0 / R

    # average poloidal magnetic field
    Bpave = eqt.global_quantities.ip * (4.0 * pi * 1e-7) / eqt.global_quantities.length_pol

    # beta_tor
    eqt.global_quantities.beta_tor = trapz(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* eqt.profiles_1d.pressure) / (Btvac^2 / 2.0 / 4.0 / pi / 1e-7) / eqt.profiles_1d.volume[end]

    # beta_pol
    eqt.global_quantities.beta_pol = trapz(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* eqt.profiles_1d.pressure) / (Bpave^2 / 2.0 / 4.0 / pi / 1e-7) / eqt.profiles_1d.volume[end]

    # beta_normal
    ip = eqt.global_quantities.ip / 1e6
    eqt.global_quantities.beta_normal = eqt.global_quantities.beta_tor / abs(ip / a / Btvac) * 100

    # rho_tor_norm
    rho = sqrt.(abs.(eqt.profiles_1d.phi ./ (pi * B0)))
    eqt.profiles_1d.rho_tor_norm = rho / rho[end]

    # phi
    eqt.profiles_2d[1].phi = Interpolations.CubicSplineInterpolation(eqt.profiles_1d.psi, eqt.profiles_1d.phi, extrapolation_bc=Interpolations.Line()).(eqt.profiles_2d[1].psi)

    function flxAvg(input, ll, fluxexpansion, int_fluxexpansion_dl)
        return trapz(ll, input .* fluxexpansion) / int_fluxexpansion_dl
    end

    # gm2: <∇ρ²/R²>
    RHO = Interpolations.CubicSplineInterpolation(eqt.profiles_1d.psi, rho, extrapolation_bc=Interpolations.Line()).(eqt.profiles_2d[1].psi)
    RHO_interpolant = Interpolations.CubicSplineInterpolation((r, z), RHO)
    for k in 1:length(eqt.profiles_1d.psi)
        tmp = [Interpolations.gradient(RHO_interpolant, PR[k][j], PZ[k][j]) for j in 1:length(PR[k])]
        dPHI2 = [j[1].^2.0 .+ j[2].^2.0 for j in tmp]
        eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k].^2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
    end

    return eqt
end

"""
    flux_surface(eqt, psi_level)

returns r,z coordiates of flux surface at given psi_level
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real)
    return flux_surface(eqt, psi_level, true)
end

function flux_surface(eqt::equilibrium__time_slice, psi_level::Real, closed::Bool)

    # handle on axis value as the first flux surface
    if psi_level == eqt.profiles_1d.psi[1]
        psi_level = eqt.profiles_1d.psi[2]
    # handle boundary by finding accurate lcfs psi
    elseif psi_level == eqt.profiles_1d.psi[end]
        psi_level = find_psi_boundary(eqt)
    end

    cl = Contour.contour(eqt.profiles_2d[1].grid.dim1,
                         eqt.profiles_2d[1].grid.dim2,
                         eqt.profiles_2d[1].psi,
                         psi_level)

    if ! closed
        prpz = []
        for line in Contour.lines(cl)
            push!(prpz, Contour.coordinates(line))
        end
        return prpz
    end

    for line in Contour.lines(cl)
        pr, pz = Contour.coordinates(line)
        # ignore flux surfaces that do not close
        if closed && ((pr[1] != pr[end]) || (pz[1] != pz[end]))
            continue
        end 
        # only consider flux surfaces that contain magnetic axis
        if (! closed) || (PolygonOps.inpolygon((eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z), StaticArrays.SVector.(pr, pz)) == 1)
            return pr, pz, psi_level
        end
    end
    return [], [], psi_level
end

"""
    function find_psi_boundary(eqt; precision=1e-3)

Find psi value of the last closed flux surface
"""
function find_psi_boundary(eqt; precision=1e-3)
    psirange = [eqt.global_quantities.psi_axis, eqt.global_quantities.psi_boundary + 0.5 * (eqt.global_quantities.psi_boundary - eqt.global_quantities.psi_axis)]
    for k in 1:100
        psimid = (psirange[1] + psirange[end]) / 2.0
        pr, pz = IMAS.flux_surface(eqt, psimid, true)
        if length(pr) > 0
            psirange[1] = psimid
            if abs(psirange[end] - psirange[1]) < precision
                return psimid
            end
        else
            psirange[end] = psimid
        end
    end
end

"""
    radial_build_start_end_radii(rb::IMAS.radial_build, layer_id::Union{Int,String})

Return start and end radii of center stack layer given its name
"""
function radial_build_start_end_radii(rb::IMAS.radial_build, layer_id::Union{Int,String})
    layer_start = 0
    layer_end = nothing
    for l in rb.center_stack
        if isa(layer_id, String) && l.name == layer_id
            layer_end = layer_start + l.thickness
            break
        elseif isa(layer_id, Int) && l.index == layer_id
            layer_end = layer_start + l.thickness
            break
        end
        layer_start = layer_start + l.thickness
    end
    if layer_end === nothing
        error("Did not find radial_build.center_stack layer $layer_id")
    end
    return layer_start, layer_end
end


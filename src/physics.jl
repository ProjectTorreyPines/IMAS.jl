import Equilibrium: cocos
import Interpolations
import Contour
import StaticArrays
import PolygonOps
import Optim

function flux_surfaces(eqt::equilibrium__time_slice)
    cc = cocos(3) # for now hardcoded to 3 because testing for 

    Br_interpolant = Interpolations.interpolate((eqt.profiles_2d[1].grid.dim1, eqt.profiles_2d[1].grid.dim2), eqt.profiles_2d[1].b_field_r, Interpolations.Gridded(Interpolations.Linear()))
    Bz_interpolant = Interpolations.interpolate((eqt.profiles_2d[1].grid.dim1, eqt.profiles_2d[1].grid.dim2), eqt.profiles_2d[1].b_field_z, Interpolations.Gridded(Interpolations.Linear()))

    eqt.profiles_1d.elongation = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.triangularity_lower = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.triangularity_upper = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.r_inboard = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.r_outboard = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.q = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.dvolume_dpsi = zero(eqt.profiles_1d.psi)
    for (k, psi_level) in enumerate(eqt.profiles_1d.psi)
        # trace flux surface
        pr, pz = flux_surface(eqt, psi_level)

        # geometry
        tmp = flux_geo(pr, pz)
        eqt.profiles_1d.elongation[k] = tmp["elongation"]
        eqt.profiles_1d.triangularity_upper[k] = tmp["triangularity_upper"]
        eqt.profiles_1d.triangularity_lower[k] = tmp["triangularity_lower"]
        eqt.profiles_1d.r_outboard[k] = tmp["r_outboard"]
        eqt.profiles_1d.r_inboard[k] = tmp["r_inboard"]
        eqt.boundary.elongation_upper = tmp["elongation_upper"]
        eqt.boundary.elongation_lower = tmp["elongation_lower"]

        dl = tmp["dl"]

        Br = Br_interpolant(pr, pz)
        Bz = Bz_interpolant(pr, pz)

        Bp_abs = sqrt.(Br.^2.0+Bz.^2.0)

        Bp = (Bp_abs
        .* cc.sigma_rhotp * cc.sigma_RpZ
        .* sign.((pz .- eqt.global_quantities.magnetic_axis.z) .* Br
              .- (pr .- eqt.global_quantities.magnetic_axis.r) .* Bz))

        fluxexpansion_dl = dl ./ Bp_abs
        int_fluxexpansion_dl = sum(fluxexpansion_dl)

        # flux-surface averaging function
        function flxAvg(input)
            return sum(fluxexpansion_dl * input) / int_fluxexpansion_dl
        end

        avg=Dict()
        avg["1/R^2"] = flxAvg(1.0./pr.^2)
        avg["Bp"] = flxAvg(Bp)

        eqt.profiles_1d.dvolume_dpsi[k] = (
            cc.sigma_rhotp
            * cc.sigma_Bp
            * sign(avg["Bp"])
            * int_fluxexpansion_dl
            * (2.0 * pi) ^ (1.0 - cc.exp_Bp)
        )

        # eqt.profiles_1d.q[k]=(
        #     cc.sigma_rhotp
        #     *cc.sigma_Bp
        #     *eqt.profiles_1d.dvolume_dpsi[k]
        #     *eqt.profiles_1d.f[k]
        #     *avg["1/R^2"]
        #     / ((2 * pi) ^ (2.0 - cc.exp_Bp))

        eqt.profiles_1d.q[k] = eqt.global_quantities.ip = cc.sigma_rhotp * sum(dl .* Bp) / (4e-7 * pi)

    end
    # special handling of on-axis
    eqt.profiles_1d.triangularity_upper[1] = 0.0
    eqt.profiles_1d.triangularity_lower[1] = 0.0
    eqt.profiles_1d.elongation[1] = eqt.profiles_1d.elongation[2]
end

"""
    flux_surface(eqt, psi_level)

returns r,z coordiates of flux surface at given psi_level
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real)
    if psi_level == eqt.profiles_1d.psi[1]
        psi_level = eqt.profiles_1d.psi[1] * 0.9 + eqt.profiles_1d.psi[2] * 0.1
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
            return pr, pz
        end
    end
end


"""
    fluxGeo(inputR::Vector{Real}, inputZ::Vector{Real})::Dict

Recturns dictionary with geometric properties of a given flux surface
"""
function flux_geo(inputR::Vector{T} where T <: Real, inputZ::Vector{T} where T <: Real)::Dict

    # concatenate inputs 3 times to avoid bound errors in minimization
    if inputR[1] == inputR[2]
        inputR = inputR[1:end - 1]
        inputR = inputZ[1:end - 1]
    end
    inputRclose = vcat(inputR, inputR, inputR)
    inputZclose = vcat(inputZ, inputZ, inputZ)

    # These are the extrema indices
    _, imaxr = findmax(inputR)
    _, iminr = findmin(inputR)
    _, imaxz = findmax(inputZ)
    _, iminz = findmin(inputZ)

    t = 1:length(inputRclose)
    
    interpR = Interpolations.CubicSplineInterpolation(t, inputRclose)
    optR = Optim.optimize(x -> -interpR(x), length(inputR) + imaxr - 1, length(inputR) + imaxr + 1)
    tmaxr = Optim.minimizer(optR)
    optR = Optim.optimize(x -> interpR(x), length(inputR) + iminr - 1, length(inputR) + iminr + 1)
    tminr = Optim.minimizer(optR)

    interpZ = Interpolations.CubicSplineInterpolation(t, inputZclose)
    optZ = Optim.optimize(x -> -interpZ(x), length(inputZ) + imaxz - 1, length(inputZ) + imaxz + 1)
    tmaxz = Optim.minimizer(optZ)
    optZ = Optim.optimize(x -> interpZ(x), length(inputZ) + iminz - 1, length(inputZ) + iminz + 1)
    tminz = Optim.minimizer(optZ)

    r_at_max_z, max_z = interpR(tmaxz), interpZ(tmaxz)
    r_at_min_z, min_z = interpR(tminz), interpZ(tminz)
    z_at_max_r, max_r = interpZ(tmaxr), interpR(tmaxr)
    z_at_min_r, min_r = interpZ(tminr), interpR(tminr)

    dl = vcat([0], sqrt.(diff(inputR).^2 + diff(inputZ).^2))

    geo = Dict()
    geo["dl"] = dl
    geo["r_at_max_z"] = r_at_max_z
    geo["r_at_min_z"] = r_at_min_z
    geo["z_at_max_r"] = z_at_max_r
    geo["z_at_min_r"] = z_at_min_r
    geo["max_z"] = max_z
    geo["min_z"] = min_z
    geo["r_outboard"] = max_r
    geo["r_inboard"] = min_r
    geo["R"] = 0.5 * (max_r + min_r)
    geo["Z"] = 0.5 * (max_z + min_z)
    geo["a"] = 0.5 * (max_r - min_r)
    geo["eps"] = geo["a"] / geo["R"]
    geo["per"] = sum(dl)
    geo["surfArea"] = 2 * pi * sum(inputR .* dl)
    geo["elongation"] = 0.5 * ((max_z - min_z) / geo["a"])
    geo["elongation_upper"] = (max_z - z_at_max_r) / geo["a"]
    geo["elongation_lower"] = (z_at_max_r - min_z) / geo["a"]
    geo["triangularity_upper"] = (geo["R"] - r_at_max_z) / geo["a"]
    geo["triangularity_lower"] = (geo["R"] - r_at_min_z) / geo["a"]
    geo["triangularity"] = 0.5 * (geo["triangularity_upper"] + geo["triangularity_lower"])
    geo["zoffset"] = z_at_max_r

    # NOTE: lonull, upnull, squareness, centroid, zeta have not been translated from OMFIT fluxGeo

    return geo
end
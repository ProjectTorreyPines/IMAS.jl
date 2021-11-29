import Equilibrium:cocos
import Interpolations
import Contour
import StaticArrays
import PolygonOps
import Optim

"""
    flux_surfaces(eq::equilibrium; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities in the equilibrium IDS
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eq::equilibrium; upsample_factor::Int=1)
    for time_index in 1:length(eq.time_slice)
        flux_surfaces(eq.time_slice[time_index]; upsample_factor)
    end
    return eq
end

"""
    flux_surfaces(eqt::equilibrium__time_slice; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eqt::equilibrium__time_slice; upsample_factor::Int=1)
    R0 = eqt.boundary.geometric_axis.r
    B0 = eqt.profiles_1d.f[end] / R0
    flux_surfaces(eqt, B0, R0; upsample_factor)
end

"""
    flux_surfaces(eqt::equilibrium__time_slice, B0::Real, R0::Real; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice, B0 and R0
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eqt::equilibrium__time_slice, B0::Real, R0::Real; upsample_factor::Int=1)
    cc = cocos(11)

    r_upsampled = r = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length=length(eqt.profiles_2d[1].grid.dim1))
    z_upsampled = z = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length=length(eqt.profiles_2d[1].grid.dim2))
    PSI_interpolant = Interpolations.CubicSplineInterpolation((r, z), eqt.profiles_2d[1].psi)
    PSI_upsampled = eqt.profiles_2d[1].psi

    # upsampling for high-resolution r,z flux surface coordinates
    if upsample_factor>1
        r_upsampled = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length=length(eqt.profiles_2d[1].grid.dim1)*upsample_factor)
        z_upsampled = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length=length(eqt.profiles_2d[1].grid.dim2)*upsample_factor)
        PSI_upsampled = PSI_interpolant(r_upsampled,z_upsampled)
    end

    # Br and Bz evaluated through spline gradient
    Br_vector_interpolant = (x,y) -> [cc.sigma_RpZ*Interpolations.gradient(PSI_interpolant, x[k], y[k])[2]/x[k]/(2*pi)^cc.exp_Bp for k in 1:length(x)]
    Bz_vector_interpolant = (x,y) -> [-cc.sigma_RpZ*Interpolations.gradient(PSI_interpolant, x[k], y[k])[1]/x[k]/(2*pi)^cc.exp_Bp for k in 1:length(x)]

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

            a = (eqt.profiles_1d.r_outboard[2] - eqt.profiles_1d.r_inboard[2]) / 100.0
            b = eqt.profiles_1d.elongation[1] * a

            t = range(0, 2 * pi, length=17)
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
            pr, pz, psi_level = flux_surface(r_upsampled, z_upsampled, PSI_upsampled, eqt.profiles_1d.psi, eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z, psi_level0, true)
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
        Br = Br_vector_interpolant(pr, pz)
        Bz = Bz_vector_interpolant(pr, pz)
        Bp_abs = sqrt.(Br.^2.0 .+ Bz.^2.0)
        Bp = (Bp_abs
            .* cc.sigma_rhotp * cc.sigma_RpZ
            .* sign.((pz .- eqt.global_quantities.magnetic_axis.z) .* Br
                  .- (pr .- eqt.global_quantities.magnetic_axis.r) .* Bz))

        # flux expansion
        ll = cumsum(vcat(0.0, sqrt.(diff(pr).^2 + diff(pz).^2)))
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
    rho_meters = rho[end]
    eqt.profiles_1d.rho_tor_norm = rho ./ rho_meters

    # phi 2D
    eqt.profiles_2d[1].phi = Interpolations.CubicSplineInterpolation(eqt.profiles_1d.psi, eqt.profiles_1d.phi, extrapolation_bc=Interpolations.Line()).(eqt.profiles_2d[1].psi)

    function flxAvg(input, ll, fluxexpansion, int_fluxexpansion_dl)
        return trapz(ll, input .* fluxexpansion) / int_fluxexpansion_dl
    end

    # rho 2D in meters
    RHO = sqrt.(abs.(eqt.profiles_2d[1].phi ./ (pi * B0)))

    # gm2: <∇ρ²/R²>
    if true
        RHO_interpolant = Interpolations.CubicSplineInterpolation((r, z), RHO)
        for k in 1:length(eqt.profiles_1d.psi)
            tmp = [Interpolations.gradient(RHO_interpolant, PR[k][j], PZ[k][j]) for j in 1:length(PR[k])]
            dPHI2 = [j[1].^2.0 .+ j[2].^2.0 for j in tmp]
            eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k].^2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
        end
    else
        dRHOdR,dRHOdZ = gradient(RHO,collect(r),collect(z))
        dPHI2_interpolant = Interpolations.CubicSplineInterpolation((r, z), dRHOdR.^2.0.+dRHOdZ.^2.0)
        for k in 1:length(eqt.profiles_1d.psi)
            dPHI2 = dPHI2_interpolant.(PR[k],PZ[k])
            eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k].^2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
        end
    end

    # fix quantities on axis
    for quantity in [:gm2]
        eqt.profiles_1d.gm2[1] = Interpolations.CubicSplineInterpolation(eqt.profiles_1d.psi[2:end], getproperty(eqt.profiles_1d,quantity)[2:end], extrapolation_bc=Interpolations.Line()).(eqt.profiles_1d.psi[1])
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
    flux_surface(eqt::equilibrium__time_slice, psi_level::Real, closed::Bool)

returns r,z coordiates of open or closed flux surface at given psi_level
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real, closed::Bool)
    dim1=eqt.profiles_2d[1].grid.dim1
    dim2=eqt.profiles_2d[1].grid.dim2
    PSI=eqt.profiles_2d[1].psi
    psi=eqt.profiles_1d.psi
    r0 = eqt.global_quantities.magnetic_axis.r
    z0 =eqt.global_quantities.magnetic_axis.z
    flux_surface(dim1, dim2, PSI, psi, r0, z0, psi_level, closed)
end

function flux_surface(dim1::Union{AbstractVector,AbstractRange},
                      dim2::Union{AbstractVector,AbstractRange},
                      PSI::AbstractArray,
                      psi::Union{AbstractVector,AbstractRange},
                      r0::Real,
                      z0::Real,
                      psi_level::Real,
                      closed::Bool)
    # handle on axis value as the first flux surface
    if psi_level == psi[1]
        psi_level = psi[2]
    # handle boundary by finding accurate lcfs psi
    elseif psi_level == psi[end]
        psi__boundary_level = find_psi_boundary(dim1, dim2, PSI, psi, r0, z0; raise_error_on_not_open=false)
        if psi__boundary_level !== nothing
            if (abs(psi__boundary_level-psi_level)<abs(psi[end]-psi[end-1]))
                psi_level=psi__boundary_level
            end
        end
    end

    cl = Contour.contour(dim1, dim2, PSI, psi_level)

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
        if (! closed) || (PolygonOps.inpolygon((r0, z0), StaticArrays.SVector.(pr, pz)) == 1)
            return pr, pz, psi_level
        end
    end
    return [], [], psi_level
end

"""
    find_psi_boundary(eqt; precision=1e-6, raise_error_on_not_open=true)

Find psi value of the last closed flux surface
"""
function find_psi_boundary(eqt; precision=1e-6, raise_error_on_not_open=true)
    dim1=eqt.profiles_2d[1].grid.dim1
    dim2=eqt.profiles_2d[1].grid.dim2
    PSI=eqt.profiles_2d[1].psi
    psi=eqt.profiles_1d.psi
    r0 = eqt.global_quantities.magnetic_axis.r
    z0 =eqt.global_quantities.magnetic_axis.z
    find_psi_boundary(dim1, dim2, PSI, psi, r0, z0; precision, raise_error_on_not_open)
end

function find_psi_boundary(dim1, dim2, PSI, psi, r0, z0; precision, raise_error_on_not_open)
    psirange_init = [psi[1]*0.9+psi[end]*0.1, psi[end] + 0.5 * (psi[end] - psi[1])]

    pr, pz = flux_surface(dim1, dim2, PSI, psi, r0, z0, psirange_init[1], true)
    if length(pr)==0
        error("Flux surface at ψ=$(psirange_init[1]) is not closed")
    end

    pr, pz = flux_surface(dim1, dim2, PSI, psi, r0, z0, psirange_init[end], true)
    if length(pr)>0
        if raise_error_on_not_open
            error("Flux surface at ψ=$(psirange_init[end]) is not open")
        else
            return nothing
        end
    end

    psirange = deepcopy(psirange_init)
    for k in 1:100
        psimid = (psirange[1] + psirange[end]) / 2.0
        pr, pz = flux_surface(dim1, dim2, PSI, psi, r0, z0, psimid, true)
        # closed flux surface
        if length(pr) > 0
            psirange[1] = psimid
            if (abs(psirange[end] - psirange[1])/abs(psirange[end] + psirange[1])/2.0) < precision
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
    radial_build_radii(rb::IMAS.radial_build)

Return list of radii in the radial build
"""
function radial_build_radii(rb::IMAS.radial_build)
    layers_radii = Real[]
    layer_start = 0
    for l in rb.layer
        push!(layers_radii, layer_start)
        layer_start = layer_start + l.thickness
    end
    push!(layers_radii, layer_start)
    return layers_radii
end

"""
    get_radial_build(rb::IMAS.radial_build;
                     type::Union{Nothing,Int}=nothing,
                     name::Union{Nothing,String}=nothing,
                     identifier::Union{Nothing,UInt,Int}=nothing,
                     hfs::Union{Nothing,Int}=nothing,
                     return_only_one=true )

Select layer(s) in radial build based on a series of selection criteria
"""
function get_radial_build(rb::IMAS.radial_build;
                          type::Union{Nothing,Int}=nothing,
                          name::Union{Nothing,String}=nothing,
                          identifier::Union{Nothing,UInt,Int}=nothing,
                          hfs::Union{Nothing,Int,Array}=nothing,
                          return_only_one=true,
                          return_index=false,
                          raise_error_on_missing=true
                          )

    if isa(hfs,Int)
        hfs=[hfs]
    end
    
    valid_layers = []
    for (k,l) in enumerate(rb.layer)
        if (name===nothing || l.name == name!) && (type===nothing || l.type == type) && (identifier===nothing || l.identifier == identifier) && (hfs===nothing || l.hfs in hfs)
            if return_index
                push!(valid_layers, k)
            else
                push!(valid_layers, l)
            end
        end
    end
    if length(valid_layers)==0
        if raise_error_on_missing
            error("Did not find radial_build.layer layer name:$name type:$type identifier:$identifier hfs:$hfs")
        else
            return nothing
        end
    end
    if return_only_one
        if length(valid_layers)==1
            return valid_layers[1]
        else
            error("Found multiple layers that satisfy name:$name type:$type identifier:$identifier hfs:$hfs")
        end
    else
        return valid_layers
    end
end

"""
    structures_mask(rb::IMAS.radial_build, resolution::Int=257, border::Real=10/resolution)

return rmask, zmask, mask of structures that are not vacuum
"""
function structures_mask(rb::IMAS.radial_build; resolution::Int=257, border_fraction::Real=0.1, one_is_for_vacuum::Bool=false)
    border = maximum(rb.layer[end].outline.r)*border_fraction
    xlim = [0.0,maximum(rb.layer[end].outline.r)+border]
    ylim = [minimum(rb.layer[end].outline.z)-border,maximum(rb.layer[end].outline.z)+border]
    rmask = range(xlim[1], xlim[2], length=resolution)
    zmask = range(ylim[1], ylim[2], length=resolution * Int(round((ylim[2] - ylim[1]) / (xlim[2] - xlim[1]))))
    mask = ones(length(rmask), length(zmask))

    valid = true
    for layer in vcat(rb.layer[end],rb.layer)
        if layer.type == -1
            valid = false
        end
        if valid && ! is_missing(layer.outline,:r)
            outline = StaticArrays.SVector.(layer.outline.r,layer.outline.z)
            if ! is_missing(layer,:material) && layer.material == "vacuum"
                for (kr, rr) in enumerate(rmask)
                    for (kz, zz) in enumerate(zmask)
                        if PolygonOps.inpolygon((rr, zz), outline) == 1
                            mask[kr,kz] = 0.0
                        end
                    end
                end
            else
                for (kr, rr) in enumerate(rmask)
                    for (kz, zz) in enumerate(zmask)
                        if PolygonOps.inpolygon((rr, zz), outline) == 1
                            mask[kr,kz] = 1.0
                        end
                    end
                end
            end
        end
    end
    rlim_oh = IMAS.get_radial_build(rb,type=1).start_radius
    for (kr, rr) in enumerate(rmask)
        for (kz, zz) in enumerate(zmask)
            if rr<rlim_oh
                mask[kr,kz] = 1.0
            end
        end
    end
    if one_is_for_vacuum
        return rmask, zmask, 1.0 .- mask
    else
        return rmask, zmask, mask
    end
end
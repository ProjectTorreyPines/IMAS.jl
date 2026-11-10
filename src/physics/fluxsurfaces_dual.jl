# ================================================================== #
# Flux surfaces of an equilibrium whose data carry derivatives (duals) #
# ================================================================== #
# The contouring, the axis and boundary searches and the refinement of the extrema are iterative
# geometric algorithms: they run on value parts, exactly as for a Float64 equilibrium. What they find is
# defined by conditions on ψ(R,Z):
#
#   - a point of a flux surface:   ψ(x) = ψ_level
#   - magnetic axis and X-points:  ∇ψ(x) = 0
#   - outermost/innermost and top/bottom points of a surface:  ψ(x) = ψ_level  and  ∂ψ/∂Z = 0  (or ∂ψ/∂R = 0)
#   - last closed surface:         ψ_level = ψ at the X-point or at the wall contact point that limits it
#
# so when ψ changes by δψ(R,Z) and the level by δψ_level each of these points moves by an amount that
# follows from the implicit function theorem. The points are given that displacement as their derivative
# part, and every flux-surface average and geometric quantity is then computed from them by the same code
# as for a Float64 equilibrium, which propagates the derivatives.

# dual number with zero value and the derivative part of x
_derivative_part(x::D) where {D<:Dual} = x - ForwardDiff.value(x)

# ψ gradient at a point that itself carries derivatives: first-order expansion around its value
function _gradient_at(itp, R::D, Z::D, H::Matrix{Float64}) where {D<:Dual}
    Rv = ForwardDiff.value(R)
    Zv = ForwardDiff.value(Z)
    g = _gradient(itp, Rv, Zv)
    _hessian_values!(H, itp, Rv, Zv)
    dR = R - Rv
    dZ = Z - Zv
    return (g[1] + H[1, 1] * dR + H[1, 2] * dZ, g[2] + H[2, 1] * dR + H[2, 2] * dZ)
end

function _hessian_values!(H::Matrix{Float64}, itp, R::Float64, Z::Float64)
    HD = FI.hessian(itp, (R, Z))
    for j in 1:2, i in 1:2
        H[i, j] = force_float64(HD[i, j])
    end
    return H
end

# Point (R, Z) of the traced contour ψ = level.
# The contouring places its points on the edges of the grid, by linear interpolation of ψ between the two
# nodes of the edge. The same formula with dual ψ and level gives the point with its derivative. Returns
# `nothing` when (R, Z) is not such a point.
function _contour_point_on_edge(r::AbstractVector{Float64}, z::AbstractVector{Float64}, PSI::Matrix{D}, R::Float64, Z::Float64, level::D) where {D<:Dual}
    ir = clamp(searchsortedlast(r, R), 1, length(r) - 1)
    jz = clamp(searchsortedlast(z, Z), 1, length(z) - 1)
    dr = r[ir+1] - r[ir]
    dz = z[jz+1] - z[jz]
    tol = 1E-8
    for i in (ir, ir + 1) # on a grid line R = r[i]: the point moves in Z
        if abs(R - r[i]) <= tol * dr
            t = (level - PSI[i, jz]) / (PSI[i, jz+1] - PSI[i, jz])
            Zd = z[jz] + t * dz
            if abs(ForwardDiff.value(Zd) - Z) <= 1E-6 * dz
                return (D(r[i]), Zd)
            end
        end
    end
    for j in (jz, jz + 1) # on a grid line Z = z[j]: the point moves in R
        if abs(Z - z[j]) <= tol * dz
            t = (level - PSI[ir, j]) / (PSI[ir+1, j] - PSI[ir, j])
            Rd = r[ir] + t * dr
            if abs(ForwardDiff.value(Rd) - R) <= 1E-6 * dr
                return (Rd, D(z[j]))
            end
        end
    end
    return nothing
end

# point (R, Z) of the surface ψ = level, displaced along ∇ψ so that it stays on the surface
function _implicit_contour_point(itp, R::Float64, Z::Float64, δlevel::D) where {D<:Dual}
    ψ, g = _value_gradient(itp, R, Z)
    gR = ForwardDiff.value(g[1])
    gZ = ForwardDiff.value(g[2])
    s = (δlevel - _derivative_part(ψ)) / (gR^2 + gZ^2)
    return (R + s * gR, Z + s * gZ)
end

# stationary point of ψ (magnetic axis, X-point)
function _implicit_stationary_point(itp, R::Float64, Z::Float64, H::Matrix{Float64}, ::Type{D}) where {D<:Dual}
    g = _gradient(itp, R, Z)
    _hessian_values!(H, itp, R, Z)
    det = H[1, 1] * H[2, 2] - H[1, 2] * H[2, 1]
    δg1 = _derivative_part(D(g[1]))
    δg2 = _derivative_part(D(g[2]))
    return (R - (H[2, 2] * δg1 - H[1, 2] * δg2) / det, Z - (-H[2, 1] * δg1 + H[1, 1] * δg2) / det)
end

# extremum in R (dir = :R, where ∂ψ/∂Z = 0) or in Z (dir = :Z, where ∂ψ/∂R = 0) of the surface ψ = level
function _implicit_extremum(itp, R::Float64, Z::Float64, dir::Symbol, δlevel::D, H::Matrix{Float64}) where {D<:Dual}
    ψ, g = _value_gradient(itp, R, Z)
    _hessian_values!(H, itp, R, Z)
    c = dir === :R ? 2 : 1 # component of the gradient that vanishes
    a11, a12 = ForwardDiff.value(g[1]), ForwardDiff.value(g[2])
    a21, a22 = H[c, 1], H[c, 2]
    b1 = δlevel - _derivative_part(ψ)
    b2 = -_derivative_part(D(g[c]))
    det = a11 * a22 - a12 * a21
    return (R + (a22 * b1 - a12 * b2) / det, Z + (-a21 * b1 + a11 * b2) / det)
end

"""
    flux_surfaces(eqt::equilibrium__time_slice{D}, wall_r, wall_z) where {D<:ForwardDiff.Dual}

Flux surface averaged and geometric quantities of an equilibrium time slice whose data carry derivatives.

The value part of every result is that of the Float64 calculation. The derivative parts follow from the
conditions that define the axis, the boundary, the flux surfaces and their extrema (implicit function
theorem), and are then propagated through the flux-surface averages. The wall is machine geometry and is
taken as given. The positions of strike points and the outlines of the separatrices are diagnostics of
the geometry and are returned without derivative parts.
"""
function flux_surfaces(eqt::equilibrium__time_slice{D}, wall_r::AbstractVector{T2}, wall_z::AbstractVector{T2}) where {D<:ForwardDiff.Dual,T2<:Real}
    wr = Float64[force_float64(x) for x in wall_r]
    wz = Float64[force_float64(x) for x in wall_z]

    # ---- Float64 calculation: values of everything, and of the points that define the geometry ----
    # (the copy of the time slice lives in a data dictionary of its own; the inputs that are expressions
    #  of other parts of the original data dictionary are evaluated there and stored as data)
    dd_v = dd{Float64}()
    resize!(dd_v.equilibrium.time_slice, 1)
    eqt_v = dd_v.equilibrium.time_slice[1]
    fill!(eqt_v, eqt)
    dd_v.equilibrium.time = [eqt.time]
    dd_v.global_time = eqt.time
    eqt_v.global_quantities.vacuum_toroidal_field.r0 = force_float64(eqt.global_quantities.vacuum_toroidal_field.r0)
    eqt_v.global_quantities.vacuum_toroidal_field.b0 = force_float64(eqt.global_quantities.vacuum_toroidal_field.b0)
    for field in (:psi, :f, :pressure, :dpressure_dpsi, :f_df_dpsi)
        setproperty!(eqt_v.profiles_1d, field, Float64[force_float64(x) for x in getproperty(eqt.profiles_1d, field)])
    end
    xpoints_in = [(xp.r, xp.z) for xp in eqt_v.boundary.x_point]
    psi_in = eqt.profiles_1d.psi
    flux_surfaces(eqt_v, wr, wz)

    eqt1d = eqt.profiles_1d
    eqt1d_v = eqt_v.profiles_1d
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    eqt2d_v = findfirst(:rectangular, eqt_v.profiles_2d)
    r, z, itp_v = ψ_interpolant(eqt2d_v)
    itp = ψ_interpolant(r, z, eqt2d.psi).PSI_interpolant # dual ψ on the plain grid
    H = Matrix{Float64}(undef, 2, 2)

    psi_sign = sign(eqt1d_v.psi[end] - eqt1d_v.psi[1])
    R0 = eqt.global_quantities.vacuum_toroidal_field.r0
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

    for field in (:psi_boundary, :psi_axis, :q_95, :q_axis, :q_min)
        empty!(eqt.global_quantities, field)
    end

    # ---- magnetic axis ----
    RAv = eqt_v.global_quantities.magnetic_axis.r
    ZAv = eqt_v.global_quantities.magnetic_axis.z
    RA, ZA = _implicit_stationary_point(itp, RAv, ZAv, H, D)
    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = RA, ZA
    psi_axis = D(itp(RAv, ZAv)) # ∇ψ = 0 at the axis: its displacement does not change ψ there to first order

    # ---- last closed flux surface ----
    psi_b_v = eqt_v.boundary.psi
    rescaled = eqt_v.boundary_separatrix.psi !== eqt_v.boundary.psi
    if !isempty(wr) || eqt.global_quantities.free_boundary == 1
        # the surface is limited by an X-point or by the wall: its ψ follows ψ at that point
        best = (Inf, RAv, ZAv)
        for xp in eqt_v.boundary.x_point
            best = min(best, (abs(itp_v(xp.r, xp.z) - psi_b_v), xp.r, xp.z))
        end
        if best[1] > 1E-3 * abs(psi_b_v - force_float64(psi_axis)) && !isempty(wr)
            best = (Inf, RAv, ZAv)
            for k in eachindex(wr)
                k2 = k == length(wr) ? 1 : k + 1
                for t in range(0.0, 1.0, 21)
                    Rw = wr[k] + t * (wr[k2] - wr[k])
                    Zw = wz[k] + t * (wz[k2] - wz[k])
                    if r[1] <= Rw <= r[end] && z[1] <= Zw <= z[end]
                        best = min(best, (abs(itp_v(Rw, Zw) - psi_b_v), Rw, Zw))
                    end
                end
            end
        end
        # The contours are traced cell by cell on the grid. In the cell that holds an X-point the two
        # branches of the contour are connected according to the mean of ψ at the four nodes of the cell,
        # so that mean, not the saddle value of the smooth interpolant, is the level at which the
        # contour opens up.
        psi_limit = D(itp(best[2], best[3]))
        ic = clamp(searchsortedlast(r, best[2]), 1, length(r) - 1)
        jc = clamp(searchsortedlast(z, best[3]), 1, length(z) - 1)
        cell = (Inf, 0, 0)
        for i in max(1, ic - 2):min(length(r) - 1, ic + 2), j in max(1, jc - 2):min(length(z) - 1, jc + 2)
            cell_mean = (eqt2d_v.psi[i, j] + eqt2d_v.psi[i+1, j] + eqt2d_v.psi[i+1, j+1] + eqt2d_v.psi[i, j+1]) / 4
            cell = min(cell, (abs(cell_mean - psi_b_v), i, j))
        end
        if cell[1] <= 1E-4 * abs(psi_b_v - force_float64(psi_axis))
            i, j = cell[2], cell[3]
            psi_limit = (eqt2d.psi[i, j] + eqt2d.psi[i+1, j] + eqt2d.psi[i+1, j+1] + eqt2d.psi[i, j+1]) / 4
        end
        # the search stops a small, fixed fraction of the way inside the limiting value
        fraction = (psi_b_v - force_float64(psi_axis)) / (force_float64(psi_limit) - force_float64(psi_axis))
        psi_boundary = psi_b_v + _derivative_part(psi_axis + fraction * (psi_limit - psi_axis))
    else
        psi_boundary = psi_in[end]
    end
    if rescaled
        eqt1d.psi = (psi_in .- psi_in[1]) ./ (psi_in[end] - psi_in[1]) .* (psi_boundary - psi_axis) .+ psi_axis
    end
    eqt.boundary.psi = psi_boundary
    eqt.boundary_separatrix.psi = eqt_v.boundary_separatrix.psi + _derivative_part(psi_boundary)

    for item in (:b_field_average, :b_field_max, :b_field_min, :elongation, :triangularity_lower, :triangularity_upper,
        :squareness_lower_inner, :squareness_lower_outer, :squareness_upper_inner, :squareness_upper_outer, :r_inboard, :r_outboard,
        :q, :surface, :dvolume_dpsi, :j_tor, :gm1, :gm2, :gm4, :gm5, :gm8, :gm9, :gm10, :fsa_bp, :trapped_fraction)
        setproperty!(eqt1d, item, zeros(D, size(eqt1d.psi)))
    end

    # ---- flux surfaces: traced on values, then given their displacement ----
    surfaces_v = trace_surfaces(eqt1d_v.psi, eqt1d_v.f, r, z, eqt2d_v.psi, itp_v, RAv, ZAv, wr, wz; refine_extrema=true, xpoints=xpoints_in)
    N = length(surfaces_v)
    surfaces = Vector{FluxSurface{D}}(undef, N)
    f = eqt1d.f
    for k in N:-1:1
        sv = surfaces_v[k]
        np = length(sv.r)
        pr = Vector{D}(undef, np)
        pz = Vector{D}(undef, np)
        if k == 1
            # artificial surface around the axis, scaled from the second one
            pr .= (surfaces[2].r .- RA) ./ 100.0 .+ RA
            pz .= (surfaces[2].z .- ZA) ./ 100.0 .+ ZA
        else
            level = eqt1d.psi[k]
            δlevel = _derivative_part(level)
            on_edges = 0
            for i in 1:np
                point = _contour_point_on_edge(r, z, eqt2d.psi, sv.r[i], sv.z[i], level)
                if point === nothing
                    point = _implicit_contour_point(itp, sv.r[i], sv.z[i], δlevel)
                else
                    on_edges += 1
                end
                pr[i], pz[i] = point
            end
            if k == 2 && on_edges == 0
                # the innermost surface could not be traced and was scaled from the next one
                pr .= (surfaces[k+1].r .- RA) ./ 2 .+ RA
                pz .= (surfaces[k+1].z .- ZA) ./ 2 .+ ZA
            end
        end

        Br = Vector{D}(undef, np)
        Bz = Vector{D}(undef, np)
        Bp = Vector{D}(undef, np)
        Btot = Vector{D}(undef, np)
        fluxexpansion = Vector{D}(undef, np)
        for i in 1:np
            g = _gradient_at(itp, pr[i], pz[i], H)
            inv_twopi_r = 1.0 / (2π * pr[i])
            Br[i] = g[2] * inv_twopi_r
            Bz[i] = -g[1] * inv_twopi_r
            Bp2 = Br[i]^2 + Bz[i]^2
            Bp_abs = sqrt(Bp2)
            Bp[i] = Bp_abs * sign((pz[i] - ZA) * Br[i] - (pr[i] - RA) * Bz[i])
            Btot[i] = sqrt(Bp2 + (f[k] / pr[i])^2)
            fluxexpansion[i] = 1.0 / Bp_abs
        end
        ll = arc_length(pr, pz)
        int_fluxexpansion_dl = trapz(ll, fluxexpansion)

        if k == 1
            frac = 0.01
            s2 = surfaces[2]
            ext = ((s2.r_at_max_z - RA) * frac + RA, (s2.max_z - ZA) * frac + ZA, (s2.r_at_min_z - RA) * frac + RA, (s2.min_z - ZA) * frac + ZA,
                (s2.z_at_max_r - ZA) * frac + ZA, (s2.max_r - RA) * frac + RA, (s2.z_at_min_r - ZA) * frac + ZA, (s2.min_r - RA) * frac + RA)
        else
            δlevel = _derivative_part(eqt1d.psi[k])
            r_at_max_z, max_z = _implicit_extremum(itp, sv.r_at_max_z, sv.max_z, :Z, δlevel, H)
            r_at_min_z, min_z = _implicit_extremum(itp, sv.r_at_min_z, sv.min_z, :Z, δlevel, H)
            max_r, z_at_max_r = _implicit_extremum(itp, sv.max_r, sv.z_at_max_r, :R, δlevel, H)
            min_r, z_at_min_r = _implicit_extremum(itp, sv.min_r, sv.z_at_min_r, :R, δlevel, H)
            ext = (r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
        end

        surfaces[k] = FluxSurface(eqt1d.psi[k], pr, pz, ext..., Br, Bz, Bp, Btot, ll, fluxexpansion, int_fluxexpansion_dl)
    end

    _flux_surfaces_from_traced!(eqt, eqt2d, surfaces, psi_sign, R0, B0)

    # ---- X-points, strike points, separatrices: values from the Float64 calculation ----
    fill!(eqt.boundary.x_point, eqt_v.boundary.x_point)
    for xp in eqt.boundary.x_point
        xp.r, xp.z = _implicit_stationary_point(itp, force_float64(xp.r), force_float64(xp.z), H, D)
    end
    fill!(eqt.boundary.strike_point, eqt_v.boundary.strike_point)
    psi_separatrix = eqt.boundary_separatrix.psi
    fill!(eqt.boundary_separatrix, eqt_v.boundary_separatrix)
    eqt.boundary_separatrix.psi = psi_separatrix
    fill!(eqt.boundary_secondary_separatrix, eqt_v.boundary_secondary_separatrix)

    # the value part of every result is that of the Float64 calculation
    replace_values!(eqt, eqt_v)

    return eqt
end

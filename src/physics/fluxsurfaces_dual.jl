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

# magnetic axis of a ψ interpolant whose data carry derivatives: search on values, then the displacement
# of the stationary point
function _dual_magnetic_axis(itp, psi_sign, rguess::Float64, zguess::Float64)
    s = force_float64(psi_sign)
    res = Optim.optimize(x -> force_float64(itp(x[1], x[2])) * s, [rguess, zguess], Optim.Newton())
    D = typeof(itp(rguess, zguess))
    H = Matrix{Float64}(undef, 2, 2)
    RA, ZA = _implicit_stationary_point(itp, res.minimizer[1], res.minimizer[2], H, D)
    return (RA=RA, ZA=ZA)
end

# levels of the last closed and first open flux surfaces found by the bisection of `find_psi_boundary`
function _dual_psi_boundaries(dimR, dimZ, PSI::Matrix{D}, itp, psi_axis, psimid, psiend, fw_r, fw_z) where {D<:Dual}
    r = Float64[force_float64(x) for x in dimR]
    z = Float64[force_float64(x) for x in dimZ]
    wr = Float64[force_float64(x) for x in fw_r]
    wz = Float64[force_float64(x) for x in fw_z]
    last_closed = _implicit_psi_boundary(r, z, PSI, itp, D(psi_axis), force_float64(psimid), force_float64(psiend), wr, wz)
    first_open = force_float64(psiend) + _derivative_part(last_closed)
    return (last_closed=last_closed, first_open=first_open)
end

# Level of the last closed flux surface as a dual number, given the bracket found by the bisection:
# `psi_b_v` is the last level whose contour is closed and `psi_open_v` the first one whose contour is open.
#
# The contours are traced cell by cell on the grid, with linear interpolation of ψ along the edges of the
# cells. Such a contour changes topology only when its level passes the value of ψ at a grid node, or, in
# a cell where ψ has a saddle, the mean of ψ at the four nodes of the cell (which decides how the two
# branches are connected). The level at which the contour opens up is therefore that node value or cell
# mean, and it is the one that lies inside the bracket. When there is none, the surface is limited by
# the wall, and the level follows ψ at the point of the wall that it touches. The bisection stops a
# small, fixed fraction of the way inside that limit.
# Saddle point of ψ (X-point) whose flux is the level `psi_b_v` to within `rtol` of the flux between axis and
# boundary: Newton iterations on ∇ψ = 0 from the grid nodes near that level where the gradient is smallest.
function _saddle_at_level(r::AbstractVector{Float64}, z::AbstractVector{Float64}, PSI::Matrix{D}, itp, psi_a_v::Float64, psi_b_v::Float64; rtol::Float64=1E-2) where {D<:Dual}
    span = abs(psi_b_v - psi_a_v)
    dr = r[2] - r[1]
    dz = z[2] - z[1]
    candidates = Tuple{Float64,Int,Int}[]
    for j in 2:length(z)-1, i in 2:length(r)-1
        if abs(ForwardDiff.value(PSI[i, j]) - psi_b_v) <= 5 * rtol * span
            gR = (ForwardDiff.value(PSI[i+1, j]) - ForwardDiff.value(PSI[i-1, j])) / (2dr)
            gZ = (ForwardDiff.value(PSI[i, j+1]) - ForwardDiff.value(PSI[i, j-1])) / (2dz)
            push!(candidates, (hypot(gR, gZ), i, j))
        end
    end
    sort!(candidates)
    H = zeros(2, 2)
    best = nothing
    for (_, i, j) in candidates[1:min(8, end)]
        R, Z = r[i], z[j]
        ok = false
        for _ in 1:30
            g = _gradient(itp, R, Z)
            g1 = force_float64(g[1])
            g2 = force_float64(g[2])
            _hessian_values!(H, itp, R, Z)
            det = H[1, 1] * H[2, 2] - H[1, 2] * H[2, 1]
            det < 0 || break
            dR = -(H[2, 2] * g1 - H[1, 2] * g2) / det
            dZ = -(-H[2, 1] * g1 + H[1, 1] * g2) / det
            R += dR
            Z += dZ
            (abs(R - r[i]) <= 2dr && abs(Z - z[j]) <= 2dz) || break
            if hypot(dR / dr, dZ / dz) < 1E-10
                ok = true
                break
            end
        end
        ok || continue
        miss = abs(force_float64(itp(R, Z)) - psi_b_v)
        if miss <= rtol * span && (best === nothing || miss < best[1])
            best = (miss, R, Z)
        end
    end
    return best === nothing ? nothing : (best[2], best[3])
end

function _implicit_psi_boundary(
    r::AbstractVector{Float64}, z::AbstractVector{Float64}, PSI::Matrix{D}, itp, psi_axis::D, psi_b_v::Float64, psi_open_v::Float64,
    wall_r::AbstractVector{Float64}, wall_z::AbstractVector{Float64}; xpoint::Bool=true) where {D<:Dual}

    psi_a_v = ForwardDiff.value(psi_axis)
    lo, hi = minmax(psi_b_v, psi_open_v)
    pad = 1E-3 * (hi - lo)
    lo -= pad
    hi += pad
    psi_limit = nothing
    if xpoint
        # diverted plasma: the last closed surface follows the flux at the X-point. The level found by the
        # search differs from it by an offset set by the contouring grid, which is piecewise in the inputs;
        # the derivative is that of the X-point flux, which is the trend over any finite change of the inputs.
        saddle = _saddle_at_level(r, z, PSI, itp, psi_a_v, psi_b_v)
        if saddle !== nothing
            psi_limit = D(itp(saddle[1], saddle[2]))
        end
    end
    for j in eachindex(z), i in eachindex(r)
        psi_limit === nothing || break
        if lo <= ForwardDiff.value(PSI[i, j]) <= hi
            psi_limit = PSI[i, j]
            break
        end
    end
    if psi_limit === nothing
        for j in 1:length(z)-1, i in 1:length(r)-1
            m = (PSI[i, j] + PSI[i+1, j] + PSI[i+1, j+1] + PSI[i, j+1]) / 4
            if lo <= ForwardDiff.value(m) <= hi
                psi_limit = m
                break
            end
        end
    end
    if psi_limit === nothing && !isempty(wall_r)
        touch = (Inf, 0.0, 0.0)
        for k in eachindex(wall_r)
            k2 = k == length(wall_r) ? 1 : k + 1
            for t in range(0.0, 1.0, 21)
                Rw = wall_r[k] + t * (wall_r[k2] - wall_r[k])
                Zw = wall_z[k] + t * (wall_z[k2] - wall_z[k])
                if r[1] <= Rw <= r[end] && z[1] <= Zw <= z[end]
                    touch = min(touch, (abs(force_float64(itp(Rw, Zw)) - psi_b_v), Rw, Zw))
                end
            end
        end
        psi_limit = D(itp(touch[2], touch[3]))
    end
    if psi_limit === nothing
        return D(psi_b_v)
    end
    fraction = (psi_b_v - psi_a_v) / (ForwardDiff.value(psi_limit) - psi_a_v)
    return psi_b_v + _derivative_part(psi_axis + fraction * (psi_limit - psi_axis))
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
        psi_boundary = _implicit_psi_boundary(collect(r), collect(z), eqt2d.psi, itp, psi_axis, psi_b_v, eqt_v.boundary_separatrix.psi, wr, wz)
        # when the level follows an X-point the points of the last surface are displaced on the smooth ψ(R,Z)
        # rather than along the grid edges, where the linear interpolation is poor next to the X-point
        smooth_boundary = _saddle_at_level(collect(r), collect(z), eqt2d.psi, itp, ForwardDiff.value(psi_axis), psi_b_v) !== nothing
    else
        psi_boundary = psi_in[end]
        smooth_boundary = false
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
                point = (k == N && smooth_boundary) ? nothing : _contour_point_on_edge(r, z, eqt2d.psi, sv.r[i], sv.z[i], level)
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

# The contouring routine takes its tolerances in the number type of the grid, with defaults built from
# `eps` of that type, which is a plain float also for dual numbers. This method supplies them in the
# right type when grid and data are dual numbers.
function IMASutils.contour_from_midplane!(
    x_cache::Vector{D}, y_cache::Vector{D}, values::Matrix{D}, x_coords::AbstractVector{D}, y_coords::AbstractVector{D},
    level::D, xaxis::D, yaxis::D, vaxis::D;
    atol::Real=eps(Float64), rtol::Real=sqrt(eps(Float64)), dynamic_resize::Bool=true) where {D<:Dual}
    return invoke(IMASutils.contour_from_midplane!,
        Tuple{Vector{T1},Vector{T1},Matrix{T2},AbstractVector{T1},AbstractVector{T1},T2,T1,T1,T2} where {T1<:Real,T2<:Real},
        x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis; atol=D(atol), rtol=D(rtol), dynamic_resize)
end

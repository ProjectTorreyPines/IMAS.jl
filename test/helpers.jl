# Shared helpers for the IMAS test suite
# (asymmetric synthetic boundaries and transformed copies of the sample DIII-D equilibrium)
using IMAS

# Luce superellipse boundary with prescribed per-quadrant squareness
# T.C. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009
function luce_boundary(; R0=1.7, a=0.6, Zrm=0.0, κu=1.8, κl=1.9, δu=0.4, δl=0.5,
    ζ=(ou=0.10, iu=-0.20, il=0.25, ol=-0.05), n=401)
    nexp(ζ) = -log(2) / log(1 / sqrt(2) + ζ * (1 - 1 / sqrt(2)))
    pr = Float64[]
    pz = Float64[]
    for θ in range(0, 2π; length=n)[1:end-1]
        if θ < π / 2          # upper outer
            t = θ; A = a * (1 + δu); B = κu * a; nn = nexp(ζ.ou); Rc = R0 - a * δu; sx = +1; sz = +1
        elseif θ < π          # upper inner
            t = π - θ; A = a * (1 - δu); B = κu * a; nn = nexp(ζ.iu); Rc = R0 - a * δu; sx = -1; sz = +1
        elseif θ < 3π / 2     # lower inner
            t = θ - π; A = a * (1 - δl); B = κl * a; nn = nexp(ζ.il); Rc = R0 - a * δl; sx = -1; sz = -1
        else                  # lower outer
            t = 2π - θ; A = a * (1 + δl); B = κl * a; nn = nexp(ζ.ol); Rc = R0 - a * δl; sx = +1; sz = -1
        end
        x = A * abs(cos(t))^(2 / nn)
        y = B * abs(sin(t))^(2 / nn)
        push!(pr, Rc + sx * x)
        push!(pz, Zrm + sz * y)
    end
    return pr, pz
end

"""
    fresh_sample()

The sample DIII-D equilibrium as read from file (no flux_surfaces call yet) and its first wall.
flux_surfaces is not idempotent at the 1e-4 level (a second call starts from the refined boundary psi),
so invariance tests must run it once on a fresh copy and once on a fresh transformed copy.
"""
function fresh_sample()
    samp = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")
    dd = IMAS.json2imas(joinpath(samp, "D3D_eq_ods.json"); show_warnings=false)
    return dd, IMAS.first_wall(dd.wall)
end

function sample_eqt_with_wall()
    samp = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")
    dd = IMAS.json2imas(joinpath(samp, "D3D_eq_ods.json"); show_warnings=false)
    fw = IMAS.first_wall(dd.wall)
    eqt = dd.equilibrium.time_slice[]
    IMAS.flux_surfaces(eqt, fw.r, fw.z)
    _, _, PSI = IMAS.ψ_interpolant(eqt.profiles_2d)
    return dd, eqt, fw, PSI
end

"""
    detached_copy(eqt)

Deep copy of an equilibrium time slice that keeps its parent `dd` (so that expressions such as
vacuum_toroidal_field.b0 can still be evaluated), or of the bare time slice if it has no parent.
"""
const _keep_alive = Any[] # IDS parents are WeakRefs: keep the copied dd alive for as long as the test process runs

function detached_copy(eqt)
    dd = IMAS.top_dd(eqt)
    if dd === nothing
        return deepcopy(eqt)
    end
    dd2 = deepcopy(dd)
    push!(_keep_alive, dd2)
    return dd2.equilibrium.time_slice[IMAS.index(eqt)]
end

function shifted_eqt_z(eqt, dz)
    e2 = detached_copy(eqt)
    e2d = findfirst(:rectangular, e2.profiles_2d)
    e2d.grid.dim2 = e2d.grid.dim2 .+ dz
    e2.global_quantities.magnetic_axis.z += dz
    e2.boundary.outline.z = e2.boundary.outline.z .+ dz
    # stored (not expression) geometric axes must move too
    if IMAS.hasdata(e2.profiles_1d.geometric_axis, :z)
        e2.profiles_1d.geometric_axis.z = e2.profiles_1d.geometric_axis.z .+ dz
    end
    if IMAS.hasdata(e2.boundary.geometric_axis, :z)
        e2.boundary.geometric_axis.z += dz
    end
    for sp in e2.boundary.strike_point
        sp.z += dz
    end
    for xp in e2.boundary.x_point
        xp.z += dz
    end
    return e2
end

"""
    reflected_eqt_z(eqt, fw)

Up/down reflection of an equilibrium time slice and of a first wall (r, z) about Z = 0.
psi(R, Z) -> psi(R, -Z); x-points, strike points, axis and boundary outline mirrored.
"""
function reflected_eqt_z(eqt, fw)
    e2 = detached_copy(eqt)
    e2d = findfirst(:rectangular, e2.profiles_2d)
    e2d.grid.dim2 = reverse(-e2d.grid.dim2)
    e2d.psi = e2d.psi[:, end:-1:1]
    for name in (:b_field_r, :b_field_z, :b_field_tor, :phi, :j_tor, :j_parallel)
        if !isempty(e2d, name)
            v = getproperty(e2d, name)[:, end:-1:1]
            setproperty!(e2d, name, name == :b_field_r ? -v : v)
        end
    end
    e2.global_quantities.magnetic_axis.z = -e2.global_quantities.magnetic_axis.z
    e2.boundary.outline.z = -e2.boundary.outline.z
    if IMAS.hasdata(e2.profiles_1d.geometric_axis, :z)
        e2.profiles_1d.geometric_axis.z = -e2.profiles_1d.geometric_axis.z
    end
    if IMAS.hasdata(e2.boundary.geometric_axis, :z)
        e2.boundary.geometric_axis.z = -e2.boundary.geometric_axis.z
    end
    for sp in e2.boundary.strike_point
        sp.z = -sp.z
    end
    for xp in e2.boundary.x_point
        xp.z = -xp.z
    end
    return e2, (r=copy(fw.r), z=-fw.z)
end

"""
    negated_psi_eqt(eqt)

Copy of an equilibrium time slice with the sign of the poloidal flux flipped everywhere it is stored
(equivalent to reversing the plasma current direction, as far as psi-based geometry is concerned).
"""
function negated_psi_eqt(eqt)
    e2 = detached_copy(eqt)
    e2d = findfirst(:rectangular, e2.profiles_2d)
    e2d.psi = -e2d.psi
    for name in (:b_field_r, :b_field_z)
        if IMAS.hasdata(e2d, name)
            setproperty!(e2d, name, -getproperty(e2d, name))
        end
    end
    e2.profiles_1d.psi = -e2.profiles_1d.psi
    for ids in (e2.boundary, e2.boundary_separatrix, e2.boundary_secondary_separatrix)
        if IMAS.hasdata(ids, :psi)
            ids.psi = -ids.psi
        end
    end
    for name in (:psi_axis, :psi_boundary)
        if IMAS.hasdata(e2.global_quantities, name)
            setproperty!(e2.global_quantities, name, -getproperty(e2.global_quantities, name))
        end
    end
    return e2
end

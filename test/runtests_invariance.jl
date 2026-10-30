using IMAS
using Test

include("helpers.jl")

# ---------------------------------------------------------------------------------------------------------------------
# Structural invariance tests for the shape / flux-surface pipeline.
#
# Every function here is fed an ASYMMETRIC input and then the same input transformed by
#   * a vertical translation (Z -> Z + dz),
#   * an up/down reflection (Z -> -Z), and
#   * a reversal of the point order.
# The outputs must transform the way their labels say: "upper" and "lower" quantities swap under reflection,
# nothing depends on the offset, and orientation does not matter. A swapped label, a mirror about the wrong axis
# or a Z = 0 reference changes the answer and fails these tests.
# ---------------------------------------------------------------------------------------------------------------------

@testset "invariance" begin

@testset "curvature invariances" begin
    θ = range(0, 2π; length=181)
    pr = 1.7 .+ 0.6 .* cos.(θ .+ 0.35 .* sin.(θ))
    pz = 0.3 .+ 1.1 .* sin.(θ) .+ 0.1 .* sin.(2θ)        # up/down asymmetric, centred off Z = 0
    pr[end], pz[end] = pr[1], pz[1]
    c = IMAS.curvature(pr, pz)
    @test IMAS.curvature(pr, pz .+ 1.3) ≈ c                     # translation
    @test IMAS.curvature(pr, -pz) ≈ -c                          # reflection flips the orientation, hence the sign
    @test IMAS.curvature(reverse(pr), reverse(pz)) ≈ -reverse(c) # reversal flips the sign, point by point
    # open path: same statements
    co = IMAS.curvature(pr[1:100], pz[1:100])
    @test IMAS.curvature(pr[1:100], pz[1:100] .+ 1.3) ≈ co
    @test IMAS.curvature(pr[1:100], -pz[1:100]) ≈ -co
end

@testset "fluxsurface_extrema / miller_R_a_κ_δ_ζ / luce_squareness invariances" begin
    ζ = (ou=0.10, iu=-0.20, il=0.25, ol=-0.05)
    kw = (R0=1.7, a=0.6, κu=1.8, κl=1.9, δu=0.4, δl=0.5, ζ=ζ)
    pr, pz = luce_boundary(; kw..., Zrm=0.0)
    g0 = IMAS.miller_R_a_κ_δ_ζ(pr, pz)

    # translation: nothing but the Z extrema moves
    dz = 0.73
    e0 = IMAS.fluxsurface_extrema(pr, pz)
    e1 = IMAS.fluxsurface_extrema(pr, pz .+ dz)
    @test e1.imaxr == e0.imaxr && e1.iminr == e0.iminr && e1.imaxz == e0.imaxz && e1.iminz == e0.iminz
    @test e1.max_z ≈ e0.max_z + dz && e1.min_z ≈ e0.min_z + dz
    @test e1.z_at_max_r ≈ e0.z_at_max_r + dz && e1.z_at_min_r ≈ e0.z_at_min_r + dz
    @test e1.r_at_max_z == e0.r_at_max_z && e1.r_at_min_z == e0.r_at_min_z && e1.max_r == e0.max_r && e1.min_r == e0.min_r
    g1 = IMAS.miller_R_a_κ_δ_ζ(pr, pz .+ dz)
    for k in keys(g0)
        @test getproperty(g1, k) ≈ getproperty(g0, k)
    end

    # reflection: upper <-> lower, everything else unchanged
    gr = IMAS.miller_R_a_κ_δ_ζ(pr, -pz)
    @test gr.R0 ≈ g0.R0 && gr.a ≈ g0.a && gr.κ ≈ g0.κ
    @test gr.δu ≈ g0.δl && gr.δl ≈ g0.δu
    @test isapprox(gr.ζou, g0.ζol; atol=2e-3) && isapprox(gr.ζol, g0.ζou; atol=2e-3)
    @test isapprox(gr.ζiu, g0.ζil; atol=2e-3) && isapprox(gr.ζil, g0.ζiu; atol=2e-3)
    # and the reflected shape is really the mirrored parameter set
    prr, pzr = luce_boundary(; kw..., κu=kw.κl, κl=kw.κu, δu=kw.δl, δl=kw.δu, ζ=(ou=ζ.ol, iu=ζ.il, il=ζ.iu, ol=ζ.ou))
    grr = IMAS.miller_R_a_κ_δ_ζ(prr, pzr)
    for k in keys(g0)
        @test isapprox(getproperty(grr, k), getproperty(gr, k); atol=2e-3)
    end

    # orientation reversal: unchanged
    gv = IMAS.miller_R_a_κ_δ_ζ(reverse(pr), reverse(pz))
    for k in keys(g0)
        @test isapprox(getproperty(gv, k), getproperty(g0, k); atol=1e-12)
    end

    # named access and positional access agree
    ext = IMAS.fluxsurface_extrema(pr, pz)
    @test (ext.imaxr, ext.iminr, ext.imaxz, ext.iminz) == (ext[1], ext[2], ext[3], ext[4])
    @test (ext.r_at_max_z, ext.max_z, ext.r_at_min_z, ext.min_z, ext.z_at_max_r, ext.max_r, ext.z_at_min_r, ext.min_r) == Tuple(ext)[5:12]
    lz = IMAS.luce_squareness(pr, pz, ext.r_at_max_z, ext.max_z, ext.r_at_min_z, ext.min_z, ext.z_at_max_r, ext.max_r, ext.z_at_min_r, ext.min_r)
    @test keys(lz) == (:zetaou, :zetaol, :zetail, :zetaiu)
    @test isapprox(lz.zetaou, ζ.ou; atol=2e-3) && isapprox(lz.zetaol, ζ.ol; atol=2e-3)
    @test isapprox(lz.zetaiu, ζ.iu; atol=2e-3) && isapprox(lz.zetail, ζ.il; atol=2e-3)
end

@testset "boundary_shape reflection and round trip" begin
    kw = (a=0.6, eps=0.35, kapu=1.9, kapl=1.7, delu=0.4, dell=0.5, zetaou=-0.1, zetaiu=-0.3, zetail=-0.2, zetaol=-0.15)
    zoffset = 0.21
    b = IMAS.boundary_shape(; kw..., zoffset)
    # round trip through the Miller fit must reproduce every parameter
    g = IMAS.miller_R_a_κ_δ_ζ(b.r, b.z)
    @test isapprox(g.a, kw.a; atol=1e-6)
    @test isapprox(g.R0, kw.a / kw.eps; atol=1e-6)
    @test isapprox(g.κ, (kw.kapu + kw.kapl) / 2; atol=1e-6)
    @test isapprox(g.δu, kw.delu; atol=1e-6)
    @test isapprox(g.δl, kw.dell; atol=1e-6)
    @test isapprox(g.ζou, kw.zetaou; atol=2e-3)
    @test isapprox(g.ζiu, kw.zetaiu; atol=2e-3)
    @test isapprox(g.ζil, kw.zetail; atol=2e-3)
    @test isapprox(g.ζol, kw.zetaol; atol=2e-3)
    @test isapprox(maximum(b.z), zoffset + kw.a * kw.kapu; atol=1e-9)
    @test isapprox(minimum(b.z), zoffset - kw.a * kw.kapl; atol=1e-9)

    # swapping the upper and lower parameters mirrors the outline about zoffset
    kwm = (a=kw.a, eps=kw.eps, kapu=kw.kapl, kapl=kw.kapu, delu=kw.dell, dell=kw.delu,
        zetaou=kw.zetaol, zetaiu=kw.zetail, zetail=kw.zetaiu, zetaol=kw.zetaou)
    bm = IMAS.boundary_shape(; kwm..., zoffset)
    gm = IMAS.miller_R_a_κ_δ_ζ(bm.r, 2 * zoffset .- bm.z)
    for k in keys(g)
        @test isapprox(getproperty(gm, k), getproperty(g, k); atol=2e-3)
    end
    # point by point: the mirrored outline is the same set of points (quadrant order differs)
    pts(r, z) = sort(collect(zip(round.(r; digits=9), round.(z; digits=9))))
    @test pts(bm.r, 2 * zoffset .- bm.z) == pts(b.r, b.z)
end

@testset "flux_surfaces translation invariance" begin
    dd, fw = fresh_sample()
    dz = 1.7
    e2 = shifted_eqt_z(dd.equilibrium.time_slice[], dz)
    eqt = dd.equilibrium.time_slice[]
    IMAS.flux_surfaces(eqt, fw.r, fw.z)
    IMAS.flux_surfaces(e2, fw.r, fw.z .+ dz)
    p0 = eqt.profiles_1d
    p2 = e2.profiles_1d
    for name in (:psi, :phi, :rho_tor, :rho_tor_norm, :volume, :area, :surface, :q, :elongation,
        :triangularity_upper, :triangularity_lower,
        :squareness_upper_outer, :squareness_upper_inner, :squareness_lower_outer, :squareness_lower_inner,
        :r_inboard, :r_outboard, :gm1, :gm2, :gm5, :gm8, :gm9, :b_field_max, :b_field_min, :trapped_fraction, :j_tor)
        @test isapprox(getproperty(p2, name), getproperty(p0, name); rtol=1e-6, atol=1e-9)
    end
    @test e2.boundary.outline.r ≈ eqt.boundary.outline.r
    @test e2.boundary.outline.z ≈ eqt.boundary.outline.z .+ dz
    @test e2.global_quantities.magnetic_axis.z ≈ eqt.global_quantities.magnetic_axis.z + dz
    @test e2.global_quantities.ip ≈ eqt.global_quantities.ip
    @test e2.global_quantities.volume ≈ eqt.global_quantities.volume
    @test length(e2.boundary.x_point) == length(eqt.boundary.x_point)
    @test [x.z for x in e2.boundary.x_point] ≈ [x.z for x in eqt.boundary.x_point] .+ dz
    @test [x.r for x in e2.boundary.x_point] ≈ [x.r for x in eqt.boundary.x_point]
    @test length(e2.boundary.strike_point) == length(eqt.boundary.strike_point)
    @test [s.z for s in e2.boundary.strike_point] ≈ [s.z for s in eqt.boundary.strike_point] .+ dz
    @test e2.boundary_secondary_separatrix.outline.z ≈ eqt.boundary_secondary_separatrix.outline.z .+ dz
end

@testset "flux_surfaces up/down reflection" begin
    dd, fw = fresh_sample()
    e2, fw2 = reflected_eqt_z(dd.equilibrium.time_slice[], fw)
    eqt = dd.equilibrium.time_slice[]
    IMAS.flux_surfaces(eqt, fw.r, fw.z)
    IMAS.flux_surfaces(e2, fw2.r, fw2.z)
    p0 = eqt.profiles_1d
    p2 = e2.profiles_1d
    # unchanged
    for name in (:psi, :phi, :rho_tor, :rho_tor_norm, :volume, :area, :surface, :q, :elongation,
        :r_inboard, :r_outboard, :gm1, :gm2, :gm5, :gm8, :gm9, :b_field_max, :b_field_min, :trapped_fraction, :j_tor)
        @test isapprox(getproperty(p2, name), getproperty(p0, name); rtol=1e-5, atol=1e-8)
    end
    # upper <-> lower
    @test isapprox(p2.triangularity_upper, p0.triangularity_lower; atol=1e-6)
    @test isapprox(p2.triangularity_lower, p0.triangularity_upper; atol=1e-6)
    @test isapprox(p2.squareness_upper_outer, p0.squareness_lower_outer; atol=1e-6)
    @test isapprox(p2.squareness_lower_outer, p0.squareness_upper_outer; atol=1e-6)
    @test isapprox(p2.squareness_upper_inner, p0.squareness_lower_inner; atol=1e-6)
    @test isapprox(p2.squareness_lower_inner, p0.squareness_upper_inner; atol=1e-6)
    @test e2.global_quantities.magnetic_axis.z ≈ -eqt.global_quantities.magnetic_axis.z
    @test e2.global_quantities.ip ≈ eqt.global_quantities.ip
    @test e2.global_quantities.volume ≈ eqt.global_quantities.volume
    @test length(e2.boundary.x_point) == length(eqt.boundary.x_point)
    @test isapprox([x.z for x in e2.boundary.x_point], -[x.z for x in eqt.boundary.x_point]; atol=1e-6)
    @test isapprox([x.r for x in e2.boundary.x_point], [x.r for x in eqt.boundary.x_point]; atol=1e-6)
    @test sort([s.z for s in e2.boundary.strike_point]) ≈ sort(-[s.z for s in eqt.boundary.strike_point])
end

@testset "find_psi_* translation and reflection invariance" begin
    dd, eqt, fw, PSI = sample_eqt_with_wall()
    dz = 1.7
    e2 = shifted_eqt_z(eqt, dz)
    e3, fw3 = reflected_eqt_z(eqt, fw)
    for (name, f) in (
        (:find_psi_separatrix, (e, wr, wz) -> IMAS.find_psi_separatrix(e)),
        (:find_psi_2nd_separatrix, (e, wr, wz) -> IMAS.find_psi_2nd_separatrix(e)),
        (:find_psi_max, (e, wr, wz) -> IMAS.find_psi_max(e)),
        (:find_psi_wall_omp, (e, wr, wz) -> IMAS.find_psi_wall_omp(e, wr, wz)),
        (:find_psi_tangent_omp, (e, wr, wz) -> IMAS.find_psi_tangent_omp(e, wr, wz, IMAS.ψ_interpolant(findfirst(:rectangular, e.profiles_2d))[3])),
        (:find_psi_last_diverted, (e, wr, wz) -> IMAS.find_psi_last_diverted(e, wr, wz, IMAS.ψ_interpolant(findfirst(:rectangular, e.profiles_2d))[3])),
    )
        v0 = f(eqt, fw.r, fw.z)
        v2 = f(e2, fw.r, fw.z .+ dz)
        v3 = f(e3, fw3.r, fw3.z)
        @testset "$name" begin
            @test all(isapprox.(collect(values(v2)), collect(values(v0)); rtol=1e-12, atol=1e-12))
            @test all(isapprox.(collect(values(v3)), collect(values(v0)); rtol=1e-6, atol=1e-9))
        end
    end
end

@testset "sol translation invariance" begin
    dd, fw = fresh_sample()
    dz = 1.7
    e2 = shifted_eqt_z(dd.equilibrium.time_slice[], dz)
    eqt = dd.equilibrium.time_slice[]
    IMAS.flux_surfaces(eqt, fw.r, fw.z)
    IMAS.flux_surfaces(e2, fw.r, fw.z .+ dz)
    S0 = IMAS.sol(eqt, fw.r, fw.z; levels=10)
    S2 = IMAS.sol(e2, fw.r, fw.z .+ dz; levels=10)
    for key in keys(S0)
        @test length(S2[key]) == length(S0[key])
        for (l0, l2) in zip(S0[key], S2[key])
            @test l2.psi ≈ l0.psi
            @test l2.r ≈ l0.r
            @test l2.z ≈ l0.z .+ dz
            @test l2.Bp ≈ l0.Bp
            @test l2.s ≈ l0.s
            @test l2.midplane_index == l0.midplane_index
            @test l2.strike_angles ≈ l0.strike_angles
            @test l2.grazing_angles ≈ l0.grazing_angles
            @test l2.wall_index == l0.wall_index
        end
    end
end

@testset "find_levels_from_P psi sign invariance" begin
    # psi(R) at the outboard midplane is made non-monotonic OUTSIDE the wall (as a coil close to the wall would),
    # then the SOL discretization must not depend on the sign of psi
    dd, fw = fresh_sample()
    eqt = dd.equilibrium.time_slice[]
    IMAS.flux_surfaces(eqt, fw.r, fw.z)
    e2d = findfirst(:rectangular, eqt.profiles_2d)
    ZA = eqt.global_quantities.magnetic_axis.z
    r_wall = maximum(fw.r)
    rmax = e2d.grid.dim1[end]
    @test rmax > r_wall + 0.05
    Rb = (r_wall + rmax) / 2
    w = (rmax - r_wall) / 4
    psi_wall = IMAS.find_psi_wall_omp(eqt, fw.r, fw.z)
    A = -3 * (psi_wall - eqt.boundary_separatrix.psi)
    for (i, R) in enumerate(e2d.grid.dim1), (j, Z) in enumerate(e2d.grid.dim2)
        e2d.psi[i, j] += A * exp(-((R - Rb)^2 + (Z - ZA)^2) / w^2)
    end
    _, _, PSI = IMAS.ψ_interpolant(e2d)
    psi_omp = PSI.(e2d.grid.dim1, ZA)
    @test !issorted(psi_omp) # the bump makes psi(R) turn around beyond the wall
    eN = negated_psi_eqt(eqt)
    _, _, PSIN = IMAS.ψ_interpolant(findfirst(:rectangular, eN.profiles_2d))
    r = collect(range(0.0, 0.3; length=100))
    q = 1e7 .* exp.(-r ./ 0.01)
    lev, R, p = IMAS.find_levels_from_P(eqt, fw.r, fw.z, PSI, r, q, 12)
    levN, RN, pN = IMAS.find_levels_from_P(eN, fw.r, fw.z, PSIN, r, q, 12)
    @test all(R .<= r_wall)
    @test RN ≈ R
    @test levN ≈ -lev
    @test pN ≈ p
end

end # invariance

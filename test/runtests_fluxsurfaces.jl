using IMAS
using Test

include("helpers.jl")

do_plot = false
if do_plot
    using Plots
end

@testset "flux_surfaces" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; show_warnings = false)

    dd_orig = deepcopy(dd)

    fw = IMAS.first_wall(dd.wall)
    IMAS.flux_surfaces(dd.equilibrium, fw.r, fw.z)

    ids_orig = dd_orig.equilibrium
    ids = dd.equilibrium

    diff(ids_orig, ids; tol = 1E-1)
end

@testset "fluxsurface_extrema" begin
    # explicit polyline; min_z ties at idx 1 and 5 -> first index wins (matches findmin)
    pr = [1.0, 2.0, 1.5, 0.5, 1.0]
    pz = [0.0, 0.5, 1.0, 0.3, 0.0]
    #     imaxr iminr imaxz iminz  r@maxz max_z r@minz min_z  z@maxr max_r z@minr min_r
    @test Tuple(IMAS.fluxsurface_extrema(pr, pz)) ==
          (2, 4, 3, 1, 1.5, 1.0, 1.0, 0.0, 0.5, 2.0, 0.3, 0.5)

    # equivalence with findmax/findmin on a closed D-shaped polyline
    θ = range(0, 2π; length = 257)
    pr2 = collect(1.7 .+ 0.6 .* cos.(θ .+ 0.3 .* sin.(θ)))
    pz2 = collect(1.8 .* sin.(θ))
    pr2[end], pz2[end] = pr2[1], pz2[1]
    g = IMAS.fluxsurface_extrema(pr2, pz2)
    @test (g[1], g[2], g[3], g[4]) == (findmax(pr2)[2], findmin(pr2)[2], findmax(pz2)[2], findmin(pz2)[2])
    @test (g[6], g[8], g[10], g[12]) == (maximum(pz2), minimum(pz2), maximum(pr2), minimum(pr2))

    @test_throws DimensionMismatch IMAS.fluxsurface_extrema([1.0, 2.0], [0.0])
end

@testset "luce_squareness quadrant labels" begin
    # four DISTINCT squareness values so that a label swap between quadrants is detected
    ζ = (ou=0.10, iu=-0.20, il=0.25, ol=-0.05)
    pr, pz = luce_boundary(; ζ)
    g = IMAS.miller_R_a_κ_δ_ζ(pr, pz)
    @test isapprox(g.ζou, ζ.ou; atol=2e-3)
    @test isapprox(g.ζol, ζ.ol; atol=2e-3)
    @test isapprox(g.ζil, ζ.il; atol=2e-3)
    @test isapprox(g.ζiu, ζ.iu; atol=2e-3)
    @test isapprox(g.δu, 0.4; atol=1e-6)
    @test isapprox(g.δl, 0.5; atol=1e-6)
end

@testset "boundary_shape zoffset" begin
    # the lower half must be mirrored about Z = zoffset, so a nonzero zoffset is a pure vertical translation
    kw = (a=0.6, eps=0.35, kapu=1.9, kapl=1.7, delu=0.4, dell=0.5, zetaou=-0.1, zetaiu=-0.3, zetail=-0.2, zetaol=-0.15)
    b0 = IMAS.boundary_shape(; kw..., zoffset=0.0)
    for zoffset in (-0.05, 0.12)
        b = IMAS.boundary_shape(; kw..., zoffset)
        @test b.r ≈ b0.r
        @test b.z ≈ b0.z .+ zoffset
        @test b.zref ≈ b0.zref .+ zoffset
        @test isapprox(maximum(b.z), zoffset + kw.a * kw.kapu; atol=1e-3)
        @test isapprox(minimum(b.z), zoffset - kw.a * kw.kapl; atol=1e-3)
        @test isapprox(b.z[argmax(b.r)], zoffset; atol=1e-3)
    end
    # same with x-point branches active (parameters from the boundary_shape docstring)
    kwx = (a=0.608, eps=0.374, kapu=1.920, kapl=1.719, delu=0.769, dell=0.463, zetaou=-0.155, zetaiu=-0.255, zetail=-0.174, zetaol=-0.227)
    for (upnull, lonull) in ((true, false), (false, true), (true, true))
        b0 = IMAS.boundary_shape(; kwx..., zoffset=0.0, upnull, lonull)
        b = IMAS.boundary_shape(; kwx..., zoffset=0.1, upnull, lonull)
        @test b.r ≈ b0.r
        @test b.z ≈ b0.z .+ 0.1
    end
end

# ---------------------------------------------------------------------------------------------------------------------
# Label / axis / argument-order checks on deliberately asymmetric cases (axis off the grid centre, nonzero vertical
# offsets), so that a swapped label or a mirror about the wrong axis changes the answer.


@testset "find_psi_separatrix labels" begin
    dd, eqt, fw, PSI = sample_eqt_with_wall()
    s = IMAS.find_psi_separatrix(eqt)
    isclosed(psi) = (sf = IMAS.flux_surface(eqt, psi, :encircling, Float64[], Float64[]); length(sf) == 1 && sf[1].r[1] == sf[1].r[end] && sf[1].z[1] == sf[1].z[end])
    psi_sign = sign(eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1])
    @test isclosed(s.psi_sep_closed)
    @test !isclosed(s.psi_sep_open)
    @test psi_sign * (s.psi_sep_open - s.psi_sep_closed) > 0
end

@testset "find_psi_last_diverted limited branch" begin
    dd, eqt, fw, PSI = sample_eqt_with_wall()
    eqt_lim = deepcopy(eqt)
    eqt_lim.boundary.psi = eqt_lim.boundary_separatrix.psi   # forces the limited-plasma branch
    res = IMAS.find_psi_last_diverted(eqt_lim, fw.r, fw.z, PSI)
    @test isfinite(res.psi_last_lfs) && isfinite(res.psi_first_lfs_far)
end

@testset "find_psi_last_diverted consumers" begin
    dd, eqt, fw, PSI = sample_eqt_with_wall()
    res = IMAS.find_psi_last_diverted(eqt, fw.r, fw.z, PSI)
    @test keys(res) == (:psi_last_lfs, :psi_first_open, :psi_first_lfs_far, :null_within_wall)
    @test res.null_within_wall isa Bool
    # particle_heat_flux unpacks this result; it must run and return finite heat fluxes
    # (explicit q(r) at the OMP, since the equilibrium-only sample has no power data for the default Eich width)
    r = collect(range(0.0, 0.4; length=200))
    q = 1e7 .* exp.(-r ./ 0.01)
    Rwall, Zwall, s, SOL, r, q = IMAS.mesher_heat_flux(dd; r, q)
    q_part, q_par = IMAS.particle_heat_flux(SOL, Rwall, Zwall, r, q)
    @test length(q_part) == length(Rwall)
    @test all(isfinite, q_part) && all(isfinite, q_par)
end

@testset "find_x_point! translation invariance" begin
    dd, eqt, fw, PSI = sample_eqt_with_wall()
    xp0 = [(x.r, x.z) for x in eqt.boundary.x_point]
    @test length(xp0) == 2
    # shifting everything vertically must shift the x-points and nothing else (fails when Z=0 is used instead of Z0)
    for dz in (0.5, 1.5, 2.0, 3.0)
        e2 = deepcopy(eqt)
        e2d = findfirst(:rectangular, e2.profiles_2d)
        e2d.grid.dim2 = e2d.grid.dim2 .+ dz
        e2.global_quantities.magnetic_axis.z += dz
        e2.boundary.outline.z = e2.boundary.outline.z .+ dz
        for sp in e2.boundary.strike_point
            sp.z += dz
        end
        IMAS.find_x_point!(e2, fw.r, fw.z .+ dz)
        xp = [(x.r, x.z - dz) for x in e2.boundary.x_point]
        @test length(xp) == length(xp0)
        if length(xp) == length(xp0)
            @test all(isapprox(a[1], b[1]; atol=1e-6) && isapprox(a[2], b[2]; atol=1e-6) for (a, b) in zip(xp, xp0))
        end
    end
end

@testset "find_levels_from_P method consistency" begin
    dd, eqt, fw, PSI = sample_eqt_with_wall()
    r = collect(range(0.0, 0.05; length=50))
    q = exp.(-r ./ 0.01)
    a = IMAS.find_levels_from_P(eqt, fw.r, fw.z, PSI, r, q, 10)
    b = IMAS.find_levels_from_P(eqt, dd.wall, PSI, r, q, 10)
    c = IMAS.find_levels_from_P(dd, r, q, 10)
    @test a[1] ≈ b[1]
    @test a[1] ≈ c[1]
end

# shift an equilibrium (and its derived points) vertically by dz, for translation-invariance tests

@testset "find_psi_tangent_omp / find_psi_max translation invariance" begin
    dd, eqt, fw, PSI = sample_eqt_with_wall()
    t0 = IMAS.find_psi_tangent_omp(eqt, fw.r, fw.z, PSI)
    m0 = IMAS.find_psi_max(eqt)
    @test isfinite(t0.psi_tangent_in) && isfinite(t0.psi_tangent_out) && isfinite(m0)
    for dz in (1.5, 3.0)
        e2 = shifted_eqt_z(eqt, dz)
        _, _, PSI2 = IMAS.ψ_interpolant(e2.profiles_2d)
        t = IMAS.find_psi_tangent_omp(e2, fw.r, fw.z .+ dz, PSI2)
        @test isapprox(t.psi_tangent_in, t0.psi_tangent_in; rtol=1e-6)
        @test isapprox(t.psi_tangent_out, t0.psi_tangent_out; rtol=1e-6)
        @test isapprox(IMAS.find_psi_max(e2), m0; rtol=1e-6)
    end
end

@testset "interp_rmid_at_psi monotonic truncation" begin
    R = collect(range(2.0, 2.5; length=51))
    # monotonic: every sample is used, including the last one
    f_quad = (r, z) -> (r - 2.0)^2
    itp = IMAS.interp_rmid_at_psi(f_quad, R, 0.0)
    @test isapprox(itp(0.5^2), 2.5; atol=1e-9)
    # increasing, then decreasing, then increasing again: only the first monotonic stretch may be used
    f_updown = (r, z) -> r < 2.2 ? (r - 2.0) : (r < 2.3 ? 0.2 - (r - 2.2) : 0.1 + (r - 2.3))
    itp2 = IMAS.interp_rmid_at_psi(f_updown, R, 0.0)
    @test isapprox(itp2(0.15), 2.15; atol=1e-9)
    @test isapprox(itp2(0.05), 2.05; atol=1e-9)
end

@testset "boundary_shape quadrant end points" begin
    kw = (a=0.6, eps=0.35, kapu=1.9, kapl=1.7, delu=0.4, dell=0.5, zetaou=-0.1, zetaiu=-0.3, zetail=-0.2, zetaol=-0.15)
    npts = 90
    b = IMAS.boundary_shape(; kw..., zoffset=0.0, npts)
    R0 = kw.a / kw.eps
    @test length(b.r) == 4 * npts + 1
    @test b.r[1] == b.r[end] && b.z[1] == b.z[end]                               # closed
    @test length(unique(collect(zip(b.r, b.z)))) == length(b.r) - 1              # no duplicate points besides the closing one
    @test isapprox(minimum(b.r), R0 - kw.a; atol=1e-12)                          # inner midplane point is generated
    @test isapprox(maximum(b.r), R0 + kw.a; atol=1e-12)
    @test isapprox(maximum(b.z), kw.a * kw.kapu; atol=1e-12)                     # top and bottom points are generated
    @test isapprox(minimum(b.z), -kw.a * kw.kapl; atol=1e-12)
    @test isapprox(b.r[argmax(b.z)], R0 - kw.a * kw.delu; atol=1e-12)
    @test isapprox(b.r[argmin(b.z)], R0 - kw.a * kw.dell; atol=1e-12)
    # round trip through the shape analysis recovers the input parameters
    g = IMAS.miller_R_a_κ_δ_ζ(b.r, b.z)
    @test isapprox(g.R0, R0; atol=1e-12) && isapprox(g.a, kw.a; atol=1e-12)
    @test isapprox(g.δu, kw.delu; atol=1e-12) && isapprox(g.δl, kw.dell; atol=1e-12)
    for (got, want) in ((g.ζou, kw.zetaou), (g.ζol, kw.zetaol), (g.ζiu, kw.zetaiu), (g.ζil, kw.zetail))
        @test isapprox(got, want; atol=2e-3)
    end
end

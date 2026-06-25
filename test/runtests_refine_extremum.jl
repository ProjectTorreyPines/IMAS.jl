using IMAS
using Test
import Interpolations

# _refine_extremum solves {ψ = target, ∂ψ/∂(diraxis) = 0} with a 2×2 Newton
# seeded at a rough extremum (falling back to the seed if it cannot converge),
# then validates the result and, near an X-point, recovers the genuine confined
# extremum using the magnetic axis passed as the last argument.
# Analytic test field: nested ellipses ψ = ((R-R0)/a0)^2 + (Z/b0)^2, whose
# level-L surface has extrema  max_r=(R0+a0√L, 0),  min_r=(R0-a0√L, 0),
#                              max_z=(R0, b0√L),    min_z=(R0, -b0√L).
@testset "_refine_extremum" begin
    R0, a0, b0 = 1.7, 0.5, 1.0
    psi_fun = (R, Z) -> ((R - R0) / a0)^2 + (Z / b0)^2
    r = range(1.0, 2.4, length=65)
    z = range(-1.3, 1.3, length=65)
    dc = step(r)
    PSI = [psi_fun(ri, zj) for ri in r, zj in z]
    itp = IMAS.ψ_interpolant(r, z, PSI).PSI_interpolant

    @testset "recovers all four extrema from a rough seed" begin
        L = 0.36
        s = sqrt(L)
        # (label, extremum_of, truth, seed ~half a cell off toward the interior)
        cases = (
            ("max_r", :R, (R0 + a0 * s, 0.0), (R0 + a0 * s - 0.4dc, 0.3dc)),
            ("min_r", :R, (R0 - a0 * s, 0.0), (R0 - a0 * s + 0.4dc, 0.3dc)),
            ("max_z", :Z, (R0, b0 * s), (R0 + 0.3dc, b0 * s - 0.4dc)),
            ("min_z", :Z, (R0, -b0 * s), (R0 + 0.3dc, -b0 * s + 0.4dc)),
        )
        for (label, extremum_of, truth, seed) in cases
            R, Z = IMAS._refine_extremum(itp, L, seed, extremum_of, (R0, 0.0))
            @test isapprox(R, truth[1]; atol=1e-6)
            @test isapprox(Z, truth[2]; atol=1e-6)
        end
    end

    @testset "rejects an invalid extremum_of" begin
        @test_throws ArgumentError IMAS._refine_extremum(itp, 0.36, (R0, 0.0), :bogus, (R0, 0.0))
    end

    @testset "falls back to seed when Newton cannot converge" begin
        # seeding exactly at the magnetic axis gives ∇ψ = 0 -> singular Jacobian.
        # The solver must not blow up to NaN/Inf; it returns the seed unchanged.
        seed = (R0, 0.0)
        R, Z = IMAS._refine_extremum(itp, 0.36, seed, :R, (R0, 0.0))
        @test (R, Z) == seed
    end

    @testset "refinement keeps max_z below X-point (DIII-D)" begin
        # On a diverted equilibrium the outermost (~separatrix) surface's rough
        # max_z overshoots ABOVE the upper X-point. The plain fast-path Newton
        # seeded there converges to the wrong basin (a solution of {ψ=target,
        # ∂ψ/∂R=0} above the X-point) instead of the confined-region top below it.
        # Given the magnetic axis, _refine_extremum detects the wrong branch
        # (curvature check), finds the X-point by ∇ψ=0, mirrors across it, and
        # re-solves below it -> the genuine confined max below the X-point.
        filename = joinpath(pkgdir(IMAS.IMASdd), "sample", "D3D_eq_ods.json")
        dd = IMAS.json2imas(filename; show_warnings=false)
        eqt = dd.equilibrium.time_slice[1]
        fw = IMAS.first_wall(dd.wall)
        eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
        rr, zz, itp2 = IMAS.ψ_interpolant(eqt2d)
        eqt1d = eqt.profiles_1d
        RA, ZA = eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z
        psi_axis = itp2(RA, ZA)
        pb = IMAS.find_psi_boundary(rr, zz, eqt2d.psi, psi_axis, eqt1d.psi[end], RA, ZA, fw.r, fw.z;
            PSI_interpolant=itp2, raise_error_on_not_open=false, raise_error_on_not_closed=false)
        psis = (eqt1d.psi .- eqt1d.psi[1]) ./ (eqt1d.psi[end] - eqt1d.psi[1]) .* (pb.last_closed - psi_axis) .+ psi_axis
        rough = IMAS.trace_surfaces(psis, eqt1d.f, collect(rr), collect(zz), eqt2d.psi, itp2, RA, ZA, fw.r, fw.z; refine_extrema=false)
        s = rough[end]
        xp_up = maximum(xp.z for xp in eqt.boundary.x_point)
        _, max_z = IMAS._refine_extremum(itp2, psis[end], (s.r_at_max_z, s.max_z), :Z, (RA, ZA))
        @test max_z < xp_up
    end
end

# The Newton refine is backend-agnostic via the _psi_* adapters: it must work not only with a
# FastInterpolations interpolant but also with an Interpolations.jl one (e.g. as FRESCO builds).
@testset "_refine_extremum backend-agnostic (Interpolations.jl)" begin
    @testset "analytic ellipse, IP interpolant" begin
        R0, a0, b0 = 1.7, 0.5, 1.0
        psi_fun = (R, Z) -> ((R - R0) / a0)^2 + (Z / b0)^2
        r = range(1.0, 2.4, length=65)
        z = range(-1.3, 1.3, length=65)
        PSI = [psi_fun(ri, zj) for ri in r, zj in z]
        itp = Interpolations.cubic_spline_interpolation((r, z), PSI; extrapolation_bc=Interpolations.Line())
        @test !(itp isa IMAS.FI.AbstractInterpolant)   # genuinely exercises the non-FI path
        L = 0.36; s = sqrt(L); dc = step(r)
        cases = (
            (:R, (R0 + a0 * s, 0.0), (R0 + a0 * s - 0.4dc, 0.3dc)),
            (:R, (R0 - a0 * s, 0.0), (R0 - a0 * s + 0.4dc, 0.3dc)),
            (:Z, (R0, b0 * s), (R0 + 0.3dc, b0 * s - 0.4dc)),
            (:Z, (R0, -b0 * s), (R0 + 0.3dc, -b0 * s + 0.4dc)),
        )
        for (extremum_of, truth, seed) in cases
            R, Z = IMAS._refine_extremum(itp, L, seed, extremum_of, (R0, 0.0))
            @test isapprox(R, truth[1]; atol=1e-3)
            @test isapprox(Z, truth[2]; atol=1e-3)
        end
    end

    @testset "DIII-D trace_surfaces: Interpolations matches FastInterpolations" begin
        filename = joinpath(pkgdir(IMAS.IMASdd), "sample", "D3D_eq_ods.json")
        dd = IMAS.json2imas(filename; show_warnings=false)
        eqt = dd.equilibrium.time_slice[1]
        fw = IMAS.first_wall(dd.wall)
        eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
        rr, zz, itp_fi = IMAS.ψ_interpolant(eqt2d)
        RA, ZA = eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z
        psi_axis = itp_fi(RA, ZA)
        eqt1d = eqt.profiles_1d
        fracs = collect(range(0.1, 0.95; length=16))
        psis = [psi_axis; psi_axis .+ fracs .* (eqt1d.psi[end] - psi_axis)]
        ff = fill(eqt1d.f[end], length(psis))
        # FRESCO-style: an Interpolations.jl cubic interpolant of the SAME ψ grid as itp_fi
        itp_ip = Interpolations.cubic_spline_interpolation((rr, zz), collect(eqt2d.psi); extrapolation_bc=Interpolations.Line())
        @test !(itp_ip isa IMAS.FI.AbstractInterpolant)
        r = collect(rr); z = collect(zz); PSI = collect(eqt2d.psi)
        # identical marching trace (same PSI matrix); only the refine backend differs
        s_fi = IMAS.trace_surfaces(psis, ff, r, z, PSI, itp_fi, RA, ZA, fw.r, fw.z; refine_extrema=true)
        s_ip = IMAS.trace_surfaces(psis, ff, r, z, PSI, itp_ip, RA, ZA, fw.r, fw.z; refine_extrema=true)
        @test length(s_ip) == length(s_fi)
        for k in 2:length(s_fi)
            @test isapprox(s_ip[k].max_r, s_fi[k].max_r; atol=2e-3)
            @test isapprox(s_ip[k].min_r, s_fi[k].min_r; atol=2e-3)
            @test isapprox(s_ip[k].max_z, s_fi[k].max_z; atol=2e-3)
            @test isapprox(s_ip[k].min_z, s_fi[k].min_z; atol=2e-3)
        end
    end
end

# Production tracing uses _robust_refine_extremum!: from ANY seed on the correct side of the
# axis (even one outside the separatrix or in the private flux region across an X-point) it must
# still arrive at the correct CONFINED extremum, never an off-surface / private / SOL point.
# This is what guards the KDEMO (a_eq<0) and MANTA (elongation≈0 -> sqrt DomainError) failures.
@testset "_robust_refine_extremum arrives from arbitrary/bad seeds" begin
    @testset "far correct-side seed (analytic ellipse, no X-point)" begin
        R0, a0, b0 = 1.7, 0.5, 1.0
        psi_fun = (R, Z) -> ((R - R0) / a0)^2 + (Z / b0)^2
        r = range(1.0, 2.4, length=65)
        z = range(-1.3, 1.3, length=65)
        itp = IMAS.ψ_interpolant(r, z, [psi_fun(ri, zj) for ri in r, zj in z]).PSI_interpolant
        L = 0.36; s = sqrt(L)
        # seed far OUTSIDE the L surface (ψ≈1.5≫L) but on the outboard side -> still finds max_r
        R, Z = IMAS._robust_refine_extremum(itp, L, (R0 + 0.6, 0.2), :R, (R0, 0.0))
        @test isapprox(R, R0 + a0 * s; atol=1e-4)
        @test isapprox(Z, 0.0; atol=1e-4)
        # seed far above the L surface on the upper side -> still finds max_z
        R, Z = IMAS._robust_refine_extremum(itp, L, (R0 + 0.1, 1.2), :Z, (R0, 0.0))
        @test isapprox(R, R0; atol=1e-4)
        @test isapprox(Z, b0 * s; atol=1e-4)
    end

    @testset "outside-separatrix and private-region seeds (DIII-D)" begin
        filename = joinpath(pkgdir(IMAS.IMASdd), "sample", "D3D_eq_ods.json")
        dd = IMAS.json2imas(filename; show_warnings=false)
        eqt = dd.equilibrium.time_slice[1]
        fw = IMAS.first_wall(dd.wall)
        eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
        rr, zz, itp = IMAS.ψ_interpolant(eqt2d)
        eqt1d = eqt.profiles_1d
        RA, ZA = eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z
        axis = (RA, ZA)
        psi_axis = itp(RA, ZA)
        pb = IMAS.find_psi_boundary(rr, zz, eqt2d.psi, psi_axis, eqt1d.psi[end], RA, ZA, fw.r, fw.z;
            PSI_interpolant=itp, raise_error_on_not_open=false, raise_error_on_not_closed=false)
        psis = (eqt1d.psi .- eqt1d.psi[1]) ./ (eqt1d.psi[end] - eqt1d.psi[1]) .* (pb.last_closed - psi_axis) .+ psi_axis
        rough = IMAS.trace_surfaces(psis, eqt1d.f, collect(rr), collect(zz), eqt2d.psi, itp, RA, ZA, fw.r, fw.z; refine_extrema=false)
        xpoints = [(xp.r, xp.z) for xp in eqt.boundary.x_point]
        xp_up = eqt.boundary.x_point[argmax(xp.z for xp in eqt.boundary.x_point)]
        dr, dz = step(rr), step(zz)
        lo = (first(rr), first(zz)); hi = (last(rr), last(zz))   # clamp to ψ grid box (as production does)
        N = length(rough)
        psiN(p) = (itp(p[1], p[2]) - psi_axis) / (pb.last_closed - psi_axis)
        for k in (N, N - 2, N - 5)
            s = rough[k]
            c = psis[k]
            good_r = IMAS._robust_refine_extremum(itp, c, (s.max_r, s.z_at_max_r), :R, axis, xpoints; lo, hi)
            good_z = IMAS._robust_refine_extremum(itp, c, (s.r_at_max_z, s.max_z), :Z, axis, xpoints; lo, hi)
            # bad seed A: outside the separatrix (ψ_N>1), far outboard
            outside = (s.max_r + 12dr, s.z_at_max_r)
            @test psiN(outside) > 1
            rA = IMAS._robust_refine_extremum(itp, c, outside, :R, axis, xpoints; lo, hi)
            @test isapprox(rA[1], good_r[1]; atol=5dr)
            @test psiN(rA) <= 1 + 1e-6
            # bad seed B: private region (reflect the rough max_z across the upper X-point)
            priv = (2xp_up.r - s.r_at_max_z, 2xp_up.z - s.max_z)
            rB = IMAS._robust_refine_extremum(itp, c, priv, :Z, axis, xpoints; lo, hi)
            @test isapprox(rB[2], good_z[2]; atol=5dz)
            @test rB[2] < xp_up.z
        end
    end

    # Hardest case: at ψ_N=0.999 the PRIVATE region across the upper X-point ALSO has a
    # ψ_N=0.999 point with ∂ψ/∂R=0 — a genuine solution of the SAME {ψ=c, ∂ψ/∂R=0} system.
    # Seeds at / above the X-point must NOT be trapped there; they must return the confined
    # max_z below the X-point. (This is the geometry behind the MANTA elongation≈0 failure.)
    @testset "not trapped in the private lobe at ψ_N=0.999 (DIII-D)" begin
        filename = joinpath(pkgdir(IMAS.IMASdd), "sample", "D3D_eq_ods.json")
        dd = IMAS.json2imas(filename; show_warnings=false)
        eqt = dd.equilibrium.time_slice[1]
        fw = IMAS.first_wall(dd.wall)
        eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
        rr, zz, itp = IMAS.ψ_interpolant(eqt2d)
        eqt1d = eqt.profiles_1d
        RA, ZA = eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z
        axis = (RA, ZA)
        psi_axis = itp(RA, ZA)
        pb = IMAS.find_psi_boundary(rr, zz, eqt2d.psi, psi_axis, eqt1d.psi[end], RA, ZA, fw.r, fw.z;
            PSI_interpolant=itp, raise_error_on_not_open=false, raise_error_on_not_closed=false)
        xpoints = [(xp.r, xp.z) for xp in eqt.boundary.x_point]
        xp_up = eqt.boundary.x_point[argmax(xp.z for xp in eqt.boundary.x_point)]
        lo = (first(rr), first(zz)); hi = (last(rr), last(zz))
        c = psi_axis + 0.999 * (pb.last_closed - psi_axis)   # ψ_N = 0.999, just inside the separatrix
        psiN(p) = (itp(p[1], p[2]) - psi_axis) / (pb.last_closed - psi_axis)
        # confined truth: seed just below the upper X-point
        truth = IMAS._robust_refine_extremum(itp, c, (xp_up.r, xp_up.z - 0.15), :Z, axis, xpoints; lo, hi)
        @test truth[2] < xp_up.z
        @test isapprox(psiN(truth), 0.999; atol=2e-3)
        # seeds AT / ABOVE the upper X-point (in / near the private lobe) must come back confined
        for sd in ((xp_up.r, xp_up.z + 0.05), (xp_up.r + 0.05, xp_up.z + 0.03),
                   (xp_up.r - 0.05, xp_up.z + 0.05), (xp_up.r, xp_up.z + 0.25))
            res = IMAS._robust_refine_extremum(itp, c, sd, :Z, axis, xpoints; lo, hi)
            @test res[2] < xp_up.z                         # not trapped above the X-point
            @test isapprox(res[2], truth[2]; atol=1e-4)    # the genuine confined max_z
        end
    end
end

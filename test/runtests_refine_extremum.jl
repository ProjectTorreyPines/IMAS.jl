using IMAS
using Test

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

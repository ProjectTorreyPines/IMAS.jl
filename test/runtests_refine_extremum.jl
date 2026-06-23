using IMAS
using Test

# _refine_extremum solves {ψ = target, ∂ψ/∂(diraxis) = 0} with a 2×2 Newton
# seeded at a rough extremum, falling back to the seed if it cannot converge.
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
            R, Z = IMAS._refine_extremum(itp, L, seed, extremum_of)
            @test isapprox(R, truth[1]; atol=1e-6)
            @test isapprox(Z, truth[2]; atol=1e-6)
        end
    end

    @testset "rejects an invalid extremum_of" begin
        @test_throws ArgumentError IMAS._refine_extremum(itp, 0.36, (R0, 0.0), :bogus)
    end

    @testset "falls back to seed when Newton cannot converge" begin
        # seeding exactly at the magnetic axis gives ∇ψ = 0 -> singular Jacobian.
        # The solver must not blow up to NaN/Inf; it returns the seed unchanged.
        seed = (R0, 0.0)
        R, Z = IMAS._refine_extremum(itp, 0.36, seed, :R)
        @test (R, Z) == seed
    end
end

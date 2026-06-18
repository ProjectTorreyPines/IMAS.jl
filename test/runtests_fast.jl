using IMAS
using Test

# Guards the ONLY Roots.find_zero call site in the package: critical_energy(...;
# approximate=false) in src/physics/fast.jl. Nothing else exercises Roots, so
# this test pins the result and keeps the widened `Roots = "2, 3"` compat honest
# (the value below was verified to match under both Roots 2.x and 3.x).
@testset "fast particles" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "omas_sample.json")
    dd = IMAS.json2imas(filename; show_warnings=false)
    cp1d = dd.core_profiles.profiles_1d[1]

    @testset "critical_energy Roots path (approximate=false)" begin
        ec_approx = IMAS.critical_energy(cp1d, 1, 2.0, 1.0; approximate=true)
        ec_roots = IMAS.critical_energy(cp1d, 1, 2.0, 1.0; approximate=false)  # Roots.find_zero

        # the root-find must run and return a sensible, distinct value
        @test isfinite(ec_roots) && ec_roots > 0
        @test !isapprox(ec_roots, ec_approx; rtol=1e-3)   # Roots path actually changed the answer

        # regression golden (verified to match under both Roots 2.x and 3.x)
        @test ec_roots ≈ 93662.596 rtol = 1e-5
    end
end

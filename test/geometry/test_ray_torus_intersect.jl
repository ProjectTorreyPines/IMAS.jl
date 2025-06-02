using IMAS
using Test

@testset "Test ray_torus_intersect" begin
    test_cases = [
        ([6.0,0.0,0.0], [-1.0,0.0,0.0], [1.0, 4.0], [-1.0, 1.0], [4.0 0.0 0.0], )
    ]
    for (origin, direction, r_bounds, z_bounds, expected) in test_cases
        @testset "origin=$origin, direction=$direction, r_bounds=$r_bounds, z_bounds=$z_bounds" begin
            # println("--- x0=$x0, y0=$y0, dx=$dx, dy=$dy, r_min=$r_min, r_max=$r_max ---")
            test_result = IMAS.ray_torus_intersect(origin, direction, r_bounds, z_bounds)
            println("--- Result | expected: $test_result | $expected  ---")
            @test test_result ≈ expected
        end
    end
end
using IMAS
using Test
@testset "Test toroical_intersection" begin
    test_cases = [
        ([1.0, 4.0, 4.0, 1.0, 1.0], [-1.0, -1.0, 1.0, 1.0, -1.0], [6.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [4.0,0.0,0.0])
    ]
    for (r_box, z_box, origin, direction, expected) in test_cases
        @testset "origin=$origin, direction=$direction, r_box=$r_box, z_box=$z_box" begin
            # println("--- x0=$x0, y0=$y0, dx=$dx, dy=$dy, r_min=$r_min, r_max=$r_max ---")
            t = IMAS.toroidal_intersection(r_box, z_box, origin, direction)
            test_result = origin .+ direction .* t
            println("--- Result | expected: $test_result | $expected  ---")
            @test test_result ≈ expected
        end
    end
end
using IMAS
using Test

@testset "Test solve_r_intersect" begin
    test_cases = [
        (3.0, 0.0, 0.0, 0.0, 1.0, 4.0, [[0.0 Inf];]), # ray is not moving, but starts on the grid
        (6.0, 0.0, 0.0, 0.0, 1.0, 4.0, zeros(Float64, 0)), # ray is not moving, and starts off the grid
        (6.0, 0.0, -1.0, 0.0, 1.0, 4.0, [2. 5.; 7. 10.]), # ray in x-direction 
        (0.0, 6.0, 0.0, -1.0, 2.0, 4.0, [2. 4.; 8. 10.]), # ray in y-direction
        (3.0, 0.0, -1.0, 0.0, 1.0, 4.0, [0.0 2.; 4.0 7.0]), # ray in x-direction with start in torus
        (6.0, 6.0, 0.0, -1.0, 2.0, 4.0, zeros(Float64, 0)), # ray in x-direction with offset causing no intersection
        (3.0, 6.0, -3.0/5.0, -4.0/5.0, 1.0, 2.0, [[5.0 41.0/5.0];]), # Single crossing with x-y ray
        (3.0, 8.0, -3.0/5.0, -4.0/5.0, 3.0, 4.0, [5.0 32.0/5.0; 10.0 57.0/5.0]) # double crossing with x-y ray
    ]

    for (x0, y0, dx, dy, r_min, r_max, expected) in test_cases
        @testset "x0=$x0, y0=$y0, dx=$dx, dy=$dy, r_min=$r_min, r_max=$r_max" begin
            # println("--- x0=$x0, y0=$y0, dx=$dx, dy=$dy, r_min=$r_min, r_max=$r_max ---")
            test_result = IMAS.solve_r_intersect(x0, y0, dx, dy, r_min, r_max)
            println("--- Result | expected: $test_result | $expected  ---")
            @test test_result ≈ expected
        end
    end
end
using IMAS
using Test

@testset "curvature" begin
    t = range(0, 2π; length=101)
    pr = 1.0 .+ 0.5 .* cos.(t)
    pz = 0.5 .* sin.(t)
    pr[end] = pr[1]
    pz[end] = pz[1]
    @test IMAS.is_closed_polygon(pr, pz)
    c = IMAS.curvature(pr, pz)
    @test all(c .> 0)                                  # same sign everywhere, including both endpoints
    @test maximum(c) - minimum(c) < 1e-10               # constant on a uniformly sampled circle
    co = IMAS.curvature(pr[1:50], pz[1:50])             # open polyline: interior values unchanged
    @test co[2:end-1] ≈ c[2:49]
end

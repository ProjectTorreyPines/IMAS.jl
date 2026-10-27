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

@testset "intersection crossing order" begin
    # one long segment of path 1 crosses a zigzag path 2 three times; the zigzag is listed in the opposite direction
    l1x = [0.0, 10.0]
    l1y = [0.0, 0.0]
    l2x = [8.0, 7.0, 5.0, 3.0, 2.0]
    l2y = [1.0, -1.0, 1.0, -1.0, 1.0]
    res = IMAS.intersection(l1x, l1y, l2x, l2y)
    @test length(res.crossings) == 4
    @test issorted([c[1] for c in res.crossings])                  # along path 1
    @test [ix[1] for ix in res.indexes] == [1, 1, 1, 1]
    @test [ix[2] for ix in res.indexes] == [4, 3, 2, 1]
    # reversing path 1 reverses the crossing order
    rev = IMAS.intersection(reverse(l1x), reverse(l1y), l2x, l2y)
    @test [c[1] for c in rev.crossings] == reverse([c[1] for c in res.crossings])
end

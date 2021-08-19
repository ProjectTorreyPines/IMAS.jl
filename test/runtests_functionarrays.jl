@testset "functionarrays" begin
    # begin adding parameters to simple dictionary
    x = Vector{Float64}(collect(1:4))
    a = FUSE.NumericalFDVector(x, x)
    b = FUSE.AnalyticalFDVector(x, y -> y)
    @test a .* b == a.^2
end
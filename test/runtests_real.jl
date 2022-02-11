using Revise
using IMAS
using Test

@testset "Measurements" begin
    @test force_float(1.0 ± 2.0) == 1.0
end

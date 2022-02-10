using Revise
using IMAS
using Test

@testset "Measurements" begin
    @test no_error(1.0 ± 2.0) == 1.0
end

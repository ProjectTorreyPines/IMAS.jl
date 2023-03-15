using IMAS
using Test

@testset "Extract" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASDD))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; verbose = false)

    @show(IMAS.extract(dd))
end

using IMAS
import IMASDD
using Test

@testset "digest" begin
    filename = joinpath(dirname(dirname(pathof(IMASDD))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; verbose = false)

    igest=IMAS.digest(dd)
    println(igest)
end
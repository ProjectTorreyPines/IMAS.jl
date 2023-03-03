using IMAS
import IMASDD
using Test

@testset "digest" begin
    filename = joinpath(dirname(dirname(pathof(IMASDD))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; verbose=false)

    dd.global_time = 0.0
    dd.summary.time = [0.0]
    dd.summary.line_average.n_e.value = [0.0]
    dd.summary.local.magnetic_axis.n_e.value = [0.0]
    igest = IMAS.digest(dd.summary)

    igest = IMAS.digest(dd.equilibrium)
end
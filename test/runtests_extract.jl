using IMAS
using Test

@testset "Extract" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; show_warnings = false)

    IMAS.extract(dd)
end

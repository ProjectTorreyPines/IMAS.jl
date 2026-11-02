using IMAS
using Test

@testset "symmetrize_equilibrium!" begin
    samp = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")
    dd = IMAS.json2imas(joinpath(samp, "D3D_eq_ods.json"); show_warnings=false)
    eqt = dd.equilibrium.time_slice[]
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI = IMAS.ψ_interpolant(eqt2d)
    psi_sign = sign(eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1])
    RA, ZA = IMAS.find_magnetic_axis(r, z, PSI, psi_sign)
    @test abs(ZA - (z[1] + z[end]) / 2) > 1e-3     # axis is off the grid centre, so the mirror axis matters

    IMAS.symmetrize_equilibrium!(eqt)

    r2, z2, PSI2 = IMAS.ψ_interpolant(eqt2d)
    @test isapprox((z2[1] + z2[end]) / 2, ZA; atol=1e-12)   # new grid is centred on the magnetic axis
    @test minimum(z2) >= minimum(z) && maximum(z2) <= maximum(z)   # and stays inside the original grid
    RA2, ZA2 = IMAS.find_magnetic_axis(r2, z2, PSI2, psi_sign)
    @test isapprox(ZA2, ZA; atol=1e-3)                      # the axis does not move
    for d in (0.2, 0.5)
        @test isapprox(PSI2(RA2, ZA + d), PSI2(RA2, ZA - d); rtol=1e-8)   # up-down symmetric about ZA
    end
end

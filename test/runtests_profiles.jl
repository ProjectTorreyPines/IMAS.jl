using IMAS
using Test

@testset "ion thermal pressure" begin
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d, 1)
    cp1d = dd.core_profiles.profiles_1d[1]
    cp1d.time = 0.0
    dd.global_time = 0.0
    rho = collect(range(0.0, 1.0, 11))
    cp1d.grid.rho_tor_norm = rho
    T = 1000.0 .* (1.0 .- 0.5 .* rho)
    n = 1E19 .* (1.0 .- 0.5 .* rho .^ 2)
    cp1d.electrons.temperature = T
    cp1d.electrons.density_thermal = 4.0 .* n
    resize!(cp1d.ion, 4)
    for ion in cp1d.ion
        ion.temperature = T
    end
    # 1: total density stored, thermal density derived from it (state left by quasi-neutrality on an ion without fast density)
    cp1d.ion[1].density = 2.0 .* n
    # 2: thermal density stored
    cp1d.ion[2].density_thermal = 0.1 .* n
    # 3: total and fast densities stored, thermal density is their difference
    cp1d.ion[3].density = 1.0 .* n
    cp1d.ion[3].density_fast = 0.25 .* n
    # 4: only a fast density: no thermal pressure
    cp1d.ion[4].density_fast = 0.5 .* n

    @test !IMAS.hasdata(cp1d.ion[1], :density_thermal)
    @test cp1d.ion[1].density_thermal == 2.0 .* n
    expected = IMAS.mks.e .* T .* n .* (2.0 + 0.1 + 0.75)
    @test IMAS.pressure_thermal(cp1d.ion) ≈ expected
    @test cp1d.pressure_ion_total ≈ expected
    @test cp1d.pressure_thermal ≈ expected .+ IMAS.mks.e .* T .* 4.0 .* n
    # the sum over ions equals the sum of the per-ion pressures
    @test IMAS.pressure_thermal(cp1d.ion) ≈ sum(IMAS.pressure_thermal(ion) for ion in cp1d.ion[1:3])
end

using Revise
using IMAS
using Test

@testset "expressions" begin
    ne0 = 1E20
    Te0 = 1E3
    pe0 = ne0 * Te0 * 1.60218e-19

    # here we test expressions starting from different heights in the data dictionary
    # also, we are mixing data, user-defined functions, and expressions and using
    # expressions that return vectors and scalars

    #= ========= =#

    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d, 1)
    profiles_1d = dd.core_profiles.profiles_1d[1]

    profiles_1d.grid.rho_tor_norm = range(0.0, 1.0, length = 21)
    profiles_1d.electrons.density = ne0 .* (1.0 .- profiles_1d.grid.rho_tor_norm .^ 2)
    profiles_1d.electrons.temperature = (x; _...) -> Te0 .* (1.0 .- x .^ 2)

    @test profiles_1d.electrons.pressure[1] ≈ pe0

    #= ========= =#

    core_profiles = IMAS.core_profiles()
    resize!(core_profiles.profiles_1d, 1)
    profiles_1d = core_profiles.profiles_1d[1]

    profiles_1d.grid.rho_tor_norm = range(0.0, 1.0, length = 21)
    profiles_1d.electrons.density = ne0 .* (1.0 .- profiles_1d.grid.rho_tor_norm .^ 2)
    profiles_1d.electrons.pressure = (x; _...) -> pe0 .* (1.0 .- x .^ 2)

    @test profiles_1d.electrons.temperature[1] ≈ Te0

    #= ========= =#

    profiles_1d = IMAS.core_profiles__profiles_1d()

    profiles_1d.grid.rho_tor_norm = range(0.0, 1.0, length = 21)
    profiles_1d.electrons.temperature = (x; _...) -> Te0 .* (1.0 .- x .^ 2)
    profiles_1d.electrons.pressure = (x; _...) -> pe0 .* (1.0 .- x .^ 2)

    @test profiles_1d.electrons.density[1] ≈ ne0

    #= ========= =#
    # test infinite recursion

    profiles_1d = IMAS.core_profiles__profiles_1d()

    profiles_1d.grid.rho_tor_norm = range(0.0, 1.0, length = 21)
    profiles_1d.electrons.pressure = (x; _...) -> pe0 .* (1.0 .- x .^ 2)

    @test_throws Exception profiles_1d.electrons.density[1]

    #= ========= =#
    # test expressions using scalar quantities
    time_slice = IMAS.equilibrium__time_slice()
    time_slice.profiles_1d.psi = range(0.0, 1.0, length = 11)
    time_slice.profiles_1d.volume = range(0.0, 1.0, length = 11)
    time_slice.profiles_1d.pressure = 1.0 .- range(0.0, 1.0, length = 11)
    @test time_slice.global_quantities.energy_mhd ≈ 0.75

end

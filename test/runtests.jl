using Pkg
Pkg.activate("..")
Pkg.develop("FUSE")
using FUSE
using Test

@testset "basic_dd_IO" begin
    # data = FUSE.dd();
    # data.core_profiles.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length=10)
    # data.core_profiles.profiles_1d[1].electrons.density = (1.0 .- data.core_profiles.profiles_1d[1].grid.rho_tor_norm).^2
    # @test length(data.core_profiles.profiles_1d[1].grid.rho_tor_norm) == 10
    
    cp = FUSE.core_profiles()
    cp.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length=10)
    cp.profiles_1d[1].electrons.density = (1.0 .- cp.profiles_1d[1].grid.rho_tor_norm).^2
    @test length(cp.profiles_1d[1].grid.rho_tor_norm) == 10

    # resize an array of struct
    resize!(cp.profiles_1d, 2)
    @test length(cp.profiles_1d) == 2
    cp.profiles_1d[1]=cp.profiles_1d[2]

    # working with data that is not time dependent? --> only use the relevant struct
    cp1d = FUSE.core_profiles__profiles_1d()
    cp1d.grid.rho_tor_norm = range(0, 1, length=10)
    cp1d.electrons.density = (1.0 .- cp1d.grid.rho_tor_norm).^2
    @test length(cp1d.grid.rho_tor_norm) == 10

end
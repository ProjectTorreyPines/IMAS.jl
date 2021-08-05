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
    
    crp = FUSE.core_profiles()
    crp.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length=10)
    crp.profiles_1d[1].electrons.density = (1.0 .- crp.profiles_1d[1].grid.rho_tor_norm).^2
    @test length(crp.profiles_1d[1].grid.rho_tor_norm) == 10

    # resize an array of struct
    resize!(crp.profiles_1d, 2)
    @test length(crp.profiles_1d) == 2

    # deepcopy of one arrray structure element to another
    crp.profiles_1d[2] = deepcopy(crp.profiles_1d[1])
    crp.profiles_1d[2].electrons.density .*= 2.0
    @test all(crp.profiles_1d[2].grid.rho_tor_norm .== crp.profiles_1d[1].grid.rho_tor_norm)
    @test all(crp.profiles_1d[2].electrons.density .== (crp.profiles_1d[1].electrons.density * 2.0))

    # working with data that is not time dependent? --> only use the relevant struct
    crp1d = FUSE.core_profiles__profiles_1d()
    crp1d.grid.rho_tor_norm = range(0, 1, length=10)
    crp1d.electrons.density = (1.0 .- crp1d.grid.rho_tor_norm).^2
    @test length(crp1d.grid.rho_tor_norm) == 10

end
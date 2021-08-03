using Pkg
Pkg.activate("..")
Pkg.develop("FUSE")
using FUSE
using Test

@testset "basic_dd_IO" begin
    data=FUSE.dd();
    data.core_profiles.profiles_1d[1].grid.rho_tor_norm=range(0,1,length=10)
    data.core_profiles.profiles_1d[1].electrons.density=(1.0.-data.core_profiles.profiles_1d[1].grid.rho_tor_norm).^2
    @test length(data.core_profiles.profiles_1d[1].grid.rho_tor_norm)==10
end

# @testset "single" begin
#     # Single run of FUSE
#     print("hello")
# end

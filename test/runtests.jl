using Pkg
Pkg.activate("..")
using FUSE
using Test

@testset "single" begin
    # Single run of FUSE
    print("hello")
end

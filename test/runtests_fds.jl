@testset "FDS" begin
    # instantiate and populate top-level FDS
    data = FUSE.dd();
    resize!(data.core_profiles.profiles_1d, 1)
    data.core_profiles.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length=10)
    @test length(data.core_profiles.profiles_1d[1].grid.rho_tor_norm) == 10

    # try adding some data to separate core_profiles FDS
    crp = FUSE.core_profiles()
    resize!(crp.profiles_1d, 1)
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

    # reach top FDS starting at different depths
    @test data === FUSE.top(data)
    @test data === FUSE.top(data.core_profiles.profiles_1d)
    @test data === FUSE.top(data.core_profiles.profiles_1d[1])
    @test data === FUSE.top(data.core_profiles.profiles_1d[1].grid)
#    data === FUSE.top(data.core_profiles.profiles_1d[1].grid.rho_tor_norm) # this does not work yet
    @test crp1d === FUSE.top(crp1d)
    @test crp1d === FUSE.top(crp1d.grid)
#    @test crp1d === FUSE.top(crp1d.grid.rho_tor_norm) # this does not work yet

    # add structure to an array of structures
    push!(data.core_profiles.profiles_1d, crp1d)
    @test data === FUSE.top(crp1d)

    # test fail of adding data without coordinate in FDS
    data = FUSE.dd();
    resize!(data.core_profiles.profiles_1d, 1)
    @test_throws Exception data.core_profiles.profiles_1d[1].electrons.temperature = Vector{Float64}(collect(1:10))
end

@testset "FDS_IMAS" begin
    data = FUSE.dd();
    resize!(data.core_profiles.profiles_1d, 1)

    # test f2i
    @test FUSE.f2i(data.core_profiles.profiles_1d[1].grid) == "core_profiles.profiles_1d[:].grid"

    # test i2p
    @test FUSE.i2p("core_profiles.profiles_1d[1].grid") == [:core_profiles, :profiles_1d, 1, :grid]
    @test FUSE.i2p("core_profiles.profiles_1d[:].grid") == [:core_profiles, :profiles_1d, ":", :grid]

    # test p2i
    @test FUSE.p2i([:core_profiles, :profiles_1d, 1, :grid]) == "core_profiles.profiles_1d[1].grid"
    @test FUSE.p2i([:core_profiles, :profiles_1d, ":", :grid]) == "core_profiles.profiles_1d[:].grid"

    # test imas_info
    @test FUSE.imas_info("core_profiles.profiles_1d[1]") == FUSE.imas_info("core_profiles.profiles_1d[:]")
    @test FUSE.imas_info("core_profiles.profiles_1d") == FUSE.imas_info("core_profiles.profiles_1d[:]")
    @test all([haskey(FUSE.imas_info("core_profiles.profiles_1d"), k) for k in ["coordinates","data_type","full_path","documentation"]])
    @test_throws Exception FUSE.imas_info("core_profiles.does_not_exist")

    # test coordinate of a coordinate
    coord_names, coord_values = FUSE.coordinates(data.core_profiles.profiles_1d[1].grid, :rho_tor_norm)
    @test coord_names[1] == "1...N"
    @test coord_values[1] === nothing

    # test coordinate of a 1D array (with uninitialized coordinate)
    coord_names, coord_values = FUSE.coordinates(data.core_profiles.profiles_1d[1].electrons, :temperature)
    @test coord_names[1] == "core_profiles.profiles_1d[:].grid.rho_tor_norm"
    @test coord_values[1] === missing

    # test coordinate of a 1D array (with initialized coordinate)
    data.core_profiles.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length=10)
    coord_names, coord_values = FUSE.coordinates(data.core_profiles.profiles_1d[1].electrons, :temperature)
    @test coord_names[1] == "core_profiles.profiles_1d[:].grid.rho_tor_norm"
    @test coord_values[1] === data.core_profiles.profiles_1d[1].grid.rho_tor_norm
end


@testset "JSON_IO" begin
    filename = joinpath(dirname(dirname(abspath(@__FILE__))), "sample", "sample_equilibrium_ods.json")
    data  = FUSE.json2fuse(filename; verbose=false)
    @test length(data.wall.description_2d[1].limiter.unit[1].outline.r) > 0
end

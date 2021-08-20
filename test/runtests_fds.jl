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
    @test data === FUSE.top(data; stop_at_ids=false)
    @test data === FUSE.top(data.core_profiles.profiles_1d; stop_at_ids=false)
    @test data === FUSE.top(data.core_profiles.profiles_1d[1]; stop_at_ids=false)
    @test data === FUSE.top(data.core_profiles.profiles_1d[1].grid; stop_at_ids=false)
    # @test data === FUSE.top(data.core_profiles.profiles_1d[1].grid.rho_tor_norm; stop_at_ids=false) # this does not work yet

    @test_throws Exception FUSE.top(data)
    @test data.core_profiles === FUSE.top(data.core_profiles.profiles_1d)
    @test data.core_profiles === FUSE.top(data.core_profiles.profiles_1d[1])
    @test data.core_profiles === FUSE.top(data.core_profiles.profiles_1d[1].grid)
    # @test data.core_profiles === FUSE.top(data.core_profiles.profiles_1d[1].grid.rho_tor_norm) # this does not work yet

    @test crp1d === FUSE.top(crp1d; stop_at_ids=false)
    @test crp1d === FUSE.top(crp1d.grid; stop_at_ids=false)
    # @test crp1d === FUSE.top(crp1d.grid.rho_tor_norm; stop_at_ids=false) # this does not work yet
    @test crp1d === FUSE.top(crp1d)
    @test crp1d === FUSE.top(crp1d.grid)
    # @test crp1d === FUSE.top(crp1d.grid.rho_tor_norm) # this does not work yet

    # add structure to an array of structures
    push!(data.core_profiles.profiles_1d, crp1d)
    @test data.core_profiles === FUSE.top(crp1d)

    # test fail of adding data without coordinate in FDS
    data = FUSE.dd();
    resize!(data.core_profiles.profiles_1d, 1)
    @test_throws Exception data.core_profiles.profiles_1d[1].electrons.temperature = Vector{Float64}(collect(1:10))
end

@testset "FDS_IMAS" begin
    data = FUSE.dd();
    resize!(data.core_profiles.profiles_1d, 2)

    # test f2u
    @test FUSE.f2u(data.core_profiles.profiles_1d[1].grid) == "core_profiles.profiles_1d[:].grid"
    @test FUSE.f2u(:core_profiles__profiles_1d___grid) == "core_profiles.profiles_1d[:].grid"
    @test FUSE.f2u("core_profiles__profiles_1d___grid") == "core_profiles.profiles_1d[:].grid"
    @test_throws Exception FUSE.f2u("core_profiles.profiles_1d[:].grid")

    # test i2p
    @test FUSE.i2p("core_profiles.profiles_1d[1].grid") == ["core_profiles", "profiles_1d", 1, "grid"]
    @test FUSE.i2p("core_profiles.profiles_1d[:].grid") == ["core_profiles", "profiles_1d", ":", "grid"]

    # test p2i
    @test FUSE.p2i(["core_profiles", "profiles_1d", 1, "grid"]) == "core_profiles.profiles_1d[1].grid"
    @test FUSE.p2i(["core_profiles", "profiles_1d", ":", "grid"]) == "core_profiles.profiles_1d[:].grid"
    @test_throws Exception FUSE.p2i([:core_profiles, :profiles_1d, ":", :grid])


    wall = FUSE.wall()
    resize!(wall.description_2d, 1)
    resize!(wall.description_2d[1].mobile.unit, 2)
    resize!(wall.description_2d[1].mobile.unit[2].outline, 2)
    wall__description_2d = FUSE.wall__description_2d()
    resize!(wall__description_2d.mobile.unit, 2)
    resize!(wall__description_2d.mobile.unit[2].outline, 2)

    # test f2p
    @test FUSE.f2p(wall.description_2d[1].mobile.unit[2].outline[1]) == ["wall","description_2d",1,"mobile","unit",2,"outline",1]
    @test FUSE.f2p(wall__description_2d.mobile.unit[2].outline[1]) == ["wall","description_2d",0,"mobile","unit",2,"outline",1]

    # test imas_info
    @test FUSE.imas_info("core_profiles.profiles_1d[1]") == FUSE.imas_info("core_profiles.profiles_1d[:]")
    @test FUSE.imas_info("core_profiles.profiles_1d") == FUSE.imas_info("core_profiles.profiles_1d[:]")
    @test all([haskey(FUSE.imas_info("core_profiles.profiles_1d"), k) for k in ["coordinates","data_type","full_path","documentation"]])
    @test_throws Exception FUSE.imas_info("core_profiles.does_not_exist")

    # test coordinate of a coordinate
    coords = FUSE.coordinates(data.core_profiles.profiles_1d[1].grid, :rho_tor_norm)
    @test coords[:names][1] == "1...N"
    @test coords[:values][1] === nothing

    # test coordinate of a 1D array (with uninitialized coordinate)
    coords = FUSE.coordinates(data.core_profiles.profiles_1d[1].electrons, :temperature)
    @test coords[:names][1] == "core_profiles.profiles_1d[:].grid.rho_tor_norm"
    @test coords[:values][1] === missing

    # test coordinate of a 1D array (with initialized coordinate)
    data.core_profiles.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length=10)
    data.core_profiles.profiles_1d[2].grid.rho_tor_norm = range(0, 1, length=3)

    coords = FUSE.coordinates(data.core_profiles.profiles_1d[1].electrons, :temperature)
    @test coords[:names][1] == "core_profiles.profiles_1d[:].grid.rho_tor_norm"
    @test coords[:values][1] === data.core_profiles.profiles_1d[1].grid.rho_tor_norm
    @test length(coords[:values][1]) == 10

    coords = FUSE.coordinates(data.core_profiles.profiles_1d[2].electrons, :temperature)
    @test coords[:names][1] == "core_profiles.profiles_1d[:].grid.rho_tor_norm"
    @test coords[:values][1] === data.core_profiles.profiles_1d[2].grid.rho_tor_norm
    @test length(coords[:values][1]) == 3

end


@testset "JSON_IO" begin
    filename = joinpath(dirname(dirname(abspath(@__FILE__))), "sample", "sample_equilibrium_ods.json")
    data  = FUSE.json2fuse(filename; verbose=false)
    @test length(data.wall.description_2d[1].limiter.unit[1].outline.r) > 0
end

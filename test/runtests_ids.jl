import Revise
import IMAS
using Test

@testset "IDS" begin
    # instantiate and populate top-level IDS
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d, 1)
    dd.core_profiles.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length = 10)
    @test length(dd.core_profiles.profiles_1d[1].grid.rho_tor_norm) == 10

    # try adding some dd to separate core_profiles IDS
    crp = IMAS.core_profiles()
    resize!(crp.profiles_1d, 1)
    crp.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length = 10)
    crp.profiles_1d[1].electrons.density = (1.0 .- crp.profiles_1d[1].grid.rho_tor_norm) .^ 2
    @test length(crp.profiles_1d[1].grid.rho_tor_norm) == 10

    # resize an array of struct
    resize!(crp.profiles_1d, 2)
    @test length(crp.profiles_1d) == 2

    # deepcopy of one arrray structure element to another
    crp.profiles_1d[2] = deepcopy(crp.profiles_1d[1])
    crp.profiles_1d[2].electrons.density = crp.profiles_1d[2].electrons.density .* 2.0
    @test all(crp.profiles_1d[2].grid.rho_tor_norm .== crp.profiles_1d[1].grid.rho_tor_norm)
    @test all(crp.profiles_1d[2].electrons.density .== (crp.profiles_1d[1].electrons.density * 2.0))

    # working with dd that is not time dependent? --> only use the relevant struct
    crp1d = IMAS.core_profiles__profiles_1d()
    crp1d.grid.rho_tor_norm = range(0, 1, length = 10)
    crp1d.electrons.density = (1.0 .- crp1d.grid.rho_tor_norm) .^ 2
    @test length(crp1d.grid.rho_tor_norm) == 10

    # reach top IDS starting at different depths
    @test dd === IMAS.top(dd; IDS_is_absolute_top = false)
    @test dd === IMAS.top(dd.core_profiles.profiles_1d; IDS_is_absolute_top = false)
    @test dd === IMAS.top(dd.core_profiles.profiles_1d[1]; IDS_is_absolute_top = false)
    @test dd === IMAS.top(dd.core_profiles.profiles_1d[1].grid; IDS_is_absolute_top = false)
    # @test dd === IMAS.top(dd.core_profiles.profiles_1d[1].grid.rho_tor_norm; IDS_is_absolute_top=false) # this does not work yet

    # test that the top() function stops at the IDS level by default
    @test_throws Exception IMAS.top(dd)
    @test dd.core_profiles === IMAS.top(dd.core_profiles.profiles_1d)
    @test dd.core_profiles === IMAS.top(dd.core_profiles.profiles_1d[1])
    @test dd.core_profiles === IMAS.top(dd.core_profiles.profiles_1d[1].grid)
    # @test dd.core_profiles === IMAS.top(dd.core_profiles.profiles_1d[1].grid.rho_tor_norm) # this does not work yet

    # test top() working on a sub structure
    @test crp1d === IMAS.top(crp1d; IDS_is_absolute_top = false)
    @test crp1d === IMAS.top(crp1d.grid; IDS_is_absolute_top = false)
    # @test crp1d === IMAS.top(crp1d.grid.rho_tor_norm; IDS_is_absolute_top=false) # this does not work yet
    @test crp1d === IMAS.top(crp1d)
    @test crp1d === IMAS.top(crp1d.grid)
    # @test crp1d === IMAS.top(crp1d.grid.rho_tor_norm) # this does not work yet

    # add structure to an array of structures
    push!(dd.core_profiles.profiles_1d, crp1d)
    @test dd.core_profiles === IMAS.top(crp1d)

    # test fail of adding dd without coordinate in IDS
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d, 1)
    @test_throws Exception dd.core_profiles.profiles_1d[1].electrons.temperature = Vector{Float64}(collect(1:10))

    # make sure an error is raised when trying to access missing dd
    @test_throws Exception dd.core_profiles.profiles_1d[1].j_total

    # resize! using function as condition to create new entry
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d)
    ion = resize!(dd.core_profiles.profiles_1d[].ion, ion -> (ion.z_ion == 1))
    @test ion === dd.core_profiles.profiles_1d[].ion[end]
    ion.z_ion = 1
    @test length(dd.core_profiles.profiles_1d[].ion) == 1
    ion = resize!(dd.core_profiles.profiles_1d[].ion, ion -> (ion.z_ion == 1))
    @test ion === dd.core_profiles.profiles_1d[].ion[end]
    @test length(dd.core_profiles.profiles_1d[].ion) == 1
    ion = resize!(dd.core_profiles.profiles_1d[].ion, ion -> (ion.z_ion == 2))
    @test ion === dd.core_profiles.profiles_1d[].ion[end]
    ion.z_ion = 2
    @test length(dd.core_profiles.profiles_1d[].ion) == 2
    ion = resize!(dd.core_profiles.profiles_1d[].ion, ion -> (ion.z_ion == 2))
    @test ion === dd.core_profiles.profiles_1d[].ion[end]
    @test length(dd.core_profiles.profiles_1d[].ion) == 2
    ion = resize!(dd.core_profiles.profiles_1d[].ion, 3)
    ion.z_ion = 2
    @test_throws Exception resize!(dd.core_profiles.profiles_1d[].ion, ion -> (ion.z_ion == 2))

    # resize! using pair of string and values as conditions to create new entry
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d)
    ion = resize!(dd.core_profiles.profiles_1d[].ion, "z_ion" => 1, "state[1].label" => "hello")
    @test ion === dd.core_profiles.profiles_1d[].ion[end]
    @test length(dd.core_profiles.profiles_1d[].ion) == 1
    ion = resize!(dd.core_profiles.profiles_1d[].ion, "z_ion" => 1, "state[1].label" => "hello")
    @test length(dd.core_profiles.profiles_1d[].ion) == 1
    ion = resize!(dd.core_profiles.profiles_1d[].ion, "z_ion" => 1)
    @test ion === dd.core_profiles.profiles_1d[].ion[end]
    @test length(dd.core_profiles.profiles_1d[].ion) == 1
    ion = resize!(dd.core_profiles.profiles_1d[].ion, "z_ion" => 2)
    @test ion === dd.core_profiles.profiles_1d[].ion[end]
    @test length(dd.core_profiles.profiles_1d[].ion) == 2
    ion = resize!(dd.core_profiles.profiles_1d[].ion, "z_ion" => 2, "state[1].label" => "hello")
    @test ion === dd.core_profiles.profiles_1d[].ion[end]
    @test length(dd.core_profiles.profiles_1d[].ion) == 3
    @test_throws Exception resize!(dd.core_profiles.profiles_1d[].ion, "z_ion" => 2)

end

@testset "IDS_IMAS" begin
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d, 2)

    # test f2u
    @test IMAS.f2u(dd.core_profiles.profiles_1d[1].grid) == "core_profiles.profiles_1d[:].grid"
    @test_throws MethodError IMAS.f2u(:core_profiles__profiles_1d___grid)
    @test_throws MethodError IMAS.f2u("core_profiles__profiles_1d___grid")

    # test i2p
    @test IMAS.i2p("core_profiles.profiles_1d[1].grid") == ["core_profiles", "profiles_1d", 1, "grid"]
    @test IMAS.i2p("core_profiles.profiles_1d[:].grid") == ["core_profiles", "profiles_1d", ":", "grid"]

    # test p2i
    @test IMAS.p2i(["core_profiles", "profiles_1d", 1, "grid"]) == "core_profiles.profiles_1d[1].grid"
    @test IMAS.p2i(["core_profiles", "profiles_1d", ":", "grid"]) == "core_profiles.profiles_1d[:].grid"
    @test_throws Exception IMAS.p2i([:core_profiles, :profiles_1d, ":", :grid])

    # test nested resizing
    wall = IMAS.wall()
    resize!(wall.description_2d, 1)
    resize!(wall.description_2d[1].mobile.unit, 2)
    resize!(wall.description_2d[1].mobile.unit[2].outline, 2)
    wall__description_2d = IMAS.wall__description_2d()
    resize!(wall__description_2d.mobile.unit, 2)
    resize!(wall__description_2d.mobile.unit[2].outline, 2)

    # test f2p
    @test IMAS.f2p(wall.description_2d[1].mobile.unit[2].outline[1]) == ["wall", "description_2d", 1, "mobile", "unit", 2, "outline", 1]
    @test IMAS.f2p(wall__description_2d.mobile.unit[2].outline[1]) == ["wall", "description_2d", 0, "mobile", "unit", 2, "outline", 1]

    # test imas_info
    @test IMAS.imas_info("core_profiles.profiles_1d[1]") == IMAS.imas_info("core_profiles.profiles_1d[:]")
    @test IMAS.imas_info("core_profiles.profiles_1d") == IMAS.imas_info("core_profiles.profiles_1d[:]")
    @test all([haskey(IMAS.imas_info("core_profiles.profiles_1d"), k) for k in ["coordinates", "data_type", "full_path", "documentation"]])
    @test_throws Exception IMAS.imas_info("core_profiles.does_not_exist")

    # test coordinate of a coordinate
    coords = IMAS.coordinates(dd.core_profiles.profiles_1d[1].grid, :rho_tor_norm)
    @test coords[:names][1] == "1...N"
    @test coords[:values][1] === nothing

    # test coordinate of a 1D array (with uninitialized coordinate)
    coords = IMAS.coordinates(dd.core_profiles.profiles_1d[1].electrons, :temperature)
    @test coords[:names][1] == "core_profiles.profiles_1d[:].grid.rho_tor_norm"
    @test coords[:values][1] === missing

    # test coordinate of a 1D array (with initialized coordinate)
    dd.core_profiles.profiles_1d[1].grid.rho_tor_norm = range(0, 1, length = 10)
    dd.core_profiles.profiles_1d[2].grid.rho_tor_norm = range(0, 1, length = 3)
    coords = IMAS.coordinates(dd.core_profiles.profiles_1d[1].electrons, :temperature)
    @test coords[:names][1] == "core_profiles.profiles_1d[:].grid.rho_tor_norm"
    @test coords[:values][1] === dd.core_profiles.profiles_1d[1].grid.rho_tor_norm
    @test length(coords[:values][1]) == 10
    coords = IMAS.coordinates(dd.core_profiles.profiles_1d[2].electrons, :temperature)
    @test coords[:names][1] == "core_profiles.profiles_1d[:].grid.rho_tor_norm"
    @test coords[:values][1] === dd.core_profiles.profiles_1d[2].grid.rho_tor_norm
    @test length(coords[:values][1]) == 3

    # # test change of grid on numerical array dd and function dd
    # dd = IMAS.dd();
    # resize!(dd.core_profiles.profiles_1d, 3)
    # dd.core_profiles.profiles_1d[1].grid.rho_tor_norm = range(0.0, stop=1.0, length=10)
    # dd.core_profiles.profiles_1d[2].grid.rho_tor_norm = range(0.0, stop=1.0, length=3)
    # dd.core_profiles.profiles_1d[1].electrons.temperature = 1.0 .- (dd.core_profiles.profiles_1d[1].grid.rho_tor_norm).^2
    # dd.core_profiles.profiles_1d[2].electrons.temperature = x -> 1.0 .- x.^2
    # dd.core_profiles.profiles_1d[1].grid.rho_tor_norm = range(0.0, stop=1.0, length=3)
    # @test length(dd.core_profiles.profiles_1d[1].electrons.temperature) == 3
    # dd.core_profiles.profiles_1d[2].grid.rho_tor_norm = range(0.0, stop=1.0, length=4)
    # @test length(dd.core_profiles.profiles_1d[2].electrons.temperature) == 4

    # test working with IDSvectorElement standalone or in a IDSvector
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d, 1)
    for profiles_1d in [dd.core_profiles.profiles_1d[1], IMAS.core_profiles__profiles_1d()]
        profiles_1d.grid.rho_tor_norm = range(0.0, 1.0, length = 101)
        profiles_1d.electrons.density = (x; _...) -> (1.0 .- x .^ 2) .^ 2.0
        profiles_1d.j_total = (x; _...) -> (1.0 .- x .^ 2) .^ 2.0
        @test length(profiles_1d.electrons.density) == length(profiles_1d.j_total)
    end
end


@testset "JSON_IO" begin
    filename = joinpath(dirname(dirname(abspath(@__FILE__))), "sample", "sample_equilibrium_ods.json")
    dd = IMAS.json2imas(filename; verbose = false)
    @test length(dd.wall.description_2d[1].limiter.unit[1].outline.r) > 0
end

using Revise
using IMAS
using Test

@testset "time_ids" begin
    dd = IMAS.dd()

    dd.global_time = 1010.0
    eqt = resize!(dd.equilibrium.time_slice)
    @test length(dd.equilibrium.time_slice) == 1
    n = length(dd.equilibrium.time_slice)
    eqt.global_quantities.ip = 1.0
    @test dd.equilibrium.time_slice[].global_quantities.ip == 1.0
    @test dd.equilibrium.time_slice[1010.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[n].global_quantities.ip === eqt.global_quantities.ip
    @test_throws Exception dd.equilibrium.time_slice[n+1].global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip
    @test (dd.equilibrium.time_slice[] = eqt) == eqt

    dd.global_time = 2020.0
    eqt = resize!(dd.equilibrium.time_slice)
    @test length(dd.equilibrium.time_slice) == 2
    n = length(dd.equilibrium.time_slice)
    dd.equilibrium.time_slice[].global_quantities.ip = 2.0
    @test dd.equilibrium.time_slice[].global_quantities.ip == 2.0
    @test dd.equilibrium.time_slice[2020.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[n].global_quantities.ip === eqt.global_quantities.ip
    @test_throws Exception dd.equilibrium.time_slice[n+1].global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip
    @test (dd.equilibrium.time_slice[] = eqt) == eqt

    dd.global_time = 3030.0
    eqt = resize!(dd.equilibrium.time_slice, IMAS.τ)
    @test length(dd.equilibrium.time_slice) == 3
    n = length(dd.equilibrium.time_slice)
    dd.equilibrium.time_slice[].global_quantities.ip = 3.0
    @test dd.equilibrium.time_slice[].global_quantities.ip == 3.0
    @test dd.equilibrium.time_slice[3030.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[3000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[n].global_quantities.ip === eqt.global_quantities.ip
    @test_throws Exception dd.equilibrium.time_slice[n+1].global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip
    @test (dd.equilibrium.time_slice[] = eqt) == eqt

    # dial back global time at an existing time (closest time slice is 2)
    dd.global_time = 2020.0
    eqt = dd.equilibrium.time_slice[2]
    @test length(dd.equilibrium.time_slice) == 3
    @test dd.equilibrium.time_slice[2020.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip
    @test (dd.equilibrium.time_slice[] = eqt) == eqt

    # dial back global time at an existing time (closest time slice is 2)
    dd.global_time = 2013.0
    eqt = dd.equilibrium.time_slice[2]
    @test length(dd.equilibrium.time_slice) == 3
    @test dd.equilibrium.time_slice[2020.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip
    @test_throws Exception dd.equilibrium.time_slice[] = eqt # setindex! will complain trying to enter data at an earlier time that is not in time array

    # resize! will complain if trying to resize at an earlier time
    @test_throws Exception resize!(dd.equilibrium.time_slice)

    # resize! with global time will complain operating on IDSvectors that are not time dependent
    @test_throws Exception resize!(dd.wall.description_2d)
    @test_throws Exception resize!(dd.wall.description_2d, IMAS.τ)
    @test_throws Exception resize!(dd.wall.description_2d, 1000.0)

    # add an empty time-slice
    dd.global_time = 4040.0
    eqt = resize!(dd.equilibrium.time_slice, IMAS.τ)
    @test length(dd.equilibrium.time_slice) == 4

    dd.global_time = 5050.0
    eqt = resize!(dd.equilibrium.time_slice, IMAS.τ)
    @test length(dd.equilibrium.time_slice) == 5
    @test_throws Exception dd.equilibrium.time_slice[].global_quantities.ip # elements within arrays of structures do not time-interpolate
    dd.equilibrium.time_slice[].global_quantities.ip = 5.0
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip

    # https://github.com/ProjectTorreyPines/FUSE.jl/issues/18
    dd = IMAS.dd()
    for i in 1:2
        isource = resize!(dd.core_sources.source, "identifier.index" => 3)
        resize!(isource.profiles_1d)
        @test dd.core_sources.source[1] === dd.core_sources.source[end] === isource
    end
    @test length(dd.core_sources.time) == 1

    # resize!(dd.equilibrium.time_slice) returns an empty IDS
    dd = IMAS.dd()
    eqt = resize!(dd.equilibrium.time_slice)
    eqt.global_quantities.ip = 1.0
    eqt = resize!(dd.equilibrium.time_slice)
    @test ismissing(eqt.global_quantities, :ip)
    @test !ismissing(eqt, :time)

end

@testset "time_array" begin
    dd = IMAS.dd()

    dd.global_time = 1010.0
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = 1.0) == 1.0
    @test dd.equilibrium.time == [1010.0]
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0) == 1.0
    @test dd.equilibrium.vacuum_toroidal_field.b0 == [1.0]

    dd.global_time = 2020.0
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = 2.0) == 2.0
    @test dd.equilibrium.time == [1010.0, 2020.0]
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0) == 2.0
    @test dd.equilibrium.vacuum_toroidal_field.b0 == [1.0, 2.0]

    dd.global_time = 3030.0
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0) == 2.0 # time interpolation
    @test dd.equilibrium.time == [1010.0, 2020.0]
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = 3.0) == 3.0
    @test dd.equilibrium.time == [1010.0, 2020.0, 3030.0]
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0) == 3.0
    @test dd.equilibrium.vacuum_toroidal_field.b0 == [1.0, 2.0, 3.0]

    # edit something in the past
    dd.global_time = 2020.0
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = -2.0) == -2.0

    # insert a time
    push!(dd.equilibrium.time, 4040.0)

    # test interpolation and insertion
    dd.global_time = 5050.0
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = 5.0) == 5.0
    @test dd.equilibrium.vacuum_toroidal_field.b0 == [1.0, -2.0, 3.0, 3.0, 5.0]

    # test insertion at later time of an empty field
    empty!(dd.equilibrium.vacuum_toroidal_field.b0)
    @test @ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = 5.0) == 5.0
    @test all(dd.equilibrium.vacuum_toroidal_field.b0 .=== [NaN, NaN, NaN, NaN, 5.0])

    # test write/read at times specified as arguments
    time = 2020.0
    @test IMAS.set_time_array(dd.equilibrium.vacuum_toroidal_field, :b0, time, -2.0) == -2.0
    @test IMAS.get_time_array(dd.equilibrium.vacuum_toroidal_field, :b0, time) == -2.0

end

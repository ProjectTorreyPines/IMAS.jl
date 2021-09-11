@testset "parameters_set_get" begin
    # begin adding parameters to simple dictionary
    run_parameters = Dict()

    # scalar parameter
    run_parameters[:testScalar] = IMAS.ScalarParameter(0.5, "m/s", "scalar test")
    
    # switch parameter
    options = Dict()
    options[:bla] = IMAS.SwitchOption(:bla, "m/s", "bla")
    options[:ta] = IMAS.SwitchOption(:ta, "1/s", "ta")
    options[:user] = IMAS.ScalarParameter(1.2, "1/s", "user defined")
    run_parameters[:testSwitch] = IMAS.SwitchParameter(options, :ta, "switch test")

    # test that defining a default that is not an option throws an error 
    @test_throws Exception IMAS.SwitchParameter(options, :does_not_exist, "failing switch test")
    
    # define ImasParameters for rapid access of .value property
    local_parameters = IMAS.ImasParameters(run_parameters)
    
    @test local_parameters[:testScalar] == 0.5
    # test mutability
    local_parameters[:testScalar] = 1.0
    @test local_parameters[:testScalar] == 1.0

    @test local_parameters[:testSwitch] == :ta
    # test mutability
    local_parameters[:testSwitch] = :bla
    @test local_parameters[:testSwitch] == :bla
    # test setting something that is not a valid option
    @test_throws Exception local_parameters[:testSwitch] = :does_not_exist

    # test user defined value
    local_parameters[:testSwitch] = :user
    @test local_parameters[:testSwitch] == 1.2
    local_parameters[:testSwitch] = :user => 1.0
    @test local_parameters[:testSwitch] == 1.0
    local_parameters[:testSwitch] = :user => 2.0
    @test local_parameters[:testSwitch] == 2.0
    @test_throws Exception local_parameters[:testSwitch] = 3.0
    @test_throws Exception local_parameters[:testSwitch] = :does_not_exist => 1.0

    # test access to some parameter attribute
    @test local_parameters.parameters[:testSwitch].description == "switch test"
end

@testset "parameters_access" begin
    # test the imas_parameters
    @test typeof(IMAS.imas_parameters[:PHYSICS_MODELS][:bootstrapModel]) <: IMAS.CBSGi
    IMAS.imas_parameters[:PHYSICS_MODELS][:bootstrapModel] = :user => 0.5
    @test IMAS.imas_parameters[:PHYSICS_MODELS][:bootstrapModel] == 0.5
end
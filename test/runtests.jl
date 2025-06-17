# Check if specific test files are requested via ARGS
if !isempty(ARGS)
    for testfile in ARGS
        @info "Running test file: $testfile"
        include(testfile)
    end
else
    # Default behavior: run all tests
    include("runtests_fluxsurfaces.jl")

    include("runtests_extract.jl")

    include("runtests_plot_recipes.jl")

    include("physics/particles/test_toroidal_intersection.jl")
end

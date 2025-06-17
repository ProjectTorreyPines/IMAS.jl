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
    
    include("geometry/test_ray_torus_intersect.jl")
    
    include("geometry/test_solve_r_intersect.jl")
end

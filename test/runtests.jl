# Check if specific test files are requested via ARGS
if !isempty(ARGS)
    for testfile in ARGS
        @info "Running test file: $testfile"
        include(testfile)
    end
else
    # Default behavior: run all tests
    include("runtests_refine_extremum.jl")

    include("runtests_fluxsurfaces_cubic.jl")

    include("runtests_fluxsurfaces.jl")

    include("runtests_interpolations.jl")

    include("runtests_interp1d.jl")

    include("runtests_fields.jl")

    include("runtests_fast.jl")

    include("runtests_extract.jl")

    include("runtests_plot_recipes.jl")
end

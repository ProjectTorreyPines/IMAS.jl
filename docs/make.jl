using Documenter, IMAS

# Call functions
open(joinpath(@__DIR__, "src/api.md"), "w") do f
    println(f, "# API Reference\n")
    for page in sort!(collect(keys(IMAS.document)))
        if page == :Expressions
            continue
        end
        println(f, "## $page\n")
        println(f, "```@docs")
        for item in IMAS.document[page]
            println(f, "IMAS.$item")
        end
        println(f, "```")
    end
end

open(joinpath(@__DIR__, "src/expressions.md"), "w") do f
    println(f, "# Expressions\n")
    println(f, "```@docs")
    for item in IMAS.document[:Expressions]
        println(f, "IMAS.$item")
    end
    return println(f, "```")
end

makedocs(;
    modules=[IMAS],
    format=Documenter.HTML(;
        analytics="G-65D8V8C8VQ",
        size_threshold=nothing,
        size_threshold_warn=nothing),
    sitename="IMAS",
    checkdocs=:none,
    pages=[
        "index.md",
        "api.md",
        "expressions.md",
        "License" => "license.md", "Notice" => "notice.md"],
    warnonly=true
)

# Deploy docs
# This function deploys the documentation to the gh-pages branch of the repository.
# The main documentation that will be hosted on
# https://projecttorreypines.github.io/IMAS.jl/stable
# will be built from latest release tagged with a version number.
# The development documentation that will be hosted on
# https://projecttorreypines.github.io/IMAS.jl/dev
# will be built from the latest commit on the chosen devbranch argument below.
# For testing purposes, the devbranch argument can be set to WIP branch like "docs".
# While merging with master, the devbranch argument should be set to "master".
deploydocs(;
    repo="github.com/ProjectTorreyPines/IMAS.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "v#.#"]
)

__precompile__()

module FUSE

using StructArrays

#= ==================================== =#
# FUSE data structure
#= ==================================== =#


Base.@kwdef mutable struct core_profiles_profiles_1d_grid
    rho_tor_norm :: Vector{Float64} = []
end

Base.@kwdef mutable struct core_profiles_profiles_1d_electrons
    density :: Vector{Float64} = []
    density_fast :: Vector{Float64} = []
end

Base.@kwdef struct core_profiles_profiles_1d
    electrons :: core_profiles_profiles_1d_electrons = core_profiles_profiles_1d_electrons()
    grid :: core_profiles_profiles_1d_grid = core_profiles_profiles_1d_grid()
end

Base.@kwdef struct core_profiles
    profiles_1d :: StructArray{core_profiles_profiles_1d} = StructArray(core_profiles_profiles_1d() for k in 1:1)
end

Base.@kwdef struct dd
    core_profiles :: core_profiles = core_profiles()
end

end # module

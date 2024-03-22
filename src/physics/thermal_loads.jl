mutable struct WallHeatFlux
    r::Vector{Float64}                  # R coordinate of the wall mesh                     - m
    z::Vector{Float64}                  # Z coordinate of the wall mesh                     - m
    q_wall::Vector{Float64}             # total heat flux on the wall                       - W/m2
    q_part::Vector{Float64}             # heat flux on the wall due to particles            - W/m2
    q_core_rad::Vector{Float64}         # heat flux on the wall due to core radiation       - W/m2
    q_parallel::Vector{Float64}         # parallel heat flux due to particles at the wall   - W/m2
    s::Vector{Float64}                  # wall curvilinear abscissa                         - m
end

"""
WallHeatFlux(; 
    r::Vector{Float64} = Float64[],
    z::Vector{Float64} = Float64[],
    q_wall::Vector{Float64} = Float64[],
    q_part::Vector{Float64} = Float64[],
    q_core_rad::Vector{Float64} = Float64[],
    q_parallel::Vector{Float64} = Float64[],
    s::Vector{Float64} = Float64[])

    Initializes a WallHeatFlux struct. IMAS.WallHeatFlux() returns a WallHeatFlux with all empty entries.

"""
function WallHeatFlux(; 
    r::Vector{Float64} = Float64[],
    z::Vector{Float64} = Float64[],
    q_wall::Vector{Float64} = Float64[],
    q_part::Vector{Float64} = Float64[],
    q_core_rad::Vector{Float64} = Float64[],
    q_parallel::Vector{Float64} = Float64[],
    s::Vector{Float64} = Float64[])
    return WallHeatFlux(r,z,q_wall,q_part,q_core_rad,q_parallel,s)
end

document[Symbol("Physics diagnostics")] = Symbol[]

"""
    lenght_line_of_sight(icls::IMAS.interferometer__channel___line_of_sight)

Returns lenght of a given interferometer channel line of sight
"""
function lenght_line_of_sight(icls::IMAS.interferometer__channel___line_of_sight)
    x1 = icls.first_point.r * cos(icls.first_point.phi)
    y1 = icls.first_point.r * sin(icls.first_point.phi)
    z1 = icls.first_point.z

    x2 = icls.second_point.r * cos(icls.second_point.phi)
    y2 = icls.second_point.r * sin(icls.second_point.phi)
    z2 = icls.second_point.z

    d = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)

    if !isempty(icls.third_point)
        x3 = icls.third_point.r * cos(icls.third_point.phi)
        y3 = icls.third_point.r * sin(icls.third_point.phi)
        z3 = icls.third_point.z

        d += sqrt((x3 - x2)^2 + (y3 - y2)^2 + (z3 - z2)^2)
    end

    return d
end

"""
    line_average(
        rho_interp, 
        lcfs_r::AbstractVector{<:Real}, 
        lcfs_z::AbstractVector{<:Real},
        q::AbstractVector{<:Real}, 
        rho_tor_norm::AbstractVector{<:Real},
        r1::T1, z1::T1, r2::T1, z2::T1;
        n_points::Int=100
    ) where {T1<:Real}

Low-level function that takes pre-computed ρ interpolant and LCFS boundary arrays.
This is more efficient when computing line averages for multiple channels or time-slices
with the same equilibrium.

Returns named tuple with (:line_integral, :line_average, :path_length_inside_lcfs)
"""
function line_average(
    rho_interp,
    lcfs_r::AbstractVector{<:Real}, 
    lcfs_z::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    rho_tor_norm::AbstractVector{<:Real},
    r1::T1, z1::T1, r2::T1, z2::T1;
    n_points::Int=100
) where {T1<:Real}

    T = promote_type(eltype(q), eltype(rho_tor_norm), typeof(r1), typeof(z1), typeof(r2), typeof(z2))
    
    # Create 1D interpolant for the quantity q
    q_interp = cubic_interp1d(rho_tor_norm, q)
    
    # Find intersections between line and LCFS
    intersections = intersection([r1, r2], [z1, z2], lcfs_r, lcfs_z)
    
    # Initialize integration variables
    line_integral = T(0)
    path_length_inside_lcfs = T(0)
    
    if length(intersections.crossings) >= 2
        # Sort intersection points along the line parameter
        intersection_params = T[]
        for (crossing_r, crossing_z) in intersections.crossings
            # Calculate parameter t where line intersects LCFS
            if abs(r2 - r1) > abs(z2 - z1)
                t = (crossing_r - r1) / (r2 - r1)
            else
                t = (crossing_z - z1) / (z2 - z1)
            end
            push!(intersection_params, clamp(t, T(0), T(1)))
        end
        sort!(intersection_params)
        
        # Take first and last intersection (entry and exit points)
        t_entry = intersection_params[1]
        t_exit = intersection_params[end]
        
        # Calculate entry and exit points
        r_entry = r1 + t_entry * (r2 - r1)
        z_entry = z1 + t_entry * (z2 - z1)
        r_exit = r1 + t_exit * (r2 - r1)
        z_exit = z1 + t_exit * (z2 - z1)
        
        # Create integration points along the segment inside LCFS
        s_params = range(T(0), T(1), n_points)
        r_seg = r_entry .+ (r_exit - r_entry) .* s_params
        z_seg = z_entry .+ (z_exit - z_entry) .* s_params
        
        # Calculate ρ_tor_norm along the segment
        rho_seg = rho_interp.RHO_interpolant.(r_seg, z_seg)
        
        # Interpolate quantity q along the segment
        q_seg = q_interp.(rho_seg)
        
        # Calculate path length inside LCFS
        path_length_inside_lcfs = sqrt((r_exit - r_entry)^2 + (z_exit - z_entry)^2)
        
        # Create coordinate array for integration
        s_coord = s_params .* path_length_inside_lcfs
        line_integral = trapz(s_coord, q_seg)
    end
    
    # Calculate line average
    line_average = path_length_inside_lcfs > 0 ? line_integral / path_length_inside_lcfs : T(0)
    
    return (line_integral=line_integral, line_average=line_average, path_length_inside_lcfs=path_length_inside_lcfs)
end

"""
    _get_rho_interp_and_lcfs(eqt::IMAS.equilibrium__time_slice)

Memoized helper function to extract and cache ρ interpolant and LCFS boundary.
This automatically handles caching for repeated calls with the same equilibrium.
"""
Memoize.@memoize function _get_rho_interp_and_lcfs(eqt::IMAS.equilibrium__time_slice)
    rho_interp = ρ_interpolant(eqt)
    lcfs_r = eqt.boundary.outline.r
    lcfs_z = eqt.boundary.outline.z
    return (rho_interp=rho_interp, lcfs_r=lcfs_r, lcfs_z=lcfs_z)
end

"""
    line_average(
        eqt::IMAS.equilibrium__time_slice, 
        q::AbstractVector{<:Real}, 
        rho_tor_norm::AbstractVector{<:Real},
        r1::Real, z1::Real, r2::Real, z2::Real;
        n_points::Int=100
    )

High-level convenience function with automatic memoization of expensive ρ interpolant creation.
Automatically caches and reuses interpolants for the same equilibrium across multiple calls.

Returns named tuple with (:line_integral, :line_average, :path_length_inside_lcfs)
"""
function line_average(
    eqt::IMAS.equilibrium__time_slice,
    q::AbstractVector{<:Real},
    rho_tor_norm::AbstractVector{<:Real},
    r1::T1, z1::T1, r2::T1, z2::T1;
    n_points::Int=100
) where {T1<:Real}
    
    # Get cached ρ interpolant and LCFS boundary
    rho_interp, lcfs_r, lcfs_z = _get_rho_interp_and_lcfs(eqt)
    
    # Call lower-level function
    return line_average(rho_interp, lcfs_r, lcfs_z, q, rho_tor_norm, r1, z1, r2, z2; n_points)
end

"""
    ne_line(
        equilibrium::IMAS.equilibrium,
        core_profiles::IMAS.core_profiles,
        interferometer::IMAS.interferometer;
        times::Vector{Float64}=equilibrium.time
    )

Calculate time-dependent line average electron density for all interferometer channels.

# Arguments
- `equilibrium`: equilibrium IDS containing time-dependent equilibrium data
- `core_profiles`: core_profiles IDS containing electron density profiles
- `interferometer`: interferometer IDS containing channel line-of-sight information
- `times`: time points to evaluate (defaults to equilibrium.time)

# Returns
Vector of vectors, where `result[i]` contains the time-dependent line average density 
for channel `i` at the specified time points.

# Example
```julia
# Calculate for all channels at equilibrium time points
ne_line_avg = ne_line(dd.equilibrium, dd.core_profiles, dd.interferometer)

# Calculate at specific times
ne_line_avg = ne_line(dd.equilibrium, dd.core_profiles, dd.interferometer; 
                     times=[1.0, 2.0, 3.0])
```
"""
function ne_line(
    equilibrium::IMAS.equilibrium,
    core_profiles::IMAS.core_profiles,
    interferometer::IMAS.interferometer;
    times::Vector{Float64}=equilibrium.time,
    n_points::Int=100
)
    n_channels = length(interferometer.channel)
    n_times = length(times)
    
    # Initialize result array - each channel gets a time series
    result = Vector{Vector{Float64}}(undef, n_channels)
    
    for (ch_idx, channel) in enumerate(interferometer.channel)
        # Initialize time series for this channel
        channel_ne_line = Vector{Float64}(undef, n_times)
        
        # Get line-of-sight coordinates for this channel
        r1 = channel.line_of_sight.first_point.r
        z1 = channel.line_of_sight.first_point.z
        r2 = channel.line_of_sight.second_point.r
        z2 = channel.line_of_sight.second_point.z
        
        for (t_idx, time0) in enumerate(times)
            try
                # Get equilibrium and core profiles at this time
                eqt = equilibrium.time_slice[time0]
                cp1d = core_profiles.profiles_1d[time0]
                
                # Calculate line average density
                result_tmp = line_average(
                    eqt, 
                    cp1d.electrons.density, 
                    cp1d.grid.rho_tor_norm,
                    r1, z1, r2, z2;
                    n_points
                )
                
                channel_ne_line[t_idx] = result_tmp.line_average
                
            catch e
                # Handle missing data gracefully
                channel_ne_line[t_idx] = NaN
            end
        end
        
        result[ch_idx] = channel_ne_line
    end
    
    return (time=times, channel = result)
end

@compat public line_average, ne_line
push!(document[Symbol("Physics diagnostics")], :line_average, :ne_line)

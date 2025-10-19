document[Symbol("Physics interferometer")] = Symbol[]

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
        r1::T1, z1::T1, ϕ1::T1,
        r2::T1, z2::T1, ϕ2::T1;
        n_points::Int=100
    ) where {T1<:Real}

Low-level function that takes pre-computed ρ interpolant and LCFS boundary arrays.

This is more efficient when computing line averages for multiple channels or time-slices with the same equilibrium.

Returns named tuple with (:line_integral, :line_average, :path_length)
"""
function line_average(
    rho_interp,
    lcfs_r::AbstractVector{<:Real},
    lcfs_z::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    rho_tor_norm::AbstractVector{<:Real},
    r1::T1, z1::T1, ϕ1::T1,
    r2::T1, z2::T1, ϕ2::T1;
    n_points::Int=100
) where {T1<:Real}

    T = promote_type(eltype(q), eltype(rho_tor_norm), typeof(r1), typeof(z1), typeof(r2), typeof(z2))

    @assert ϕ1 == ϕ2 "line_average cannot handle different ϕ yet"

    # Create 1D interpolant for the quantity q
    q_interp = cubic_interp1d(rho_tor_norm, q)

    # Find intersections between line and LCFS
    intersections = intersection([r1, r2], [z1, z2], lcfs_r, lcfs_z)
    @assert length(intersections.crossings) in (0, 2)

    # Calculate path length
    path_length = sqrt((r1 - r2)^2 + (z1 - z2)^2)
    path_length_inside_lcfs = T(0)

    line_integral = T(0)
    if length(intersections.crossings) == 2
        (r_entry, z_entry), (r_exit, z_exit) = intersections.crossings

        # Create integration points along the segment inside LCFS
        r_seg = range(r_entry, r_exit, n_points)
        z_seg = range(z_entry, z_exit, n_points)

        # Calculate ρ_tor_norm along the segment
        rho_seg = rho_interp.RHO_interpolant.(r_seg, z_seg)

        # Interpolate quantity q along the segment
        q_seg = q_interp.(rho_seg)

        # Calculate path length inside LCFS
        path_length_inside_lcfs = sqrt((r_exit - r_entry)^2 + (z_exit - z_entry)^2)

        # Create coordinate array for integration
        line_integral = trapz(range(0, path_length_inside_lcfs, n_points), q_seg)
    end

    return (line_integral=line_integral, line_average=line_integral / path_length, path_length=path_length)
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
        line_of_sight;
        n_points::Int=100
    )

High-level convenience function for computing line averages

Returns named tuple with (:line_integral, :line_average, :path_length)
"""
function line_average(
    eqt::IMAS.equilibrium__time_slice,
    q::AbstractVector{<:Real},
    rho_tor_norm::AbstractVector{<:Real},
    line_of_sight;
    n_points::Int=100
)
    if !isempty(line_of_sight.third_point)
        return line_average(
            eqt,
            q,
            rho_tor_norm,
            line_of_sight.first_point.r, line_of_sight.first_point.z, line_of_sight.first_point.phi,
            line_of_sight.second_point.r, line_of_sight.second_point.z, line_of_sight.second_point.phi,
            line_of_sight.third_point.r, line_of_sight.third_point.z, line_of_sight.third_point.phi;
            n_points
        )
    else
        return line_average(
            eqt,
            q,
            rho_tor_norm,
            line_of_sight.first_point.r, line_of_sight.first_point.z, line_of_sight.first_point.phi,
            line_of_sight.second_point.r, line_of_sight.second_point.z, line_of_sight.second_point.phi;
            n_points
        )
    end
end

function line_average(
    eqt::IMAS.equilibrium__time_slice,
    q::AbstractVector{<:Real},
    rho_tor_norm::AbstractVector{<:Real},
    r1::T1, z1::T1, ϕ1::T1,
    r2::T1, z2::T1, ϕ2::T1;
    n_points::Int
) where {T1<:Real}

    # Get cached ρ interpolant and LCFS boundary
    rho_interp, lcfs_r, lcfs_z = _get_rho_interp_and_lcfs(eqt)

    # Call lower-level function
    return line_average(rho_interp, lcfs_r, lcfs_z, q, rho_tor_norm, r1, z1, ϕ1, r2, z2, ϕ2; n_points)
end

function line_average(
    eqt::IMAS.equilibrium__time_slice,
    q::AbstractVector{<:Real},
    rho_tor_norm::AbstractVector{<:Real},
    r1::T1, z1::T1, ϕ1::T1,
    r2::T1, z2::T1, ϕ2::T1,
    r3::T1, z3::T1, ϕ3::T1;
    n_points::Int
) where {T1<:Real}

    # Get cached ρ interpolant and LCFS boundary
    rho_interp, lcfs_r, lcfs_z = _get_rho_interp_and_lcfs(eqt)

    # Call lower-level function
    tmp12 = line_average(rho_interp, lcfs_r, lcfs_z, q, rho_tor_norm, r1, z1, ϕ1, r2, z2, ϕ2; n_points)
    tmp23 = line_average(rho_interp, lcfs_r, lcfs_z, q, rho_tor_norm, r2, z2, ϕ2, r3, z3, ϕ3; n_points)

    return (
        line_integral=tmp12.line_integral + tmp23.line_integral,
        line_average=tmp12.line_average,
        path_length=tmp12.path_length + tmp23.path_length)
end

@compat public line_average
push!(document[Symbol("Physics interferometer")], :line_average)

"""
    ne_line(
        equilibrium::IMAS.equilibrium,
        core_profiles::IMAS.core_profiles,
        interferometer::IMAS.interferometer;
        times::Vector{Float64}=equilibrium.time
    )

Calculate time-dependent line average electron density for all interferometer channels.

Returns new interferometer IDS with synthetic data.
"""
function ne_line(
    equilibrium::IMAS.equilibrium,
    core_profiles::IMAS.core_profiles,
    interferometer::IMAS.interferometer{T};
    times::Vector{Float64}=equilibrium.time,
    n_points::Int=100
) where {T<:Real}

    n_channels = length(interferometer.channel)
    n_times = length(times)

    # Initialize result array - each channel gets a time series
    interferometer_out = IMAS.interferometer{T}()
    resize!(interferometer_out.channel, n_channels)

    for (ch_idx, (ch, chout)) in enumerate(zip(interferometer.channel, interferometer_out.channel))
        # Initialize time series for this channel
        chout.n_e_line.time = times
        chout.n_e_line.data = Vector{T}(undef, n_times)
        chout.name = ch.name

        # fill line-of-sight for this channel
        fill!(chout.line_of_sight, ch.line_of_sight)

        for (t_idx, time0) in enumerate(times)
            try
                # Get equilibrium and core profiles at this time
                eqt = equilibrium.time_slice[time0]
                cp1d = core_profiles.profiles_1d[time0]

                # Calculate line integral of the density
                density_thermal = line_average(eqt, cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, ch.line_of_sight; n_points)
                chout.n_e_line.data[t_idx] = density_thermal.line_integral

            catch e
                if typeof(e) <: IMASdd.IMASbadTime
                    # Handle missing data gracefully
                    chout.n_e_line.data[t_idx] = NaN
                else
                    rethrow(e)
                end
            end
        end
    end

    return interferometer_out
end

@compat public ne_line
push!(document[Symbol("Physics interferometer")], :ne_line)

"""
    interferometer!(
        intf::IMAS.interferometer,
        eqt::IMAS.equilibrium__time_slice,
        cp1d::IMAS.core_profiles__profiles_1d;
        time0::Float64=eqt.time,
        n_points::Int=100
    )

Populate interferometer IDS with synthetic line-averaged and line-integrated electron density measurements for a single time slice.
"""
function interferometer!(
    intf::IMAS.interferometer,
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d;
    time0::Float64=eqt.time,
    n_points::Int=100
)
    for ch in intf.channel
        # Calculate line average density for this time slice
        density_result = line_average(
            eqt,
            cp1d.electrons.density_thermal,
            cp1d.grid.rho_tor_norm,
            ch.line_of_sight;
            n_points
        )

        # Store using set_time_array
        IMAS.set_time_array(ch.n_e_line_average, :data, time0, density_result.line_average)
        IMAS.set_time_array(ch.n_e_line, :data, time0, density_result.line_integral)
    end

    return intf
end

"""
    interferometer!(
        intf::IMAS.interferometer,
        eq::IMAS.equilibrium,
        cp::IMAS.core_profiles;
        n_points::Int=100
    )

Populate interferometer IDS with synthetic line-averaged and line-integrated electron density measurements for all time slices.
"""
function interferometer!(
    intf::IMAS.interferometer,
    eq::IMAS.equilibrium,
    cp::IMAS.core_profiles;
    n_points::Int=100
)

    # remove existing data
    for channel in intf.channel
        empty!(channel.n_e)
        empty!(channel.n_e_line)
        empty!(channel.n_e_line_average)
    end

    # calculate new data for all equilibrium time slices
    for time0 in eq.time
        eqt = eq.time_slice[time0]
        cp1d = cp.profiles_1d[time0]
        interferometer!(intf, eqt, cp1d; time0, n_points)
    end

    return intf
end

@compat public interferometer!
push!(document[Symbol("Physics interferometer")], :interferometer!)

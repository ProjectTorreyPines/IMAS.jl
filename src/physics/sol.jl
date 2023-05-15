struct OpenFieldLine
    r::Vector{Float64}
    z::Vector{Float64}
    Br::Vector{Float64}
    Bz::Vector{Float64}
    Bp::Vector{Float64}
    Bt::Vector{Float64}
    pitch::Vector{Float64}
    s::Vector{Float64}
    midplane_index::Int
    strike_angles::Vector{Float64} # poloidal strike angle measured from the normal
end

@recipe function plot_ofl(ofl::OpenFieldLine)
    @series begin
        aspect_ratio --> :equal
        label --> ""
        colorbar_title := "Connection length [m]"
        line_z := ofl.s
        ofl.r, ofl.z
    end
end

@recipe function plot_OFL(OFL_hfs_lfs::Tuple{Vector{OpenFieldLine},Vector{OpenFieldLine}})
    for OFL in OFL_hfs_lfs
        @series begin
            OFL
        end
    end
end

@recipe function plot_OFL(OFL::Vector{OpenFieldLine})
    for ofl in OFL
        @series begin
            ofl
        end
    end
end

function sol(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall; levels::Int=1)
    return sol(eqt, first_wall(wall).r, first_wall(wall).z; levels)
end

"""
    sol(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{T}, wall_z::Vector{T}, levels::Int=1) where {T<:Real}

Returns vectors of hfs and lfs OpenFieldLine
"""
function sol(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{T}, wall_z::Vector{T}; levels::Int=1) where {T<:Real}
    r0, b0 = vacuum_r0_b0(eqt)

    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    ############
    r, z, PSI_interpolant = ψ_interpolant(eqt)
    r_wall_midplane, _ = intersection([R0, maximum(wall_r)], [Z0, Z0], wall_r, wall_z; as_list_of_points=false)
    psi_wall_midplane = PSI_interpolant.(r_wall_midplane, Z0)[1]
    psi__axis_level = eqt.profiles_1d.psi[1]
    psi__boundary_level = find_psi_boundary(eqt; raise_error_on_not_open=true)
    psi_sign = sign(psi__boundary_level - psi__axis_level)
    ############

    # pack points near lcfs
    levels = psi__boundary_level .+ psi_sign .* 10.0 .^ LinRange(-3, log10(abs(psi_wall_midplane - psi__boundary_level)), levels + 1)[1:end-1]

    OFL_hfs = OpenFieldLine[]
    OFL_lfs = OpenFieldLine[]
    for level in levels
        lines = flux_surface(eqt, level, false)
        for (r, z) in lines
            rr, zz, strike_angles = line_wall_2_wall(r, z, wall_r, wall_z, R0, Z0)
            if isempty(rr) || all(zz .> Z0) || all(zz .< Z0)
                continue
            end

            # add a point exactly at the midplane
            crossing_index, r_midplane, z_midplane = intersection([minimum(wall_r), maximum(wall_r)], [Z0, Z0], rr, zz; as_list_of_points=false, return_indexes=true)
            rr = [rr[1:crossing_index[1][2]]; r_midplane; rr[crossing_index[1][2]+1:end]]
            zz = [zz[1:crossing_index[1][2]]; z_midplane; zz[crossing_index[1][2]+1:end]]
            midplane_index = crossing_index[1][2] + 1

            # calculate quantities along field line
            Br, Bz = Br_Bz_vector_interpolant(PSI_interpolant, rr, zz)
            Bp = sqrt.(Br .^ 2.0 .+ Bz .^ 2.0)
            Bt = abs.(b0 .* r0 ./ rr)
            dp = sqrt.(gradient(rr) .^ 2.0 .+ gradient(zz) .^ 2.0)
            pitch = sqrt.(1.0 .+ (Bt ./ Bp) .^ 2)
            s = cumsum(pitch .* dp)
            s = abs.(s .- s[midplane_index])

            # select HFS or LFS and add line to the list
            if rr[midplane_index] < R0
                OFL = OFL_hfs
            else
                OFL = OFL_lfs
            end
            push!(OFL, OpenFieldLine(rr, zz, Br, Bz, Bp, Bt, pitch, s, midplane_index, strike_angles))
        end
    end
    return OFL_hfs, OFL_lfs
end

"""
    line_wall_2_wall(r::T, z::T, wall_r::T, wall_z::T, R0::Real, Z0::Real) where {T<:AbstractVector{<:Real}}

Returns r, z coordinates of open field line contained within wall, as well as angles of incidence at the strike locations
"""
function line_wall_2_wall(r::T, z::T, wall_r::T, wall_z::T, R0::Real, Z0::Real) where {T<:AbstractVector{<:Real}}

    indexes, crossings = intersection(r, z, wall_r, wall_z; as_list_of_points=true, return_indexes=true)
    r_z_index = [k[1] for k in indexes]
    if length(r_z_index) == 0
        return Float64[], Float64[], Float64[]

    elseif length(r_z_index) == 1
        error("line_wall_2_wall: open field line should intersect wall at least twice.
               If it does not it's likely because the equilibrium grid was too small.")
    end

    # angle of incidence
    strike_angles = intersection_angles(r, z, wall_r, wall_z, indexes)

    if length(r_z_index) == 2
        # pass

    else
        # closest midplane point (favoring low field side)
        j0 = argmin(abs.(z .- Z0) .+ (r .< R0))
        # the closest intersection point (in steps) to z=Z0
        i1 = sortperm(abs.(r_z_index .- j0))[1]
        # the intersection on the other size of the midplane
        j1 = r_z_index[i1]
        if j0 < j1
            i2 = i1 - 1
        else
            i2 = i1 + 1
        end
        i = sort([i1, i2])
        r_z_index = r_z_index[i]
        crossings = crossings[i]
        strike_angles = strike_angles[i]
    end

    rr = vcat(crossings[1][1], r[r_z_index[1]+1:r_z_index[2]], crossings[2][1])
    zz = vcat(crossings[1][2], z[r_z_index[1]+1:r_z_index[2]], crossings[2][2])

    # sort clockwise (COCOS 11)
    if atan(zz[1] - Z0, rr[1] - R0) > atan(zz[end] - Z0, rr[end] - R0)
        rr = reverse(rr)
        zz = reverse(zz)
    end

    rr, zz, strike_angles
end

"""
    Bpol(a::T, κ::T, Ip::T) where {T<:Real}

Average poloidal magnetic field magnitude
"""
function Bpol(a::T, κ::T, Ip::T) where {T<:Real}
    return (constants.μ_0 * Ip) / (2π * a * sqrt((1.0 + κ^2) / 2.0))
end

"""
    Bpol_omp(eqt::IMAS.equilibrium__time_slice)

Poloidal magnetic field magnitude evaluated at the outer midplane
"""
function Bpol_omp(eqt::IMAS.equilibrium__time_slice)
    r, z, PSI_interpolant = ψ_interpolant(eqt)
    eq1d = eqt.profiles_1d
    R_omp = eq1d.r_outboard[end]
    Z_omp = eqt.global_quantities.magnetic_axis.z
    Br, Bz = Br_Bz_vector_interpolant(PSI_interpolant, [R_omp], [Z_omp])
    return sqrt(Br[1]^2.0 + Bz[1]^2.0)
end

"""
    power_sol(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d)

Total power coming out of the SOL [W]
"""
function power_sol(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d)
    tot_src = total_sources(core_sources, cp1d)
    p_sol = tot_src.electrons.power_inside[end] + tot_src.total_ion_power_inside[end]
    if p_sol < 0.0
        return 0.0
    else
        return p_sol
    end
end

function widthSOL_loarte(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    R0, B0 = vacuum_r0_b0(eqt)
    q95 = eqt.global_quantities.q95
    Psol = power_sol(core_sources, cp1d)
    return widthSOL_loarte(B0, q95, Psol)
end

function widthSOL_loarte(dd::IMAS.dd)
    return widthSOL_loarte(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

"""
    widthSOL_loarte(B0::T, q95::T, Psol::T) where {T<:Real}

Returns midplane power decay length λ_q in meters
"""
function widthSOL_loarte(B0::T, q95::T, Psol::T) where {T<:Real}
    return 0.00265 * Psol^0.38 * B0^(-0.71) * q95^0.3
end

function widthSOL_sieglin(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eq1d = eqt.profiles_1d
    R0 = (eq1d.r_outboard[end] .+ eq1d.r_inboard[end]) / 2.0
    a = (eq1d.r_outboard[end] .- eq1d.r_inboard[end])
    Psol = power_sol(core_sources, cp1d)
    ne_ped = interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(0.95)
    return widthSOL_sieglin(R0, a, Bpol_omp(eqt), Psol, ne_ped)
end

function widthSOL_sieglin(dd::IMAS.dd)
    return widthSOL_sieglin(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

"""
    widthSOL_sieglin(R0::T, a::T, Bpol_omp::T, Psol::T, ne_ped::T) where {T<:Real}

Returns integral power decay length λ_int in meters
Eich scaling(NF 53 093031) & B. Sieglin PPCF 55 (2013) 124039
"""
function widthSOL_sieglin(R0::T, a::T, Bpol_omp::T, Psol::T, ne_ped::T) where {T<:Real}
    λ_q = widthSOL_eich(R0, a, Bpol_omp, Psol)

    # From B. Sieglin PPCF 55 (2013) 124039
    # S includes the geometrical effects of the divertor assembly itself
    ne_ped /= 1.e19
    S = 0.09 * 1E-3 * ne_ped^1.02 * Bpol_omp^-1.01

    # extrapolate S from ASDEX results
    S *= (R0 / 1.65)

    # This is a valid approximation when S/λ_q < 10
    if S / λ_q > 10
        @warn "S/λ_q = $(S/λ_q) > 10 integral power decay approximation is inaccurate"
    end
    λ_int = λ_q + 1.64 * S

    return λ_int
end

function widthSOL_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eq1d = eqt.profiles_1d
    R0 = (eq1d.r_outboard[end] .+ eq1d.r_inboard[end]) / 2.0
    a = (eq1d.r_outboard[end] .- eq1d.r_inboard[end])
    Psol = power_sol(core_sources, cp1d)
    return widthSOL_eich(R0, a, Bpol_omp(eqt), Psol)
end

function widthSOL_eich(dd::IMAS.dd)
    return widthSOL_eich(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

"""
    widthSOL_eich(R0::T, a::T, Bpol_omp::T, Psol::T) where {T<:Real}

Returns midplane power decay length λ_q in meters
Eich scaling (NF 53 093031)
"""
function widthSOL_eich(R0::T, a::T, Bpol_omp::T, Psol::T) where {T<:Real}
    if Psol < 0.0
        return 0.0
    end
    λ_q = 1.35 * 1E-3 * (Psol / 1E6)^-0.02 * R0^0.04 * Bpol_omp^-0.92 * (a / R0)^0.42
    return λ_q
end

"""
    find_strike_points(wall_outline_r::T, wall_outline_z::T, pr::T, pz::T) where {T<:AbstractVector{<:Real}}

Finds strike points and angles of incidence between two paths
"""
function find_strike_points(wall_outline_r::T, wall_outline_z::T, pr::T, pz::T) where {T<:AbstractVector{<:Real}}
    indexes, pvx, pvy = intersection(wall_outline_r, wall_outline_z, pr, pz; as_list_of_points=false, return_indexes=true)
    angles = intersection_angles(wall_outline_r, wall_outline_z, pr, pz, indexes)
    return pvx, pvy, angles
end

"""
    find_strike_points(eqt::IMAS.equilibrium__time_slice, wall_outline_r::T, wall_outline_z::T) where {T<:AbstractVector{<:Real}}

Finds equilibrium strike points and angle of incidence between wall and strike leg
"""
function find_strike_points(eqt::IMAS.equilibrium__time_slice, wall_outline_r::T, wall_outline_z::T) where {T<:AbstractVector{<:Real}}
    Rx = Float64[]
    Zx = Float64[]
    θx = Float64[]

    private = flux_surface(eqt, eqt.profiles_1d.psi[end], false)
    for (pr, pz) in private
        pvx, pvy, angles = find_strike_points(wall_outline_r, wall_outline_z, pr, pz)
        append!(Rx, pvx)
        append!(Zx, pvy)
        append!(θx, angles)
    end

    return Rx, Zx, θx
end

"""
    find_strike_points!(eqt::IMAS.equilibrium__time_slice, dv::IMAS.divertors)

Adds strike points location to equilibrium IDS and the tilt_angle_pol in the divertors IDS
"""
function find_strike_points!(eqt::IMAS.equilibrium__time_slice, dv::IMAS.divertors)
    Rx = Float64[]
    Zx = Float64[]
    θx = Float64[]

    time = eqt.time

    for divertor in dv.divertor
        for target in divertor.target
            Rx0, Zx0, θx0 = find_strike_points(eqt, target.tile[1].surface_outline.r, target.tile[1].surface_outline.z)
            # allow for strike points to miss the divertors
            if isempty(Rx0)
                continue
            end
            push!(Rx, Rx0[1])
            push!(Zx, Zx0[1])
            push!(θx, θx0[1])
            set_time_array(target.tilt_angle_pol, :data, time, θx0[1])
        end
    end

    resize!(eqt.boundary_separatrix.strike_point, length(Rx))
    for (k, strike_point) in enumerate(eqt.boundary_separatrix.strike_point)
        strike_point.r = Rx[k]
        strike_point.z = Zx[k]
    end

    return Rx, Zx, θx
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall_outline_r::T, wall_outline_z::T) where {T<:AbstractVector{<:Real}}
    Rx, Zx, θx = find_strike_points(eqt, wall_outline_r, wall_outline_z)

    resize!(eqt.boundary_separatrix.strike_point, length(Rx))
    for (k, strike_point) in enumerate(eqt.boundary_separatrix.strike_point)
        strike_point.r = Rx[k]
        strike_point.z = Zx[k]
    end

    return Rx, Zx, θx
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice, bd::IMAS.build)
    wall_outline = get_build(bd, type=_plasma_).outline
    find_strike_points!(eqt, wall_outline.r, wall_outline.z)
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall)
    wall_outline = first_wall(wall)
    if wall_outline !== missing
        return find_strike_points!(eqt, wall_outline.r, wall_outline.z)
    end
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice)
    dd = top_dd(eqt)
    wall_outline = first_wall(dd.wall)
    if wall_outline !== missing
        return find_strike_points!(eqt, wall_outline.r, wall_outline.z)
    elseif !isempty(dd.build.layer)
        return find_strike_points!(eqt, dd.build)
    end
end
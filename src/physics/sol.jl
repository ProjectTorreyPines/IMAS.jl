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

@recipe function plot_OFL(OFL::Vector{OpenFieldLine})
    for ofl in OFL
        @series begin
            ofl
        end
    end
end

function sol(eq::IMAS.equilibrium, wall::IMAS.wall)
    return sol(eq, IMAS.first_wall(wall).r, IMAS.first_wall(wall).z)
end

function sol(eq::IMAS.equilibrium, wall_r::Vector{T}, wall_z::Vector{T}) where {T<:Real}
    r0 = eq.vacuum_toroidal_field.r0
    b0 = @ddtime(eq.vacuum_toroidal_field.b0)
    return sol(eq.time_slice[], r0, b0, wall_r, wall_z)
end

"""
    sol(eqt::IMAS.equilibrium__time_slice, r0::T, b0::T, wall_r::Vector{T}, wall_z::Vector{T}) where {T<:Real}

Trace open field lines up to wall
"""
function sol(eqt::IMAS.equilibrium__time_slice, r0::T, b0::T, wall_r::Vector{T}, wall_z::Vector{T}) where {T<:Real}
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    ############
    cc = IMAS.cocos(11)
    r, z, PSI_interpolant = ψ_interpolant(eqt)
    r_wall_midplane, _ = IMAS.intersection([R0, maximum(wall_r)], [Z0, Z0], wall_r, wall_z; as_list_of_points=false)
    psi_wall_midplane = PSI_interpolant.(r_wall_midplane, Z0)[1]
    psi__axis_level = eqt.profiles_1d.psi[1]
    psi__boundary_level = IMAS.find_psi_boundary(eqt; raise_error_on_not_open=true)
    psi_sign = sign(psi__boundary_level - psi__axis_level)
    ############

    # pack points near lcfs
    levels = psi__boundary_level .+ psi_sign .* 10.0 .^ LinRange(-3, log10(abs(psi_wall_midplane - psi__boundary_level)), 22)[1:end-1]

    OFL = OpenFieldLine[]
    for level in levels
        lines = IMAS.flux_surface(eqt, level, false)
        for line in lines
            rr, zz = line_wall_2_wall(line..., wall_r, wall_z, R0, Z0)
            if isempty(rr) || all(zz .> Z0) || all(zz .< Z0)
                continue
            end
            Br, Bz = Br_Bz_vector_interpolant(PSI_interpolant, cc, rr, zz)
            Bp = sqrt.(Br .^ 2.0 .+ Bz .^ 2.0)
            Bt = abs.(b0 .* r0 ./ rr)
            dp = sqrt.(gradient(rr) .^ 2.0 .+ gradient(zz) .^ 2.0)
            pitch = sqrt.(1.0 .+ (Bt ./ Bp) .^ 2)
            s = cumsum(pitch .* dp)
            midplane_index = argmin(abs.(zz .- Z0) .+ (rr .< R0))
            s = abs.(s .- s[midplane_index])
            push!(OFL, OpenFieldLine(rr, zz, Br, Bz, Bp, Bt, pitch, s, midplane_index))
        end
    end
    return OFL
end

"""
    line_wall_2_wall(
        r::AbstractVector{T},
        z::AbstractVector{T},
        wall_r::AbstractVector{T},
        wall_z::AbstractVector{T},
        R0::Real,
        Z0::Real) where {T<:Real}

Returns r, z coordinates of open field line contained within wall
"""
function line_wall_2_wall(
    r::AbstractVector{T},
    z::AbstractVector{T},
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T},
    R0::Real,
    Z0::Real) where {T<:Real}

    indexes, crossings = IMAS.intersection(r, z, wall_r, wall_z; as_list_of_points=true, return_indexes=true)
    indexes = [k[1] for k in indexes]
    if length(indexes) == 0
        return [], []

    elseif length(indexes) == 1
        error("line_wall_2_wall: open field line should intersect wall at least twice.
               If it does not it's likely because the equilibrium grid was too small.")

    elseif length(indexes) == 2
        # pass

    else
        # closest midplane point (favoring low field side)
        j0 = argmin(abs.(z .- Z0) .+ (r .< R0))
        # the closest intersection point (in steps) to z=Z0
        i1 = sortperm(abs.(indexes .- j0))[1]
        # the intersection on the other size of the midplane
        j1 = indexes[i1]
        if j0 < j1
            i2 = i1 - 1
        else
            i2 = i1 + 1
        end
        i = sort([i1, i2])
        indexes = indexes[i]
        crossings = crossings[i]
    end

    rr = vcat(crossings[1][1], r[indexes[1]+1:indexes[2]], crossings[2][1])
    zz = vcat(crossings[1][2], z[indexes[1]+1:indexes[2]], crossings[2][2])

    # sort clockwise (COCOS 11)
    if atan(zz[1] - Z0, rr[1] - R0) > atan(zz[end] - Z0, rr[end] - R0)
        rr = reverse(rr)
        zz = reverse(zz)
    end

    rr, zz
end

"""
    Bpol(a::T, κ::T, Ip::T) where {T<:Real}

Average poloidal magnetic field magnitude
"""
function Bpol(a::T, κ::T, Ip::T) where {T<:Real}
    return (constants.μ_0 * Ip) / (2π * a * sqrt((1.0 + κ^2) / 2.0))
end

"""
    power_sol(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d)

Total power coming out of the SOL [W]
"""
function power_sol(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d)
    tot_src = IMAS.total_sources(core_sources, cp1d)
    return tot_src.electrons.power_inside[end] + tot_src.total_ion_power_inside[end]
end

function widthSOL_sieglin(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eq1d = eqt.profiles_1d
    R0 = (eq1d.r_outboard[end] .+ eq1d.r_inboard[end]) / 2.0
    a = (eq1d.r_outboard[end] .- eq1d.r_inboard[end])
    κ = eq1d.elongation[end]
    Ip = eqt.global_quantities.ip
    Psol = power_sol(core_sources, cp1d)
    ne_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(0.95)
    return widthSOL_sieglin(R0, a, κ, Ip, Psol, ne_ped)
end

function widthSOL_sieglin(dd::IMAS.dd)
    return widthSOL_sieglin(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

"""
    widthSOL_sieglin(R0::T, a::T, κ::T, Ip::T, Psol::T, ne_ped::T) where {T<:Real}

Returns integral power decay length λ_int in meters
Eich scaling(NF 53 093031) & B. Sieglin PPCF 55 (2013) 124039
"""
function widthSOL_sieglin(R0::T, a::T, κ::T, Ip::T, Psol::T, ne_ped::T) where {T<:Real}
    λ_q = widthSOL_eich(R0, a, κ, Ip, Psol)

    # From B. Sieglin PPCF 55 (2013) 124039
    # S includes the geometrical effects of the divertor assembly itself
    ne_ped /= 1.e19
    S = 0.09 * 1E-3 * ne_ped^1.02 * Bpol(a, κ, Ip)^-1.01

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
    κ = eq1d.elongation[end]
    Ip = eqt.global_quantities.ip
    Psol = power_sol(core_sources, cp1d)
    return widthSOL_eich(R0, a, κ, Ip, Psol)
end

function widthSOL_eich(dd::IMAS.dd)
    return widthSOL_eich(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

"""
    widthSOL_eich(R0::T, a::T, κ::T, Ip::T, Psol::T) where {T<:Real}

Returns midplane power decay length λ_q in meters
Eich scaling (NF 53 093031)
"""
function widthSOL_eich(R0::T, a::T, κ::T, Ip::T, Psol::T) where {T<:Real}
    if Psol < 0.0
        return 0.0
    end
    Psol /= 1E6
    λ_q = 1.35 * 1E-3 * Psol^-0.02 * R0^0.04 * Bpol(a, κ, Ip)^-0.92 * (a / R0)^0.42
    return λ_q
end
"""
    find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall_outline_r::T, wall_outline_z::T) where {T<:AbstractVector{<:Real}}

Finds equilibrium strike points and adds them to eqt.boundary_separatrix.strike_point
"""
function find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall_outline_r::T, wall_outline_z::T) where {T<:AbstractVector{<:Real}}
    Rx = Float64[]
    Zx = Float64[]

    private = IMAS.flux_surface(eqt, eqt.profiles_1d.psi[end], false)
    for (pr, pz) in private
        pvx, pvy = IMAS.intersection(wall_outline_r, wall_outline_z, pr, pz; as_list_of_points=false)
        append!(Rx, pvx)
        append!(Zx, pvy)
    end

    resize!(eqt.boundary_separatrix.strike_point, length(Rx))
    for (k, strike_point) in enumerate(eqt.boundary_separatrix.strike_point)
        strike_point.r = Rx[k]
        strike_point.z = Zx[k]
    end

    return Rx, Zx
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice, bd::IMAS.build)
    wall_outline = IMAS.get_build(bd, type=_plasma_).outline
    return find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall_outline.r, wall_outline.z)
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall)
    wall_outline = IMAS.first_wall(wall)
    if wall_outline !== missing
        return find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall_outline.r, wall_outline.z)
    end
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice)
    dd = IMAS.top_dd(eqt)
    return find_strike_points!(eqt, dd.wall)
end
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
            if isempty(rr) || all(zz.>Z0) || all(zz.<Z0)
                continue
            end
            Br, Bz = Br_Bz_vector_interpolant(PSI_interpolant, cc, rr, zz)
            Bp = sqrt.(Br .^ 2.0 .+ Bz .^ 2.0)
            Bt = abs.(b0 .* r0 ./ rr)
            dp = sqrt.(IMAS.gradient(rr) .^ 2.0 .+ IMAS.gradient(zz) .^ 2.0)
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
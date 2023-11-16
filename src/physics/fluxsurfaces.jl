using LinearAlgebra

"""
    ψ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)

Returns r, z, and ψ interpolant
"""
function ψ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    grid_type = identifier_name(eqt2d.grid_type, :rectangular)
    if grid_type == :rectangular
        r = range(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end]; length=length(eqt2d.grid.dim1))
        z = range(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end]; length=length(eqt2d.grid.dim2))
        return r, z, Interpolations.cubic_spline_interpolation((r, z), eqt2d.psi)
    else
        error("ψ_interpolant cannot handle grid of type `$grid_type`")
    end
end

function ψ_interpolant(dd::IMAS.dd)
    return ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d[])
end

"""
    Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::Array{T}, z::Array{T}) where {T<:Real}

Returns Br and Bz tuple evaluated at r and z starting from ψ interpolant
"""
function Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::Array{T}, z::Array{T}) where {T<:Real}
    # Check that r and z are the same size
    @assert size(r) == size(z)
    grad = [Interpolations.gradient(PSI_interpolant, r[idx], z[idx]) for idx in CartesianIndices(r)]
    Br = [grad[idx][2] / r[idx] / (2π) for idx in CartesianIndices(r)]
    Bz = [-grad[idx][1] / r[idx] / (2π) for idx in CartesianIndices(r)]
    return Br, Bz
end

function Br_Bz(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    Z, R = meshgrid(z, r)
    return Br_Bz(PSI_interpolant, R, Z)
end

"""
    Bp(PSI_interpolant::Interpolations.AbstractInterpolation, r::Vector{T}, z::Vector{T}) where {T<:Real}

Returns Bp evaluated at r and z starting from ψ interpolant
"""
function Bp(PSI_interpolant::Interpolations.AbstractInterpolation, r::Vector{T}, z::Vector{T}) where {T<:Real}
    Br, Bz = Br_Bz(PSI_interpolant, r, z)
    return sqrt.(Br .^ 2.0 .+ Bz .^ 2.0)
end

function Bp(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    Z, R = meshgrid(z, r)
    return Bp(PSI_interpolant, R, Z)
end

"""
    find_psi_boundary(eqt::IMAS.equilibrium__time_slice; precision::Float64=1e-6, raise_error_on_not_open::Bool=true, raise_error_on_not_closed::Bool=true)

Find psi value of the last closed flux surface
"""
function find_psi_boundary(eqt::IMAS.equilibrium__time_slice; precision::Float64=1e-6, raise_error_on_not_open::Bool=true, raise_error_on_not_closed::Bool=true)
    dim1 = eqt.profiles_2d[1].grid.dim1
    dim2 = eqt.profiles_2d[1].grid.dim2
    PSI = eqt.profiles_2d[1].psi
    psi = eqt.profiles_1d.psi
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z
    return find_psi_boundary(dim1, dim2, PSI, psi, R0, Z0; precision, raise_error_on_not_open, raise_error_on_not_closed)
end

function find_psi_boundary(dd::IMAS.dd; precision::Float64=1e-6, raise_error_on_not_open::Bool=true, raise_error_on_not_closed::Bool=true)
    return find_psi_boundary(dd.equilibrium.time_slice[]; precision, raise_error_on_not_open, raise_error_on_not_closed)
end

function find_psi_boundary(
    dim1::Union{AbstractVector{T},AbstractRange{T}},
    dim2::Union{AbstractVector{T},AbstractRange{T}},
    PSI::Matrix{T},
    psi::Union{AbstractVector{T},AbstractRange{T}},
    R0::T,
    Z0::T;
    precision::Float64=1e-6,
    raise_error_on_not_open::Bool,
    raise_error_on_not_closed::Bool) where {T<:Real}

    psirange_init = [psi[1] * 0.9 + psi[end] * 0.1, psi[end] + 0.5 * (psi[end] - psi[1])]

    # innermost tentative flux surface (which should be closed!)
    pr, pz = flux_surface(dim1, dim2, PSI, psi, R0, Z0, psirange_init[1], true)
    if isempty(pr)
        if raise_error_on_not_closed
            error("Flux surface at ψ=$(psirange_init[1]) is not closed; ψ=[$(psi[1])...$(psi[end])]")
        else
            return nothing
        end
    end

    # outermost tentative flux surface (which should be open!)
    pr, pz = flux_surface(dim1, dim2, PSI, psi, R0, Z0, psirange_init[end], true)
    if length(pr) > 0
        if raise_error_on_not_open
            error("Flux surface at ψ=$(psirange_init[end]) is not open; ψ=[$(psi[1])...$(psi[end])]")
        else
            return nothing
        end
    end

    δd = sqrt((dim1[2] - dim1[1])^2 + (dim2[2] - dim2[1])^2)
    psirange = deepcopy(psirange_init)
    for k in 1:100
        psimid = (psirange[1] + psirange[end]) / 2.0
        pr, pz = flux_surface(dim1, dim2, PSI, psi, R0, Z0, psimid, true)
        # closed flux surface
        if length(pr) > 0
            psirange[1] = psimid
            if (abs(psirange[end] - psirange[1]) / abs(psirange[end] + psirange[1]) / 2.0) < precision
                if any(abs.([(minimum(pr) - minimum(dim1)), (maximum(pr) - maximum(dim1)), (minimum(pz) - minimum(dim2)), (maximum(pz) - maximum(dim2))]) .< 2 * δd)
                    return psi[end]
                else
                    return psimid
                end
            end
            # open flux surface
        else
            psirange[end] = psimid
        end
    end

    return error("Could not find closed boundary between ψ=$(psirange_init[1]) and ψ=$(psirange_init[end])")
end

"""
    find_psi_2nd_separatrix(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation) 

    Returns ψ of second magentic separatrix
    This function is needed if eqt.boundary_secondary_separatrix is empty
"""
function find_psi_2nd_separatrix(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation) 

    null2 = find_2nd_x_point!(eqt, PSI_interpolant)
    psi_2nd_separatix = PSI_interpolant.(null2[1], null2[2]) # psi of second separatrix

    return psi_2nd_separatix
end

function find_psi_2nd_separatrix(dd::IMAS.dd)
        rr, zz, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d[])
    return find_psi_2nd_separatrix(dd.equilibrium.time_slice[],  PSI_interpolant)
end

"""
    find_psi_last_diverted(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{<:Real}, wall_z::Vector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation)

    Returns [psi_up, psi_low] of the two flux surfaces around the last diverted flux surface. 
    psi_up will be the first surface inside OFL[:lfs_far]; psi_low will be the last surface inside OFL[:lfs]
    Precision between the two is defined on the poloidal crossection area at the OMP (Psol*precision = power flowing between psi_up and psi_low ~ 0)
 """
function find_psi_last_diverted(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{<:Real}, wall_z::Vector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation; precision::Float64=1e-7)

    RA = eqt.global_quantities.magnetic_axis.r # R of magnetic axis
    ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis

    Xpoint2 = collect(find_2nd_x_point!(eqt,PSI_interpolant))
    sign_z = sign(Xpoint2[2]) # sign of Z coordinate of 2nd null
    
    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt,PSI_interpolant)         # psi second magnetic separatrix
    psi_separatrix    = find_psi_boundary(eqt; raise_error_on_not_open=true) # psi LCFS
    
    #intersect 2nd separatrix with wall, and look 
    surface = flux_surface(eqt, psi_2ndseparatrix*1.0001, false)
    r_intersect = Float64[]
    z_intersect = Float64[]
    for (r,z) in surface
        rr, zz, strike_angles = line_wall_2_wall(r, z, wall_r, wall_z, RA, ZA) 
        if isempty(rr) || all(zz .< ZA)
            continue
        end
        #save intersections with wall 
        push!(r_intersect, rr[1], rr[end])
        push!(z_intersect, zz[1], zz[end])
    end

    #take intersections above midplane
    r_intersect = r_intersect[z_intersect.>ZA]
    z_intersect = z_intersect[z_intersect.>ZA]
    order = sortperm(r_intersect)
    z_intersect = z_intersect[order]
    r_intersect = r_intersect[order]

    # check if upper null is inside the wall, by checking if upper null is left/right of the vector between the 2 (ordered) intersections
    # This is an approximation (should work except for exotic walls)
    vec  = [diff(r_intersect)[1], diff(z_intersect)[1]]
    vec2 = Xpoint2 - [r_intersect[1], z_intersect[1]]
    if vec[1]*vec2[2]-vec[2]*vec2[1]>0
        #upper null on the left = is outside wall
        null_is_inside = false
        else
        #upper null on the right = is inside wall
        null_is_inside = true
    end
    # find the two surfaces psi_up (inside OFL[:lfs_far])  and psi_low (inside OFL[:lfs]) around the last diverted flux surface
    counter_max = 50
    counter =0
    psi__axis_level = eqt.profiles_1d.psi[1] # psi value on axis 
    psi_sign = sign(psi_separatrix - psi__axis_level) # sign of the poloidal flux taking psi_axis = 0
    psi_wall = PSI_interpolant.(wall_r,wall_z)
    if psi_sign >0
        psi_up  = minimum([psi_2ndseparatrix+psi_sign, maximum(psi_wall)*0.999]) # increase value just to be sure of being inside OFL[:lfs_far]
    else
        psi_up  = maximum([psi_2ndseparatrix+psi_sign, minimum(psi_wall)*1.001]) # increase value just to be sure of being inside OFL[:lfs_far]
    end

    psi_low = psi_separatrix
    psi = (psi_up+psi_low)/2
    err = 1
    while abs(err)> precision && counter< counter_max
        surface = flux_surface(eqt, psi, false)
        for (r,z) in surface
            rr, zz, strike_angles = line_wall_2_wall(r, z, wall_r, wall_z, RA, ZA) 
            
            if isempty(rr) || all(zz .> ZA) || all(zz .< ZA)
                continue
            end

            if sign_z*zz[1]>0|| sign_z*zz[end]>0
                # psi intersects top FW -> update upper bound
                psi_up = psi
            else
                #psi intersects divertor -> update lower bound
                psi_low = psi
            end
        end
        
        # better to compute error on poloidal area between [psi_low, psi_up] (needed for accurate power balance)
        r_up  = find_r_midplane_from_ψ(eqt,psi_up)
        r_low = find_r_midplane_from_ψ(eqt,psi_low)
 
        A = π*(r_up^2-r_low^2) # annular area between r_up and r_low [m^2] 
        err = abs(A)
        psi = (psi_up+psi_low)/2
        
        counter = counter +1
    end
    
    return [psi_low, psi_up], null_is_inside # return both psi_up and psi_low to increase resolution around last diverted flux surface
end

function find_psi_last_diverted(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall, PSI_interpolant::Interpolations.AbstractInterpolation; precision::Float64=1e-7)
    return find_psi_last_diverted(eqt, first_wall(wall).r, first_wall(wall).z, PSI_interpolant; precision)
end 

function find_psi_last_diverted(dd::IMAS.dd; precision::Float64=1e-7) 
    rr, zz, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d[])
    return find_psi_last_diverted(dd.equilibrium.time_slice[], dd.wall, PSI_interpolant; precision)
end

"""
    find_ψ_from_r_midplane(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation,r::T) where {T<:Real}

    Returns ψ(r): the poloidal flux ψ of the flux surface that intersects the midplane at the coordinate r (major radius)
"""
function find_ψ_from_r_midplane(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation, r::T ) where {T<:Real}
    psi__boundary_level   = find_psi_boundary(eqt; raise_error_on_not_open=true) # psi magnetic separatrix
    psi__axis_level = eqt.profiles_1d.psi[1] # psi value on axis 
    psi_sign = sign(psi__boundary_level - psi__axis_level) # sign of the poloidal flux taking psi_axis = 0
    r_max_grid             = maximum(eqt.profiles_2d[].grid.dim1)                 # max r on grid of equilibrium
    @assert r < r_max_grid
    RA = eqt.global_quantities.magnetic_axis.r # R of magnetic axis
    @assert r > RA
    ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis
    
    psi = PSI_interpolant.(r, ZA)[1]           # compute psi at r on midplane (first guess)

    # PSI_interpolant is not as precise as intersecting flux_surface with midplane (as is done in find_r_midplane_from_ψ)
    find_psi = [psi-0.01, psi+0.01]
    r_of_find_psi = find_r_midplane_from_ψ(eqt, find_psi) 
    psi = find_psi[1] + diff(find_psi)[1]/diff(r_of_find_psi)[1]*(r-r_of_find_psi[1]) #linear interp of find_r_midplane_from_ψ
    if (psi - psi__boundary_level) <= psi_sign * psi__boundary_level*0.000001
        # this is needed in find_r_midplane_from_ψ. This makes the two function work consistently r(ψ(R)) = R
        psi = psi__boundary_level + psi_sign * psi__boundary_level*0.000001 # if psi exactly psi__boundary_level flux surfaces does not find separatrix
    end
    return psi
end

function find_ψ_from_r_midplane(dd::IMAS.dd,r::T) where {T<:Real}
    rr, zz, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d[])
    return find_ψ_from_r_midplane(dd.equilibrium.time_slice[], PSI_interpolant, r)
end

function find_ψ_from_r_midplane(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation, r::T ) where {T<:AbstractVector{<:Real}}
    # if r is a vector
    PSI = similar(r)

    for index in 1:length(PSI)
        PSI[index] = find_ψ_from_r_midplane(eqt, PSI_interpolant, r[index])
    end
    return PSI
end

function find_ψ_from_r_midplane(dd::IMAS.dd, r::T ) where {T<:AbstractVector{<:Real}}
    rr, zz, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d[])
    return find_ψ_from_r_midplane(dd.equilibrium.time_slice[], PSI_interpolant, r)
end

"""
    find_r_midplane_from_ψ(eqt::IMAS.equilibrium__time_slice,psi::T) where {T<:Real}
    or
    find_r_midplane_from_ψ(eqt::IMAS.equilibrium__time_slice, psi::T) where {T<:AbstractVector{<:Real}}

    Returns r(ψ): the r coordinate of the intersection between the flux surface of defined poloidal flux ψ and midplane
"""
function find_r_midplane_from_ψ(eqt::IMAS.equilibrium__time_slice, psi::T) where {T<:Real}
    
    psi__boundary_level = find_psi_boundary(eqt; raise_error_on_not_open=true) #psi on separatrix
    psi__axis_level = eqt.profiles_1d.psi[1] # psi value on axis 
    psi_sign = sign(psi__boundary_level - psi__axis_level) # sign of the poloidal flux taking psi_axis = 0

    RA = eqt.global_quantities.magnetic_axis.r # R of magnetic axis
    ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis
    psi = float(psi)

    if psi > psi__boundary_level # check if surface is open or closed
        is_closed = false
    else
        is_closed = true
    end

    surface   = flux_surface(eqt, psi, is_closed) #returns (r,z) of surfaces with psi;
    r_midplane = Float64[] 
    if is_closed 
        # when closed = true, flux_surface returns Tuple{Vector{Float64}, Vector{Float64}, Float64}
        ind, crossings = intersection([RA, 100], [ZA, ZA], surface[1], surface[2]) #cross surface with OMP - segment from RA to 100 m
        r_midplane = crossings[1][1]
    else
        for (r,z) in surface
            # when closed = false, flux_surface returns Vector{Tuple{Vector{Float64}, Vector{Float64}}}
            ind, crossings = intersection([RA, 100], [ZA, ZA], r, z) #cross surface with OMP - segment from RA to 100 m
            if isempty(ind) # check if intresection exists
                continue
            end
            push!(r_midplane,crossings[1][1])       # if so save value
            r_midplane = r_midplane[1] # make it Float64 instead of Vector{1,Float64}
        end
    end
    
    return r_midplane
end

function find_r_midplane_from_ψ(dd::IMAS.dd, psi::T) where {T<:Real}
    return find_r_midplane_from_ψ(dd.equilibrium.time_slice[], psi)
end

function find_r_midplane_from_ψ(eqt::IMAS.equilibrium__time_slice, psi::T) where {T<:AbstractVector{<:Real}}
    # if psi is a vector
    R_midplane = similar(psi)
    for index in 1:length(psi)
        R_midplane[index] = find_r_midplane_from_ψ(eqt,psi[index])
    end
    return R_midplane
end

function find_r_midplane_from_ψ(dd::IMAS.dd, psi::T) where {T<:AbstractVector{<:Real}}
    return find_r_midplane_from_ψ(dd.equilibrium.time_slice[], psi)
end


"""
    flux_surfaces(eq::equilibrium; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities in the equilibrium IDS
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eq::equilibrium; upsample_factor::Int=1)
    for time_index in eachindex(eq.time_slice)
        flux_surfaces(eq.time_slice[time_index]; upsample_factor)
    end
    return eq
end

"""
    flux_surfaces(eqt::equilibrium__time_slice; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eqt::equilibrium__time_slice; upsample_factor::Int=1)
    R0, B0 = vacuum_r0_b0(eqt)
    return flux_surfaces(eqt, B0, R0; upsample_factor)
end

function find_magnetic_axis!(r::AbstractVector{<:Real}, z::AbstractVector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation, psi_sign::Real)
    res = Optim.optimize(
        x -> begin
            try
                PSI_interpolant(x[1], x[2]) * psi_sign
            catch e
                if typeof(e) <: BoundsError
                    return Inf
                else
                    rethrow(e)
                end
            end
        end,
        [r[Int(round(length(r) / 2))], z[Int(round(length(z) / 2))]],
        Optim.Newton(),
        Optim.Options(; g_tol=1E-8);
        autodiff=:forward
    )
    return res.minimizer[1], res.minimizer[2]
end

"""
    flux_surfaces(eqt::equilibrium__time_slice{T}, B0::T, R0::T; upsample_factor::Int=1) where {T<:Real}

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice, B0 and R0
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eqt::equilibrium__time_slice{T}, B0::T, R0::T; upsample_factor::Int=1) where {T<:Real}
    r, z, PSI_interpolant = ψ_interpolant(eqt.profiles_2d[1])
    PSI = eqt.profiles_2d[1].psi

    # upsampling for high-resolution r,z flux surface coordinates
    if upsample_factor > 1
        r = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end]; length=length(eqt.profiles_2d[1].grid.dim1) * upsample_factor)
        z = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end]; length=length(eqt.profiles_2d[1].grid.dim2) * upsample_factor)
        PSI = PSI_interpolant(r, z)
    end

    psi_sign = sign(eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1])

    # find magnetic axis
    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = find_magnetic_axis!(r, z, PSI_interpolant, psi_sign)
    psi_axis = PSI_interpolant(eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z)
    eqt.profiles_1d.psi =
        (eqt.profiles_1d.psi .- eqt.profiles_1d.psi[1]) ./ (eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1]) .* (eqt.profiles_1d.psi[end] - psi_axis) .+ psi_axis

    for item in (
        :b_field_average,
        :b_field_max,
        :b_field_min,
        :elongation,
        :triangularity_lower,
        :triangularity_upper,
        :squareness_lower_inner,
        :squareness_lower_outer,
        :squareness_upper_inner,
        :squareness_upper_outer,
        :r_inboard,
        :r_outboard,
        :q,
        :surface,
        :dvolume_dpsi,
        :j_tor,
        :area,
        :volume,
        :gm1,
        :gm2,
        :gm4,
        :gm5,
        :gm8,
        :gm9,
        :phi,
        :trapped_fraction
    )
        setproperty!(eqt.profiles_1d, item, zeros(eltype(eqt.profiles_1d.psi), size(eqt.profiles_1d.psi)))
    end

    PR = Vector{T}[]
    PZ = Vector{T}[]
    LL = Vector{T}[]
    FLUXEXPANSION = Vector{T}[]
    INT_FLUXEXPANSION_DL = zeros(T, length(eqt.profiles_1d.psi))
    BPL = zeros(T, length(eqt.profiles_1d.psi))
    for (k, psi_level0) in reverse!(collect(enumerate(eqt.profiles_1d.psi)))

        if k == 1 # on axis flux surface is a synthetic one
            eqt.profiles_1d.elongation[1] = eqt.profiles_1d.elongation[2] - (eqt.profiles_1d.elongation[3] - eqt.profiles_1d.elongation[2])
            eqt.profiles_1d.triangularity_upper[1] = 0.0
            eqt.profiles_1d.triangularity_lower[1] = 0.0

            a = (eqt.profiles_1d.r_outboard[2] - eqt.profiles_1d.r_inboard[2]) / 100.0
            b = eqt.profiles_1d.elongation[1] * a

            t = range(0, 2π; length=17)
            pr = cos.(t) .* a .+ eqt.global_quantities.magnetic_axis.r
            pz = sin.(t) .* b .+ eqt.global_quantities.magnetic_axis.z

            # Extrema on array indices
            (imaxr, iminr, imaxz, iminz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = fluxsurface_extrema(pr, pz)

        else  # other flux surfaces
            # trace flux surface
            pr, pz, psi_level =
                flux_surface(r, z, PSI, eqt.profiles_1d.psi, eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z, psi_level0, true)
            if isempty(pr)
                # p = heatmap(r, z, PSI'; colorbar=true, aspect_ratio=:equal)
                # contour!(r, z, PSI'; color=:white, levels=100)
                # contour!(r, z, PSI'; levels=[eqt.profiles_1d.psi[end]], color=:white, lw=2)
                # display(p)
                error("IMAS: Could not trace closed flux surface $k out of $(length(eqt.profiles_1d.psi)) at ψ = $(psi_level)")
            end

            # Extrema on array indices
            (imaxr, iminr, imaxz, iminz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = fluxsurface_extrema(pr, pz)

            # accurate geometric quantities by finding geometric extrema as optimization problem
            w = 1E-4 # push away from magnetic axis
            function fx(x::AbstractVector{<:Real}, psi_level::Float64, eqt::IMAS.equilibrium__time_slice, w::Float64)
                try
                    (PSI_interpolant(x[1], x[2]) - psi_level)^2 - (x[1] - eqt.global_quantities.magnetic_axis.r)^2 * w
                catch
                    return 100
                end
            end
            function fz(x::AbstractVector{<:Real}, psi_level::Float64, eqt::IMAS.equilibrium__time_slice, w::Float64)
                try
                    (PSI_interpolant(x[1], x[2]) - psi_level)^2 - (x[2] - eqt.global_quantities.magnetic_axis.z)^2 * w
                catch
                    return 100
                end
            end
            res = Optim.optimize(x -> fx(x, psi_level, eqt, w), [max_r, z_at_max_r], Optim.Newton(), Optim.Options(; g_tol=1E-8); autodiff=:forward)
            (max_r, z_at_max_r) = (res.minimizer[1], res.minimizer[2])
            res = Optim.optimize(x -> fx(x, psi_level, eqt, w), [min_r, z_at_min_r], Optim.Newton(), Optim.Options(; g_tol=1E-8); autodiff=:forward)
            (min_r, z_at_min_r) = (res.minimizer[1], res.minimizer[2])
            if psi_level0 != eqt.profiles_1d.psi[end]
                res = Optim.optimize(x -> fz(x, psi_level, eqt, w), [r_at_max_z, max_z], Optim.Newton(), Optim.Options(; g_tol=1E-8); autodiff=:forward)
                (r_at_max_z, max_z) = (res.minimizer[1], res.minimizer[2])
                res = Optim.optimize(x -> fz(x, psi_level, eqt, w), [r_at_min_z, min_z], Optim.Newton(), Optim.Options(; g_tol=1E-8); autodiff=:forward)
                (r_at_min_z, min_z) = (res.minimizer[1], res.minimizer[2])
            end
            # p = plot(pr, pz, label = "")
            # plot!([max_r], [z_at_max_r], marker = :cicle)
            # plot!([min_r], [z_at_min_r], marker = :cicle)
            # plot!([r_at_max_z], [max_z], marker = :cicle)
            # plot!([r_at_min_z], [min_z], marker = :cicle)
            # display(p)

            # plasma boundary information
            if k == length(eqt.profiles_1d.psi)
                eqt.boundary.outline.r = pr
                eqt.boundary.outline.z = pz
            end
        end

        eqt.profiles_1d.r_outboard[k] = max_r
        eqt.profiles_1d.r_inboard[k] = min_r

        # miller geometric coefficients
        _, _, κ, δu, δl, ζou, ζol, ζil, ζiu = miller_R_a_κ_δ_ζ(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
        eqt.profiles_1d.elongation[k] = κ
        eqt.profiles_1d.triangularity_upper[k] = δu
        eqt.profiles_1d.triangularity_lower[k] = δl
        eqt.profiles_1d.squareness_lower_outer[k] = ζol
        eqt.profiles_1d.squareness_upper_outer[k] = ζou
        eqt.profiles_1d.squareness_lower_inner[k] = ζil
        eqt.profiles_1d.squareness_upper_inner[k] = ζiu

        # poloidal magnetic field (with sign)
        Br, Bz = Br_Bz(PSI_interpolant, pr, pz)
        Bp2 = Br .^ 2.0 .+ Bz .^ 2.0
        Bp_abs = sqrt.(Bp2)
        Bp = (
            Bp_abs .*
            sign.((pz .- eqt.global_quantities.magnetic_axis.z) .* Br .- (pr .- eqt.global_quantities.magnetic_axis.r) .* Bz)
        )

        # flux expansion
        dl = vcat(0.0, sqrt.(diff(pr) .^ 2 + diff(pz) .^ 2))
        ll = cumsum(dl)
        fluxexpansion = 1.0 ./ Bp_abs
        int_fluxexpansion_dl = integrate(ll, fluxexpansion)
        Bpl = integrate(ll, Bp)

        # save flux surface coordinates for later use
        pushfirst!(PR, pr)
        pushfirst!(PZ, pz)
        pushfirst!(LL, ll)
        pushfirst!(FLUXEXPANSION, fluxexpansion)
        INT_FLUXEXPANSION_DL[k] = int_fluxexpansion_dl
        BPL[k] = Bpl

        # trapped fraction
        Bt = eqt.profiles_1d.f[k] ./ pr
        Btot = sqrt.(Bp2 .+ Bt .^ 2)
        Bmin = minimum(Btot)
        Bmax = maximum(Btot)
        Bratio = Btot ./ Bmax
        avg_Btot = flxAvg(Btot, ll, fluxexpansion, int_fluxexpansion_dl)
        avg_Btot2 = flxAvg(Btot .^ 2, ll, fluxexpansion, int_fluxexpansion_dl)
        hf = flxAvg((1.0 .- sqrt.(1.0 .- Bratio) .* (1.0 .+ Bratio ./ 2.0)) ./ Bratio .^ 2, ll, fluxexpansion, int_fluxexpansion_dl)
        h = avg_Btot / Bmax
        h2 = avg_Btot2 / Bmax^2
        ftu = 1.0 - h2 / (h^2) * (1.0 - sqrt(1.0 - h) * (1.0 + 0.5 * h))
        ftl = 1.0 - h2 * hf
        eqt.profiles_1d.trapped_fraction[k] = 0.75 * ftu + 0.25 * ftl

        # Bavg
        eqt.profiles_1d.b_field_average[k] = avg_Btot

        # Bmax
        eqt.profiles_1d.b_field_max[k] = Bmax

        # Bmin
        eqt.profiles_1d.b_field_min[k] = Bmin

        # gm1 = <1/R^2>
        eqt.profiles_1d.gm1[k] = flxAvg(1.0 ./ pr .^ 2, ll, fluxexpansion, int_fluxexpansion_dl)

        # gm4 = <1/B^2>
        eqt.profiles_1d.gm4[k] = flxAvg(1.0 ./ Btot .^ 2, ll, fluxexpansion, int_fluxexpansion_dl)

        # gm5 = <B^2>
        eqt.profiles_1d.gm5[k] = avg_Btot2

        # gm8 = <R>
        eqt.profiles_1d.gm8[k] = flxAvg(pr, ll, fluxexpansion, int_fluxexpansion_dl)

        # gm9 = <1/R>
        eqt.profiles_1d.gm9[k] = flxAvg(1.0 ./ pr, ll, fluxexpansion, int_fluxexpansion_dl)

        # j_tor = <j_tor/R> / <1/R>
        eqt.profiles_1d.j_tor[k] =
            (
                -(eqt.profiles_1d.dpressure_dpsi[k] + eqt.profiles_1d.f_df_dpsi[k] * eqt.profiles_1d.gm1[k] / constants.μ_0) *
                (2π)
            ) / eqt.profiles_1d.gm9[k]

        # dvolume_dpsi
        eqt.profiles_1d.dvolume_dpsi[k] = (sign(flxAvg(Bp, ll, fluxexpansion, int_fluxexpansion_dl)) * int_fluxexpansion_dl)

        # surface area
        eqt.profiles_1d.surface[k] = 2π * sum(pr .* dl)

        # q
        eqt.profiles_1d.q[k] = eqt.profiles_1d.dvolume_dpsi[k] .* eqt.profiles_1d.f[k] .* eqt.profiles_1d.gm1[k] ./ (2π)

        # quantities calculated on the last closed flux surface
        if k == length(eqt.profiles_1d.psi)
            # perimeter
            eqt.global_quantities.length_pol = ll[end]
        end
    end

    # integral quantities
    for k in 2:length(eqt.profiles_1d.psi)
        # area
        eqt.profiles_1d.area[k] = integrate(eqt.profiles_1d.psi[1:k], eqt.profiles_1d.dvolume_dpsi[1:k] .* eqt.profiles_1d.gm9[1:k]) ./ 2π

        # volume
        eqt.profiles_1d.volume[k] = integrate(eqt.profiles_1d.psi[1:k], eqt.profiles_1d.dvolume_dpsi[1:k])

        # phi
        eqt.profiles_1d.phi[k] = integrate(eqt.profiles_1d.psi[1:k], eqt.profiles_1d.q[1:k])
    end

    # rho_tor_norm
    rho = sqrt.(abs.(eqt.profiles_1d.phi ./ (π * B0)))
    rho_meters = rho[end]
    eqt.profiles_1d.rho_tor = rho
    eqt.profiles_1d.rho_tor_norm = rho ./ rho_meters

    # phi 2D
    eqt.profiles_2d[1].phi =
        Interpolations.cubic_spline_interpolation(
            to_range(eqt.profiles_1d.psi) * psi_sign,
            eqt.profiles_1d.phi;
            extrapolation_bc=Interpolations.Line()
        ).(eqt.profiles_2d[1].psi * psi_sign)

    # rho 2D in meters
    RHO = sqrt.(abs.(eqt.profiles_2d[1].phi ./ (π * B0)))

    # gm2: <∇ρ²/R²>
    if false
        RHO_interpolant = Interpolations.cubic_spline_interpolation((r, z), RHO)
        for k in 1:length(eqt.profiles_1d.psi)
            tmp = [Interpolations.gradient(RHO_interpolant, PR[k][j], PZ[k][j]) for j in 1:length(PR[k])]
            dPHI2 = [j[1] .^ 2.0 .+ j[2] .^ 2.0 for j in tmp]
            eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k] .^ 2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
        end
    else
        dRHOdR, dRHOdZ = gradient(collect(r), collect(z), RHO)
        dPHI2_interpolant = Interpolations.cubic_spline_interpolation((r, z), dRHOdR .^ 2.0 .+ dRHOdZ .^ 2.0)
        for k in 1:length(eqt.profiles_1d.psi)
            dPHI2 = dPHI2_interpolant.(PR[k], PZ[k])
            eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k] .^ 2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
        end
    end
    eqt.profiles_1d.gm2[1] =
        Interpolations.cubic_spline_interpolation(
            to_range(eqt.profiles_1d.psi[2:end]) * psi_sign,
            eqt.profiles_1d.gm2[2:end];
            extrapolation_bc=Interpolations.Line()
        ).(eqt.profiles_1d.psi[1] * psi_sign)

    # ip
    eqt.global_quantities.ip = IMAS.integrate(eqt.profiles_1d.area, eqt.profiles_1d.j_tor)
    #eqt.global_quantities.ip = -gradient(eqt.profiles_1d.psi, eqt.profiles_1d.phi)[end] .* eqt.profiles_1d.gm2[end] .* eqt.profiles_1d.dvolume_dpsi[end] / (2π * constants.μ_0) * 4π

    # Geometric major and minor radii
    Rgeo = (eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]) / 2.0
    a = (eqt.profiles_1d.r_outboard[end] - eqt.profiles_1d.r_inboard[end]) / 2.0

    # vacuum magnetic field at the geometric center
    Btvac = B0 * R0 / Rgeo

    # average poloidal magnetic field
    Bpave = eqt.global_quantities.ip * constants.μ_0 / eqt.global_quantities.length_pol

    # li
    Bp2v = integrate(eqt.profiles_1d.psi, BPL)
    eqt.global_quantities.li_3 = 2.0 * Bp2v / Rgeo / (eqt.global_quantities.ip * constants.μ_0)^2

    # beta_tor
    avg_press = volume_integrate(eqt, eqt.profiles_1d.pressure) / eqt.profiles_1d.volume[end]
    eqt.global_quantities.beta_tor = abs(avg_press / (Btvac^2 / 2.0 / constants.μ_0))

    # beta_pol
    eqt.global_quantities.beta_pol = abs(avg_press / (Bpave^2 / 2.0 / constants.μ_0))

    # beta_normal
    ip = eqt.global_quantities.ip / 1e6
    eqt.global_quantities.beta_normal = eqt.global_quantities.beta_tor / abs(ip / a / Btvac) * 100

    # find quantities on separatrix
    find_x_point!(eqt)
    find_strike_points!(eqt)

    return eqt
end

"""
    flux_surface(eqt::equilibrium__time_slice, psi_level::Real)

Returns r,z coordiates of closed flux surface at given psi_level
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real)
    return flux_surface(eqt, psi_level, true)
end

"""
    flux_surface(eqt::equilibrium__time_slice, psi_level::Real, closed::Union{Nothing,Bool})

Returns r,z coordiates of open or closed flux surface at given psi_level

The `closed` parameter:

  - nothing: return all contours
  - true: all closed flux-surface that encircle the magnetic axis
  - false: all open flux-surfaces
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real, closed::Union{Nothing,Bool})
    dim1 = eqt.profiles_2d[1].grid.dim1
    dim2 = eqt.profiles_2d[1].grid.dim2
    PSI = eqt.profiles_2d[1].psi
    psi = eqt.profiles_1d.psi
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z
    return flux_surface(dim1, dim2, PSI, psi, R0, Z0, psi_level, closed)
end

function flux_surface(
    dim1::Union{AbstractVector{T},AbstractRange{T}},
    dim2::Union{AbstractVector{T},AbstractRange{T}},
    PSI::AbstractArray{T},
    psi::Union{AbstractVector{T},AbstractRange{T}},
    R0::T,
    Z0::T,
    psi_level::T,
    closed::Union{Nothing,Bool}) where {T<:Real}

    if psi_level == psi[1]
        # handle on axis value as the first flux surface
        psi_level = psi[2]

    elseif psi_level == psi[end]
        # handle boundary by finding accurate lcfs psi
        psi__boundary_level = find_psi_boundary(dim1, dim2, PSI, psi, R0, Z0; raise_error_on_not_open=false, raise_error_on_not_closed=false)
        if psi__boundary_level !== nothing
            if abs(psi__boundary_level - psi_level) < abs(psi[end] - psi[end-1])
                psi_level = psi__boundary_level
            end
        end
    end

    # contouring routine
    cl = Contour.contour(dim1, dim2, PSI, psi_level)

    prpz = Tuple{Vector{T},Vector{T}}[]
    if closed === nothing
        # if no open/closed check, then return all contours
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            R0 = 0.5 * (maximum(pr) + minimum(pr))
            Z0 = 0.5 * (maximum(pz) + minimum(pz))
            reorder_flux_surface!(pr, pz, R0, Z0; force_close=false)
            push!(prpz, (pr, pz))
        end
        return prpz

    elseif closed
        # look for closed flux-surface
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surface that close and contain magnetic axis
            if (pr[1] == pr[end]) && (pz[1] == pz[end]) && (PolygonOps.inpolygon((R0, Z0), collect(zip(pr, pz))) == 1)
                R0 = 0.5 * (maximum(pr) + minimum(pr))
                Z0 = 0.5 * (maximum(pz) + minimum(pz))
                reorder_flux_surface!(pr, pz, R0, Z0; force_close=false)
                return pr, pz, psi_level
            end
        end
        return [], [], psi_level

    elseif !closed
        # look for open flux-surfaces
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surfaces that do not close
            if (pr[1] != pr[end]) || (pz[1] != pz[end])
                R0 = 0.5 * (maximum(pr) + minimum(pr))
                Z0 = 0.5 * (maximum(pz) + minimum(pz))
                reorder_flux_surface!(pr, pz, R0, Z0; force_close=false)
                push!(prpz, (pr, pz))
            end
        end
        return prpz
    end
end

"""
    tweak_psi_to_match_psilcfs!(eqt::IMAS.equilibrium.time_slice{D}; ψbound::Union{Nothing,D}=nothing) where {D<:Real}

Tweak `eqt.profiles_2d[1].psi` and `eqt.profiles_1d.psi` so that contouring finds LCFS at `eqt.profiles_1d.psi[end]`.

If `ψbound !== nothing` this routine also shifts `eqt.profiles_2d[1].psi` and `eqt.profiles_1d.psi` so that `eqt.profiles_1d.psi[end] == ψbound`.
"""
function tweak_psi_to_match_psilcfs!(eqt::IMAS.equilibrium__time_slice{D}; ψbound::Union{Nothing,D}=nothing) where {D<:Real}
    eq1d = eqt.profiles_1d
    eq2d = eqt.profiles_2d[1]

    psia = eq1d.psi[1]
    psib = eq1d.psi[end]
    if ψbound === nothing
        delta_psib = 0.0
    else
        delta_psib = ψbound - psib
    end

    # retrace the last closed flux surface
    true_psib = IMAS.find_psi_boundary(eqt)
    if true_psib !== nothing
        # scale psirz so to match original psi bounds (also add delta_psib to get desired ψbound)
        @. eq2d.psi = (eq2d.psi - psia) * (psib - psia) / (true_psib - psia) + psia + delta_psib
        @. eq1d.psi = eq1d.psi + delta_psib
    end

    return nothing
end

function flxAvg(input::AbstractVector{T}, ll::AbstractVector{T}, fluxexpansion::AbstractVector{T}, int_fluxexpansion_dl::T)::T where {T<:Real}
    return integrate(ll, input .* fluxexpansion) / int_fluxexpansion_dl
end

"""
    volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Integrate quantity over volume
"""
function volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::T where {T<:Real}
    return integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what)
end

"""
    volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Cumulative integrate quantity over volume
"""
function cumlul_volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    return cumul_integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what)
end

"""
    surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Integrate quantity over surface
"""
function surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::T where {T<:Real}
    return integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what .* eqt.profiles_1d.gm9) ./ 2π
end

"""
    surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Cumulative integrate quantity over surface
"""
function cumlul_surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    return cumul_integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what .* eqt.profiles_1d.gm9) ./ 2π
end

"""
    find_x_point!(eqt::IMAS.equilibrium__time_slice)::eqt.boundary.x_point

Find X-points on the last closed flux surface
"""
function find_x_point!(eqt::IMAS.equilibrium__time_slice)::IDSvector{<:IMAS.equilibrium__time_slice___boundary__x_point}
    rlcfs, zlcfs = flux_surface(eqt, eqt.profiles_1d.psi[end], true)
    ll = sqrt((maximum(zlcfs) - minimum(zlcfs)) * (maximum(rlcfs) - minimum(rlcfs))) / 5.0
    private = flux_surface(eqt, eqt.profiles_1d.psi[end], false)
    Z0 = sum(zlcfs) / length(zlcfs)
    empty!(eqt.boundary.x_point)
    for (pr, pz) in private
        if sign(pz[1] - Z0) != sign(pz[end] - Z0)
            # open flux surface does not encicle the plasma
            continue
        elseif minimum_distance_two_shapes(pr, pz, rlcfs, zlcfs) > ll
            # secondary Xpoint far away
            continue
        elseif (sum(pz) < Z0)
            index = argmax(pz)
        elseif (sum(pz) > Z0)
            # upper private region
            index = argmin(pz)
        else
            continue
        end
        indexcfs = argmin((rlcfs .- pr[index]) .^ 2 .+ (zlcfs .- pz[index]) .^ 2)
        resize!(eqt.boundary.x_point, length(eqt.boundary.x_point) + 1)
        eqt.boundary.x_point[end].r = (pr[index] + rlcfs[indexcfs]) / 2.0
        eqt.boundary.x_point[end].z = (pz[index] + zlcfs[indexcfs]) / 2.0
    end

    r, z, PSI_interpolant = ψ_interpolant(eqt.profiles_2d[1])
    # refine x-point location
    for rz in eqt.boundary.x_point
        res = Optim.optimize(
            x -> Bp(PSI_interpolant, [rz.r + x[1]], [rz.z + x[2]])[1],
            [0.0, 0.0],
            Optim.NelderMead(),
            Optim.Options(; g_tol=1E-8)
        )
        rz.r += res.minimizer[1]
        rz.z += res.minimizer[2]
    end

    return eqt.boundary.x_point
end

"""
    find_2nd_x_point!(eqt::IMAS.equilibrium__time_slice)

Returns a tuple with coordiantes of upper X-point, which could be on the lcfs or not.
"""
function find_2nd_x_point!(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation)
    psi_separatrix = find_psi_boundary(eqt; raise_error_on_not_open=true) # psi at LCFS
      
    #first: look for second x point in equilibrium time slice
    xpoints = eqt.boundary.x_point  # list of x points from eqt
    z_xpoints = zeros(length(xpoints))
    psi_xpoints = Float64[]
    for jj in 1:length(xpoints)  # retrive Z coordinates of all x points
        z_xpoints[jj] = xpoints[jj].z
        push!(psi_xpoints, PSI_interpolant(xpoints[jj].r,xpoints[jj].z)[1]) # save psi of nulls 
    end

    c = sign(z_xpoints[argmin(abs.(psi_xpoints.-psi_separatrix))]) # find sign of Z coordinate of first null
    index = 1:length(xpoints)
    index = index[-c.*z_xpoints.>0] # index in xpoints of nulls with positive Z (if upper)
    
    if isempty(index)
        # no nulls found with Z coordinate with opposite sign of first null
        # find second null looking for Bp = 0
        ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis
        # optimization to find x-point location
        null2 = [xpoints[argmin(abs.(psi_xpoints.-psi_separatrix))].r, xpoints[argmin(abs.(psi_xpoints.-psi_separatrix))].z]
        null2[2] = -1*null2[2]+ZA # start optimization from flipped first xpoint
         res = Optim.optimize(
            x -> IMAS.Bp(PSI_interpolant, [null2[1] + x[1]], [null2[2] + x[2]])[1],
            [0.0, 0.0],
            Optim.NelderMead(),
            Optim.Options(; g_tol=1E-8)
        )
        null2[1] += res.minimizer[1]
        null2[2] += res.minimizer[2]
        return (null2[1],null2[2])

    else
        if sum(-c.*z_xpoints.>0) == 1
            # we have only one null having Z with opposite sign of first x point
            return (xpoints[index[1]].r,xpoints[index[1]].z)
        else
            # we have more than one null  having Z with opposite sign of first x point
            # choose xpoint with closest psi to separatrix
            psi = Float64[]
            for i in index
                push!(psi, PSI_interpolant(xpoints[i].r,xpoints[i].z)[1]) # find psi of nulls with Z>0 (if upper)
            end

            index = index[argmin(abs.(psi.-psi_separatrix))] # find index of null with psi closest to psi_separatrix

            return (xpoints[index].r,xpoints[index].z)
        end
    end 
end

function find_2nd_x_point!(dd::IMAS.dd)
    r, z, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d[])
    return find_2nd_x_point!(dd.equilibrium.time_slice[], PSI_interpolant)
end 

"""
    miller_R_a_κ_δ_ζ(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)

Returns R0, a, κ, δu, δl, ζou, ζol, ζil, ζiu of a contour
"""
function miller_R_a_κ_δ_ζ(pr::Vector{T}, pz::Vector{T}, r_at_max_z::T, max_z::T, r_at_min_z::T, min_z::T, z_at_max_r::T, max_r::T, z_at_min_r::T, min_r::T) where {T<:Real}
    R0 = 0.5 * (max_r + min_r)
    a = 0.5 * (max_r - min_r)
    b = 0.5 * (max_z - min_z)
    κ = b / a
    δu = (R0 - r_at_max_z) / a
    δl = (R0 - r_at_min_z) / a
    ζou, ζol, ζil, ζiu = luce_squareness(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
    return R0, a, κ, δu, δl, ζou, ζol, ζil, ζiu
end

function miller_R_a_κ_δ_ζ(pr::Vector{T}, pz::Vector{T}) where {T<:Real}
    (imaxr, iminr, imaxz, iminz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = fluxsurface_extrema(pr, pz)
    return miller_R_a_κ_δ_ζ(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
end

"""
    fluxsurface_extrema(pr::Vector{T}, pz::Vector{T}) where {T<:Real}

Returns extrema indexes and values of R,Z flux surfaces vectors:

    imaxr, iminr,
    imaxz, iminz,
    r_at_max_z, max_z,
    r_at_min_z, min_z,
    z_at_max_r, max_r,
    z_at_min_r, min_r
"""
function fluxsurface_extrema(pr::Vector{T}, pz::Vector{T}) where {T<:Real}
    _, imaxr = findmax(pr)
    _, iminr = findmin(pr)
    _, imaxz = findmax(pz)
    _, iminz = findmin(pz)
    r_at_max_z, max_z = pr[imaxz], pz[imaxz]
    r_at_min_z, min_z = pr[iminz], pz[iminz]
    z_at_max_r, max_r = pz[imaxr], pr[imaxr]
    z_at_min_r, min_r = pz[iminr], pr[iminr]
    return (imaxr, iminr, imaxz, iminz,
        r_at_max_z, max_z,
        r_at_min_z, min_z,
        z_at_max_r, max_r,
        z_at_min_r, min_r)
end

"""
    luce_squareness(pr::AbstractVector{T}, pz::AbstractVector{T}, r_at_max_z::T, max_z::T, r_at_min_z::T, min_z::T, z_at_max_r::T, max_r::T, z_at_min_r::T, min_r::T) where {T<:Real}

Squareness from: "An analytic functional form for characterization and generation of axisymmetric plasma boundaries"
T.C. Luce, Plasma Phys. Control. Fusion 55 (2013) http://dx.doi.org/10.1088/0741-3335/55/9/095009

Returns: zetaou, zetaol, zetail, zetaiu
"""
function luce_squareness(
    pr::AbstractVector{T}, pz::AbstractVector{T},
    r_at_max_z::T, max_z::T,
    r_at_min_z::T, min_z::T,
    z_at_max_r::T, max_r::T,
    z_at_min_r::T, min_r::T) where {T<:Real}

    # zetaou
    PO = (r_at_max_z, z_at_max_r)
    PE = (max_r, max_z)
    PD = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz)[2][1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaou = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetaol
    PO = (r_at_min_z, z_at_max_r)
    PE = (max_r, min_z)
    PD = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz)[2][1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaol = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetaiu
    PO = (r_at_max_z, z_at_min_r)
    PE = (min_r, max_z)
    PD = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz)[2][1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaiu = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetail
    PO = (r_at_min_z, z_at_min_r)
    PE = (min_r, min_z)
    PD = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz)[2][1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetail = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    return zetaou, zetaol, zetail, zetaiu
end


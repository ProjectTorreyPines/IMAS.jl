document[Symbol("Physics transport")] = Symbol[]

"""
    profile_from_z_transport(
        profile_old::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        transport_grid::AbstractVector{<:Real},
        z_transport_grid::AbstractVector{<:Real},
        rho_ped::Real=0.0)

Updates profile_old with the scale lengths given by z_transport_grid

If `rho_ped > transport_grid[end]` then scale-length is linearly interpolated between `transport_grid[end]` and `rho_ped`

if `rho_ped < transport_grid[end]` then scale-length then boundary condition is at `transport_grid[end]`
"""
function profile_from_z_transport(
    profile_old::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    transport_grid::AbstractVector{<:Real},
    z_transport_grid::AbstractVector{<:Real},
    rho_ped::Real=0.0)

    transport_indices = [argmin_abs(rho, rho_x) for rho_x in transport_grid]
    index_ped = argmin_abs(rho, rho_ped)
    index_last = transport_indices[end]
    if index_ped > index_last
        z_old = calc_z(rho, profile_old, :backward)
        transport_indices = vcat(1, transport_indices, index_ped)
        z_transport_grid = vcat(0.0, z_transport_grid, z_old[index_ped])
    else
        transport_indices = vcat(1, transport_indices)
        z_transport_grid = vcat(0.0, z_transport_grid)
    end

    z = interp1d(transport_indices, z_transport_grid).(1:index_last)

    profile_new = similar(profile_old)
    profile_new[index_last:end] .= @views profile_old[index_last:end]
    profile_new[1:index_last] .= @views integ_z(rho[1:index_last], -z, profile_new[index_last])

    return profile_new
end

@compat public profile_from_z_transport
push!(document[Symbol("Physics transport")], :profile_from_z_transport)

"""
    profile_from_rotation_shear_transport(
        profile_old::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        transport_grid::AbstractVector{<:Real},
        rotation_shear_transport_grid::AbstractVector{<:Real},
        rho_ped::Real=0.0)

Updates rotation profile using the rotation shear given by rotation_shear_transport_grid

We integrate dω/dr to get ω(r).

If `rho_ped > transport_grid[end]` then rotation shear is linearly interpolated between `transport_grid[end]` and `rho_ped`

If `rho_ped < transport_grid[end]` then boundary condition is at `transport_grid[end]`
"""
function profile_from_rotation_shear_transport(
    profile_old::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    transport_grid::AbstractVector{<:Real},
    rotation_shear_transport_grid::AbstractVector{<:Real},
    rho_ped::Real=0.0)

    transport_indices = [IMAS.argmin_abs(rho, rho_x) for rho_x in transport_grid]
    index_ped = IMAS.argmin_abs(rho, rho_ped)
    index_last = transport_indices[end]

    if index_ped > index_last
        # If rho_ped is beyond transport_grid[end], extend interpolation
        dw_dr_old = gradient(rho, profile_old; method=:backward)
        transport_indices = vcat(1, transport_indices, index_ped)
        rotation_shear_transport_grid = vcat(0.0, rotation_shear_transport_grid, dw_dr_old[index_ped])
    else
        # If rho_ped is within transport_grid, use transport_grid[end] as boundary
        transport_indices = vcat(1, transport_indices)
        rotation_shear_transport_grid = vcat(0.0, rotation_shear_transport_grid)
    end

    # Interpolate rotation shear to the integration range
    dw_dr = IMAS.interp1d(transport_indices, rotation_shear_transport_grid).(1:index_last)

    # Set up the profile
    profile_new = similar(profile_old)
    profile_new[index_last:end] .= @views profile_old[index_last:end]

    # Integrate rotation shear from boundary inward to axis
    # Integration: ω(rho) = ω(rho_boundary) - ∫_{rho}^{rho_boundary} (dω/dr') dr'
    profile_new[index_last] = profile_old[index_last]
    for i in (index_last-1):-1:1
        # Trapezoidal integration: negative because we integrate inward
        dw_avg = (dw_dr[i] + dw_dr[i+1]) / 2.0
        dr = rho[i+1] - rho[i]
        profile_new[i] = profile_new[i+1] - dw_avg * dr
    end

    return profile_new
end

@compat public profile_from_rotation_shear_transport
push!(document[Symbol("Physics transport")], :profile_from_rotation_shear_transport)

"""
    total_fluxes(
        core_transport::IMAS.core_transport{T},
        cp1d::IMAS.core_profiles__profiles_1d,
        rho_total_fluxes::AbstractVector{<:Real};
        time0::Float64) where {T<:Real}

Sums up all the fluxes and returns it as a core_transport.model IDS
"""
function total_fluxes(
    core_transport::IMAS.core_transport{T},
    cp1d::IMAS.core_profiles__profiles_1d,
    rho_total_fluxes::AbstractVector{<:Real};
    time0::Float64) where {T<:Real}
    total_flux1d = core_transport__model___profiles_1d{T}()
    return total_fluxes!(total_flux1d, core_transport, cp1d, rho_total_fluxes; time0)
end

function total_fluxes!(
    total_flux1d::IMAS.core_transport__model___profiles_1d{T},
    core_transport::IMAS.core_transport{T},
    cp1d::IMAS.core_profiles__profiles_1d,
    rho_total_fluxes::AbstractVector{<:Real};
    time0::Float64) where {T<:Real}

    total_flux1d.grid_flux.rho_tor_norm = rho_total_fluxes

    skip_flux_list = [:unknown, :unspecified, :combined]
    index_to_name = index_2_name(core_transport.model)

    # initialize ions
    for ion in cp1d.ion
        resize!(total_flux1d.ion, "element[1].a" => ion.element[1].a, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label; wipe=false)
    end
    for model in core_transport.model
        if !isempty(model.profiles_1d) && time0 >= model.profiles_1d[1].time
            model1d = model.profiles_1d[time0]
            for ion in model1d.ion
                resize!(total_flux1d.ion, "element[1].a" => ion.element[1].a, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label; wipe=false)
            end
        end
    end

    # defines paths to fill
    paths = Union{Tuple{Symbol,Symbol},Tuple{Symbol,Int64,Symbol},Tuple{Symbol}}[]
    push!(paths, (:electrons, :energy))
    push!(paths, (:electrons, :particles))
    for k in eachindex(total_flux1d.ion)
        push!(paths, (:ion, k, :energy))
        push!(paths, (:ion, k, :particles))
    end
    push!(paths, (:momentum_tor,))
    push!(paths, (:total_ion_energy,))

    # zero out total_flux
    for path in paths
        ids1 = goto(total_flux1d, path)
        if hasdata(ids1, :flux) && length(getfield(ids1, :flux)) == size(rho_total_fluxes)
            fill!(getproperty(ids1, :flux), zero(T))
        else
            setproperty!(ids1, :flux, zeros(T, size(rho_total_fluxes)))
        end
    end

    if !isempty(rho_total_fluxes)
        for model in core_transport.model
            # Make sure we don't double count a specific flux type
            if index_to_name[model.identifier.index] ∈ skip_flux_list
                @debug "skipped model.identifier.index = $(model.identifier.index) [$(index_to_name[model.identifier.index])]"
                continue
            end
            push!(skip_flux_list, index_to_name[model.identifier.index])
            # ignore models that start after time0
            if nearest_causal_time(model.profiles_1d, time0; bounds_error=false).out_of_bounds
                continue
            end

            m1d = model.profiles_1d[time0]
            x = m1d.grid_flux.rho_tor_norm
            x_1 = argmin_abs(rho_total_fluxes, x[1])
            x_2 = argmin_abs(rho_total_fluxes, x[end])

            for path in paths
                ids1 = try
                    goto(m1d, path)
                catch e
                    if isa(e, InterruptException)
                        rethrow(e)
                    else
                        continue
                    end
                end
                if !hasdata(ids1, :flux)
                    continue
                end
                ids2 = goto(total_flux1d, path)
                y = getproperty(ids1, :flux)
                if rho_total_fluxes == x
                    getproperty(ids2, :flux)[x_1:x_2] .+= y
                else
                    getproperty(ids2, :flux)[x_1:x_2] .+= interp1d(x, y, :linear).(rho_total_fluxes[x_1:x_2])
                end
            end
        end
    end

    return total_flux1d
end

"""
    total_fluxes(dd::IMAS.DD, rho_total_fluxes::AbstractVector{<:Real}=dd.core_profiles.profiles_1d[].grid.rho_tor_norm; time0::Float64=dd.global_time)
"""
function total_fluxes(dd::IMAS.DD, rho_total_fluxes::AbstractVector{<:Real}=dd.core_profiles.profiles_1d[].grid.rho_tor_norm; time0::Float64=dd.global_time)
    return total_fluxes(dd.core_transport, dd.core_profiles.profiles_1d[], rho_total_fluxes; time0)
end


"""
    _diffusivity_terms(dd::IMAS.DD{T}; ne::Vector{T}=T[], Te::Vector{T}=T[], Ti::Vector{T}=T[], ω::Vector{T}=T[]) where {T<:Real}

Geometry, profile, and conducted-power terms shared by `calculate_diffusivities` and `pedestal_diffusivities`.

Everything is returned on the `total_sources` grid `rho` (= `rho_tor_norm`). The logarithmic gradients
`dln*dr` are taken with respect to `rho` and are positive where a profile decreases outwards.
"""
function _diffusivity_terms(dd::IMAS.DD{T}; ne::Vector{T}=T[], Te::Vector{T}=T[], Ti::Vector{T}=T[], ω::Vector{T}=T[]) where {T<:Real}

    cs1d = IMAS.total_sources(dd)
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d
    rho_cp = cs1d.grid.rho_tor_norm
    rho_eq = eqt1d.rho_tor_norm

    # Volumes and surfaces
    surf = IMAS.interp1d(rho_eq, eqt1d.surface).(rho_cp)
    dVdρ = IMAS.interp1d(rho_eq, eqt1d.dvolume_drho_tor).(rho_cp)

    # used below into real-space gradients. Uses <|∇ρ|²> ≈ <|∇ρ|>² because `gm3` is not among the
    Leff = @. eqt1d.rho_tor[end] * dVdρ / surf
    isfinite(Leff[1]) || (Leff[1] = Leff[2])             # surf and dVdρ both vanish on axis

    # Default profiles if not provided
    ne = isempty(ne) ? cp1d.electrons.density_thermal : ne
    Te = isempty(Te) ? cp1d.electrons.temperature : Te
    Ti = isempty(Ti) ? cp1d.ion[1].temperature : Ti
    ω = isempty(ω) ? cp1d.rotation_frequency_tor_sonic : ω

    zeff = cp1d.zeff
    ni = @. ne / zeff                                    # main-ion density

    # Logarithmic gradients
    dlntedr = .-IMAS.calc_z(rho_cp, Te, :backward)
    dlnnedr = .-IMAS.calc_z(rho_cp, ne, :backward)
    dlntidr = .-IMAS.calc_z(rho_cp, Ti, :backward)

    # Volume-integrated fluxes, and the conducted part of the power (total minus the 5/2⋅Γ⋅T convected part)
    Γe = cs1d.electrons.particles_inside                 # [particles/s]
    Γi = @. Γe / zeff                                    # main-ion flux, consistent with `ni`
    Qe_cond = @. cs1d.electrons.power_inside - 2.5 * IMAS.mks.e * Te * Γe
    Qi_cond = @. cs1d.total_ion_power_inside - 2.5 * IMAS.mks.e * Ti * Γi

    # momentum pressure
    mi = cp1d.ion[1].element[1].a * IMAS.mks.m_p
    Rmaj = eqt.boundary.geometric_axis.r               # major radius
    dωdr = .-IMAS.calc_z(rho_cp, ω, :backward) .* ω
    pφ = dωdr .* ne .* mi .* Rmaj .^ 2

    return (; rho=rho_cp, surf, Leff, ne, ni, Te, Ti, zeff, dlnnedr, dlntedr, dlntidr,
        Γe, Qe_cond, Qi_cond, pφ, torque=cs1d.torque_tor_inside)
end

"""
    calculate_diffusivities(dd::IMAS.DD{T}; ne::Vector{T}=T[], Te::Vector{T}=T[], Ti::Vector{T}=T[], ω::Vector{T}=T[]) where {T<:Real}

Compute transport particle, heat, and rotation effective diffusivities `(Dn, χe, χi, χφ)` [m²/s] from desired profiles.

The coefficients are inverted from the flux-surface-averaged transport equation

    S_inside(ρ) = -Y ⟨|∇ρ_tor|²⟩ (dV/dρ_tor) dX/dρ_tor

where `S_inside` is a volume-integrated source, `X` the driving profile and `Y` the diffusivity. Gradients
are evaluated against `rho_tor_norm` and converted to real-space gradients with the effective radial length
`Leff = 1/⟨|∇ρ_tor_norm|⟩ = ρ_tor_boundary (dV/dρ_tor) / surface` [m], which relies on ⟨|∇ρ|²⟩ ≈ ⟨|∇ρ|⟩².

`χe` and `χi` follow the conduction-only convention `χ = -q_cond / (n dT/dr)` used by edge codes: the
convected `5/2 Γ T` part of the power is subtracted from the numerator, and the denominator holds the
temperature gradient rather than the pressure gradient.
"""
function calculate_diffusivities(dd::IMAS.DD{T}; ne::Vector{T}=T[], Te::Vector{T}=T[], Ti::Vector{T}=T[], ω::Vector{T}=T[]) where {T<:Real}

    trm = _diffusivity_terms(dd; ne, Te, Ti, ω)

    # Diffusivities
    Dn = @. trm.Leff * trm.Γe / (trm.surf * trm.dlnnedr * trm.ne)
    χe = @. trm.Leff * trm.Qe_cond / (trm.surf * trm.ne * IMAS.mks.e * trm.Te * trm.dlntedr)
    χi = @. trm.Leff * trm.Qi_cond / (trm.surf * trm.ni * IMAS.mks.e * trm.Ti * trm.dlntidr)
    χφ = @. trm.Leff * trm.torque / (trm.surf * trm.pφ)

    # Patch singular boundary values
    Dn[1] = Dn[2]
    χe[1] = χe[2]
    χi[1] = χi[2]
    χφ[1] = χφ[2]

    return Dn, χe, χi, χφ
end

@compat public calculate_diffusivities
push!(document[Symbol("Physics transport")], :calculate_diffusivities)

"""
    _pedestal_average(rho::AbstractVector{<:Real}, y::AbstractVector{<:Real}, index::AbstractVector{Int})

Average `y` over the samples `index`, dropping non-finite and non-positive entries.

Returns `NaN` when nothing usable survives, leaving it to the caller to decide on a fallback.
"""
function _pedestal_average(rho::AbstractVector{<:Real}, y::AbstractVector{<:Real}, index::AbstractVector{Int})
    keep = [k for k in index if isfinite(y[k]) && y[k] > 0.0]
    if isempty(keep)
        return NaN
    elseif length(keep) == 1
        return float(y[keep[1]])
    end
    x = @views rho[keep]
    return trapz(x, @views y[keep]) / (x[end] - x[1])
end

"""
    pedestal_diffusivities(dd::IMAS.DD{T}; rho_ped::Float64=NaN, ne::Vector{T}=T[], Te::Vector{T}=T[], Ti::Vector{T}=T[]) where {T<:Real}

Cross-field transport coefficients [m²/s] averaged over the pedestal region `[rho_ped, 1.0]`, for use as
boundary/edge transport inputs (e.g. the `D_perp`/`chi_perp` of a SOL or edge-plasma code).

Returns `(; D_perp, chi_e, chi_i, chi_perp, rho_ped)`.

`chi_perp` is a *single* heat diffusivity for both species, built to conserve the total conducted heat flux
(`chi_perp = (Qe_cond + Qi_cond) / (surface ⟨|∇ρ|⟩ e (ne dTe/dρ + ni dTi/dρ))`) rather than by naively
averaging `chi_e` and `chi_i` — this matches edge codes that run a common `kye = kyi`.

`rho_ped` defaults to `dd.summary.local.pedestal.position.rho_tor_norm`, falling back to a
[`pedestal_finder`](@ref) fit of the electron density when the summary is not filled.

Any coefficient for which no finite, positive sample exists in the pedestal comes back as `NaN`.
"""
function pedestal_diffusivities(dd::IMAS.DD{T}; rho_ped::Float64=NaN, ne::Vector{T}=T[], Te::Vector{T}=T[], Ti::Vector{T}=T[]) where {T<:Real}

    trm = _diffusivity_terms(dd; ne, Te, Ti)

    if isnan(rho_ped)
        if hasdata(dd.summary.local.pedestal.position, :rho_tor_norm)
            rho_ped = @ddtime(dd.summary.local.pedestal.position.rho_tor_norm)
        else
            cp1d = dd.core_profiles.profiles_1d[]
            rho_ped = 1.0 - pedestal_finder(cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm).width
        end
    end

    # pedestal samples, always keeping at least the two outermost points
    index = findall(x -> x >= rho_ped, trm.rho)
    if length(index) < 2
        index = collect(length(trm.rho)-1:length(trm.rho))
    end

    Dn = @. trm.Leff * trm.Γe / (trm.surf * trm.dlnnedr * trm.ne)
    χe = @. trm.Leff * trm.Qe_cond / (trm.surf * trm.ne * IMAS.mks.e * trm.Te * trm.dlntedr)
    χi = @. trm.Leff * trm.Qi_cond / (trm.surf * trm.ni * IMAS.mks.e * trm.Ti * trm.dlntidr)

    # single heat diffusivity conserving the total conducted heat flux
    χ = @. trm.Leff * (trm.Qe_cond + trm.Qi_cond) /
           (trm.surf * IMAS.mks.e * (trm.ne * trm.Te * trm.dlntedr + trm.ni * trm.Ti * trm.dlntidr))

    return (;
        D_perp=_pedestal_average(trm.rho, Dn, index),
        chi_e=_pedestal_average(trm.rho, χe, index),
        chi_i=_pedestal_average(trm.rho, χi, index),
        chi_perp=_pedestal_average(trm.rho, χ, index),
        rho_ped)
end

@compat public pedestal_diffusivities
push!(document[Symbol("Physics transport")], :pedestal_diffusivities)

@compat public total_fluxes
push!(document[Symbol("Physics transport")], :total_fluxes)

"""
    rho_q_surface(rho::AbstractVector{<:Real}, q::AbstractVector{<:Real}, q_res::Real)

Returns the outermost rho where abs(q) crosses q_res (linearly interpolated), or `nothing` if no crossing is found
"""
function rho_q_surface(rho::AbstractVector{<:Real}, q::AbstractVector{<:Real}, q_res::Real)
    if abs(q[end]) == q_res
        return rho[end]
    end
    for k in length(q)-1:-1:1
        q1 = abs(q[k]) - q_res
        q2 = abs(q[k+1]) - q_res
        if q1 == 0.0
            return rho[k]
        elseif q1 * q2 < 0.0
            return rho[k] + (rho[k+1] - rho[k]) * q1 / (q1 - q2)
        end
    end
    return nothing
end

@compat public rho_q_surface
push!(document[Symbol("Physics transport")], :rho_q_surface)

"""
    island_window(rho::AbstractVector{T}, rho_res::Real, island_width::Real, w_c::Real) where {T<:Real}

Radial shape [0, 1] of the transport enhancement due to a magnetic island of full width `island_width` centered at `rho_res`

    W(ρ) = f₀(w) * ½ [tanh((ρ - ρ_res + w/2) / δ) - tanh((ρ - ρ_res - w/2) / δ)]

  - `δ = w_c / 2`: the edges are smoothed over the critical island width `w_c`, over which perpendicular and parallel transport compete
  - `f₀(w) = w² / (w² + w_c²)`: islands smaller than `w_c` do not flatten the profiles (heuristic)
"""
function island_window(rho::AbstractVector{T}, rho_res::Real, island_width::Real, w_c::Real) where {T<:Real}
    @assert w_c > 0.0 "w_c must be > 0.0"
    if island_width <= 0.0
        return zeros(T, length(rho))
    end
    w = island_width
    δ = w_c / 2
    f0 = w^2 / (w^2 + w_c^2)
    return @. f0 * 0.5 * (tanh((rho - rho_res + w / 2) / δ) - tanh((rho - rho_res - w / 2) / δ))
end

@compat public island_window
push!(document[Symbol("Physics transport")], :island_window)

"""
    tearing_mode_transport!(
        m1d::IMAS.core_transport__model___profiles_1d,
        cp1d::IMAS.core_profiles__profiles_1d,
        rho_transport::AbstractVector{<:Real},
        window::AbstractVector{<:Real};
        chi_e::Real,
        chi_i::Real,
        D::Real,
        chi_tor::Real,
        R0::Real)

Fills `m1d` with the fluxes at `rho_transport` due to ad hoc diffusivities `χ(ρ) = χ_island W(ρ)` inside magnetic islands,
where `W(ρ)` is the radial shape of the islands evaluated on the `cp1d` grid (see `island_window`).

This mimics the short-circuiting of radial transport by parallel transport along the island field lines:
when `χ_island` is much larger than the turbulent diffusivity, flux matching drives the gradients inside the island
below the turbulence threshold and the profiles flatten (the bootstrap hole then follows from the flattened pressure).

The fluxes follow the same convention as `calculate_diffusivities`, with `Leff = (dV/dρ) / surface`:

    q_e = χ_e e n_e T_e (-dlnT_e/dρ) / Leff        [W/m²]
    q_i = Σ_k χ_i e n_k T_k (-dlnT_k/dρ) / Leff
    Γ_s = D n_s (-dln n_s/dρ) / Leff                [1/m²/s]
    Π   = χ_tor m_i n_e R0² (-dω/dρ) / Leff          [N/m]

with `ρ = rho_tor_norm`, temperatures in eV, and the diffusivities in m²/s.
"""
function tearing_mode_transport!(
    m1d::IMAS.core_transport__model___profiles_1d,
    cp1d::IMAS.core_profiles__profiles_1d,
    rho_transport::AbstractVector{<:Real},
    window::AbstractVector{<:Real};
    chi_e::Real,
    chi_i::Real,
    D::Real,
    chi_tor::Real,
    R0::Real)

    rho = cp1d.grid.rho_tor_norm
    gridpoints = [argmin_abs(rho, rho_x) for rho_x in rho_transport]

    dVdρ = gradient(rho, cp1d.grid.volume)[gridpoints]
    Leff = dVdρ ./ cp1d.grid.surface[gridpoints]
    W = window[gridpoints]

    # conductive flux of a profile X with diffusivity Y: Y * X * (-dlnX/dρ) / Leff (calc_z returns dlnX/dρ)
    flux(Y, X) = .-Y .* W .* X[gridpoints] .* calc_z(rho, X, :backward)[gridpoints] ./ Leff

    m1d.grid_flux.rho_tor_norm = rho_transport

    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature
    m1d.electrons.energy.flux = flux(chi_e, Te) .* ne[gridpoints] .* mks.e
    m1d.electrons.particles.flux = flux(D, ne)

    qi = zeros(eltype(m1d.electrons.energy.flux), length(gridpoints))
    for ion in cp1d.ion
        mion = resize!(m1d.ion, "element[1].a" => ion.element[1].a, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label)
        if !ismissing(ion, :density_thermal) && all(ion.density_thermal .> 0.0)
            mion.particles.flux = flux(D, ion.density_thermal)
            if !ismissing(ion, :temperature)
                mion.energy.flux = flux(chi_i, ion.temperature) .* ion.density_thermal[gridpoints] .* mks.e
                qi .+= mion.energy.flux
            else
                mion.energy.flux = zero(qi)
            end
        else
            mion.particles.flux = zero(qi)
            mion.energy.flux = zero(qi)
        end
    end
    m1d.total_ion_energy.flux = qi

    if !ismissing(cp1d, :rotation_frequency_tor_sonic)
        ω = cp1d.rotation_frequency_tor_sonic
        mi = cp1d.ion[1].element[1].a * mks.m_p
        dωdρ = gradient(rho, ω)[gridpoints]
        m1d.momentum_tor.flux = @. -chi_tor * W * mi * ne[gridpoints] * R0^2 * dωdρ / Leff
    else
        m1d.momentum_tor.flux = zero(qi)
    end

    return m1d
end

@compat public tearing_mode_transport!
push!(document[Symbol("Physics transport")], :tearing_mode_transport!)
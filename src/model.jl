using FUSE


#= ============== =#
# initialization #
#= ============== =#

profiles_1d = FUSE.core_profiles__profiles_1d()

κ = plasma_parameters[:κ]
St = plasma_parameters[:St]
Sn = plasma_parameters[:Sn]
Sj = plasma_parameters[:Sj]
Zeff = plasma_parameters[:Zeff]

profiles_1d.grid.rho_tor_norm = range(0.0, 1.0, length=101)
profiles_1d.j_total = (x;_...) -> (1.0 .- x.^2).^Sj
profiles_1d.electrons.density = (x;_...) -> (1.0 .- x.^2).^Sn
profiles_1d.electrons.temperature = (x;_...) -> (1.0 .- x.^2).^St

println(profiles_1d.electrons.pressure)



bootstrapCoefficient = FUSE.collisionless_bootstrap(physics_models[:bootstrapModel], κ, St, Sn, Sj, Zeff)
println(bootstrapCoefficient)
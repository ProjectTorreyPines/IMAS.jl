#using Trapz
expressions = Dict{String,Function}()

expressions["core_profiles.profiles_1d[:].electrons.pressure"] =
    (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density * 1.60218e-19

expressions["core_profiles.profiles_1d[:].electrons.density"] =
    (rho_tor_norm; electrons, _...) -> electrons.pressure ./ (electrons.temperature * 1.60218e-19)

expressions["core_profiles.profiles_1d[:].electrons.temperature"] =
    (rho_tor_norm; electrons, _...) -> electrons.pressure ./ (electrons.density * 1.60218e-19)

expressions["equilibrium.time_slice[:].global_quantities.energy_mhd"] =
    (;time_slice, _...) -> 3 / 2 * trapz(time_slice.profiles_1d.volume, time_slice.profiles_1d.pressure)
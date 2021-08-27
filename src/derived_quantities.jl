derived_quantities = Dict{String,Function}()

derived_quantities["core_profiles.profiles_1d[:].electrons.pressure"] = (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density * 1.60218e-19
derived_quantities["core_profiles.profiles_1d[:].electrons.density"] = (rho_tor_norm; electrons, _...) -> electrons.pressure ./ (electrons.temperature * 1.60218e-19)
derived_quantities["core_profiles.profiles_1d[:].electrons.temperature"] = (rho_tor_norm; electrons, _...) -> electrons.pressure ./ (electrons.density * 1.60218e-19)

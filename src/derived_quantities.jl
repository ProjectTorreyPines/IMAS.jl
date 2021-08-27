derived_quantities = Dict{String,Function}()

derived_quantities["core_profiles.profiles_1d[:].electrons.pressure"] = (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density * 1.60218e-19;

for (k, v) in derived_quantities
    k1 = replace(replace(k, "." => "__"), "[:]" => "_")
    derived_quantities[k1] = pop!(derived_quantities, k)
end
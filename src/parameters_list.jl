# begin adding parameters to simple dictionary
imas_parameters = Dict()

#= ============== =#
#  PHYSICS_MODELS  #
#= ============== =#
imas_parameters[:PHYSICS_MODELS] = Dict()

# bootstrapModel
options = Dict()
options[:pomphrey] = SwitchOption(CBSPomphrey(), "", "Pomphrey, N. Bootstrap dependence on plasma profile parameters. PPPL, 1992.")
options[:gi] = SwitchOption(CBSGi(), "", "Gi et al., Fus. Eng. Design 89 2709 (2014)")
options[:wilson] = SwitchOption(CBSWilson(), "", "Wilson et al., Nucl. Fusion 32 257 (1992)")
options[:user] = ScalarParameter(0.7, "", "User-defined constant")
imas_parameters[:PHYSICS_MODELS][:bootstrapModel] = SwitchParameter(options, :gi, "Bootstrap current efficiency model")

# this should really go under a 0D data structure
#= ================= =#
#  PLASMA_PARAMETERS  #
#= ================= =#

imas_parameters[:PLASMA_PARAMETERS] = Dict()

# profiles shapes
imas_parameters[:PLASMA_PARAMETERS][:Sn] = ScalarParameter(1.0, "", "Shape of density profile (1-x^2)^Sn")
imas_parameters[:PLASMA_PARAMETERS][:St] = ScalarParameter(1.0, "", "Shape of temperature profile (1-x^2)^St")
imas_parameters[:PLASMA_PARAMETERS][:Sj] = ScalarParameter(1.0, "", "Shape of current density profile (1-x^2)^Sj")

# profiles scales
imas_parameters[:PLASMA_PARAMETERS][:Te0] = ScalarParameter(1.0E3, "eV", "electron temperature on axis")
imas_parameters[:PLASMA_PARAMETERS][:ne0] = ScalarParameter(1E19, "m^-3", "electron density on axis")

# plasma composition
imas_parameters[:PLASMA_PARAMETERS][:Zeff] = ScalarParameter(2, "", "plasma effective charge")

# plasma geometry
imas_parameters[:PLASMA_PARAMETERS][:R0] = ScalarParameter(1.7, "m", "Major radius")
imas_parameters[:PLASMA_PARAMETERS][:aspect_ratio] = ScalarParameter(2.7, "", "plasma aspect ratio")
imas_parameters[:PLASMA_PARAMETERS][:elongation] = ScalarParameter(1.9, "", "plasma elongation")

#= ================= =#

# define ImasParameters for rapid access of .value property
for (k, v) in imas_parameters
    imas_parameters[k] = ImasParameters(v)
end

plasma_parameters = imas_parameters[:PLASMA_PARAMETERS]
physics_models = imas_parameters[:PHYSICS_MODELS]
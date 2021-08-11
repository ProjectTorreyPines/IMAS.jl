# begin adding parameters to simple dictionary
fuse_parameters = Dict()

#= ============== =#
#  PHYSICS_MODELS  #
#= ============== =#
fuse_parameters[:PHYSICS_MODELS] = Dict()

# bootstrapModel
options = Dict()
options[:pomphrey] = SwitchOption(CBSPomphrey(), "", "Pomphrey, N. Bootstrap dependence on plasma profile parameters. PPPL, 1992.")
options[:gi] = SwitchOption(CBSGi(), "", "Gi et al., Fus. Eng. Design 89 2709 (2014)")
options[:wilson] = SwitchOption(CBSWilson(), "", "Wilson et al., Nucl. Fusion 32 257 (1992)")
options[:user] = ScalarParameter(0.7, "", "User-defined constant")
fuse_parameters[:PHYSICS_MODELS][:bootstrapModel] = SwitchParameter(options, :gi, "Bootstrap current efficiency model")

# this should really go under a 0D data structure
#= ================= =#
#  PLASMA_PARAMETERS  #
#= ================= =#

fuse_parameters[:PLASMA_PARAMETERS] = Dict()

# Sn
fuse_parameters[:PLASMA_PARAMETERS][:Sn] = ScalarParameter(1.0, "", "Shape of density profile (1-x^2)^Sn")

# St
fuse_parameters[:PLASMA_PARAMETERS][:St] = ScalarParameter(1.0, "", "Shape of temperature profile (1-x^2)^St")

# Sj
fuse_parameters[:PLASMA_PARAMETERS][:Sj] = ScalarParameter(1.0, "", "Shape of current density profile (1-x^2)^Sj")

# define FuseParameters for rapid access of .value property
for (k, v) in fuse_parameters
    fuse_parameters[k] = FuseParameters(v)
end

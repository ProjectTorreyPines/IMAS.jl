# begin adding parameters to simple dictionary
fuse_parameters = Dict()

# bootstrapModel
options = Dict()
options[:pomphrey] = SwitchOption(CBSPomphrey(), "", "Pomphrey, N. Bootstrap dependence on plasma profile parameters. PPPL, 1992.")
options[:gi] = SwitchOption(CBSGi(), "", "Gi et al., Fus. Eng. Design 89 2709 (2014)")
options[:wilson] = SwitchOption(CBSWilson(), "", "Wilson et al., Nucl. Fusion 32 257 (1992)")
options[:user] = ScalarParameter(0.7, "", "User-defined constant")
fuse_parameters[:bootstrapModel] = SwitchParameter(options, :gi, "Bootstrap current efficiency model")


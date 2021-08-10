# begin adding parameters to simple dictionary
fuse_parameters = Dict()

# bootstrapModel
options = Dict()
options[:pomphrey] = SwitchOption(:pomphrey, "", "pomphrey")
options[:gi] = SwitchOption(:gi, "", "gi")
options[:wilson] = SwitchOption(:wilson, "", "wilson")
options[:user] = ScalarParameter(0.7, "", "User-defined constant")
fuse_parameters[:bootstrapModel] = SwitchParameter(options, :gi, "Bootstrap current efficiency model")

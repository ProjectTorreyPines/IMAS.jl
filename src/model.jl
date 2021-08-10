using FUSE

epsilon = 1.8
St=1.0
Sn=1.0
Sj=1.0
effectiveZ=2.0

bootstrapCoefficient = FUSE.collisionless_bootstrap(FUSE.fuse_parameters[:bootstrapModel], epsilon, St, Sn, Sj, effectiveZ)
println(bootstrapCoefficient)
using FUSE

epsilon = 1.8
St = FUSE.fuse_parameters[:PLASMA_PARAMETERS][:St]
Sn = FUSE.fuse_parameters[:PLASMA_PARAMETERS][:Sn]
Sj = FUSE.fuse_parameters[:PLASMA_PARAMETERS][:Sj]
effectiveZ = 2.0

bootstrapCoefficient = FUSE.collisionless_bootstrap(FUSE.fuse_parameters[:PHYSICS_MODELS][:bootstrapModel], epsilon, St, Sn, Sj, effectiveZ)
println(bootstrapCoefficient)
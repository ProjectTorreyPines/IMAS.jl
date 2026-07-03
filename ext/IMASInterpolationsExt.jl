module IMASInterpolationsExt

import IMAS
import Interpolations

# Interpolations.jl backend for the generic interpolant helpers in IMAS (src/physics/fields.jl).
# Active only when Interpolations is loaded. It uses positional coords and returns SVector/SMatrix,
# which index like the FastInterpolations NTuple/Matrix the IMAS math expects.
@inline IMAS._gradient(itp::Interpolations.AbstractInterpolation, r, z) =
    Interpolations.gradient(itp, r, z)

@inline IMAS._value_gradient(itp::Interpolations.AbstractInterpolation, r, z) =
    (itp(r, z), Interpolations.gradient(itp, r, z))

@inline IMAS._hessian!(H, itp::Interpolations.AbstractInterpolation, r, z) =
    (H .= Interpolations.hessian(itp, r, z); H)

end

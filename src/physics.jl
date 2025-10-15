import CoordinateConventions: cocos, COCOS
import Interpolations
import Contour
import StaticArrays
import PolygonOps
import Optim
import NLsolve
import Statistics
import DataInterpolations: DataInterpolations, ExtrapolationType
import IMASutils: IMASutils, trapz, cumtrapz, cumtrapz!, mirror_bound, argmin_abs
import PeriodicTable: elements
import MillerExtendedHarmonic: MXH, MXH!, flat_coeffs, reorder_flux_surface!
import Memoize
import Roots
using RecipesBase

include(joinpath("physics", "constants.jl"))
include(joinpath("physics", "outline.jl"))
include(joinpath("physics", "equilibrium.jl"))
include(joinpath("physics", "build.jl"))
include(joinpath("physics", "tf.jl"))
include(joinpath("physics", "currents.jl"))
include(joinpath("physics", "fields.jl"))
include(joinpath("physics", "fluxsurfaces.jl"))
include(joinpath("physics", "rf.jl"))
include(joinpath("physics", "neoclassical.jl"))
include(joinpath("physics", "profiles.jl"))
include(joinpath("physics", "sol.jl"))
include(joinpath("physics", "nuclear.jl"))
include(joinpath("physics", "sources.jl"))
include(joinpath("physics", "fast.jl"))
include(joinpath("physics", "transport.jl"))
include(joinpath("physics", "pedestal.jl"))
include(joinpath("physics", "radiation.jl"))
include(joinpath("physics", "collisions.jl"))
include(joinpath("physics", "boundary.jl"))
include(joinpath("physics", "particles.jl"))
include(joinpath("physics", "pf_active.jl"))
include(joinpath("physics", "technology.jl"))
include(joinpath("physics", "thermal_loads.jl"))
include(joinpath("physics", "diagnostics.jl"))
include(joinpath("physics", "interferometer.jl"))
include(joinpath("physics", "magnetics.jl"))

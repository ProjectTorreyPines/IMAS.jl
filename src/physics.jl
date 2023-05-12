import CoordinateConventions: cocos, COCOS
import Interpolations
import Contour
import StaticArrays
import PolygonOps
import Optim
import NumericalIntegration: integrate, cumul_integrate
import PeriodicTable: elements
import MillerExtendedHarmonic: MXH, MXH!, flat_coeffs, reorder_flux_surface!
import Memoize
import Roots
using RecipesBase

@enum BuildLayerType::Int _plasma_ = -1 _gap_ _oh_ _tf_ _shield_ _blanket_ _wall_ _vessel_ _cryostat_ _divertor_
@enum BuildLayerSide::Int _lfs_ = -1 _lhfs_ _hfs_ _in_ _out_
@enum BuildLayerShape::Int _offset_ _negative_offset_ _convex_hull_ _princeton_D_exact_ _princeton_D_ _princeton_D_scaled_ _rectangle_ _double_ellipse_ _triple_arc_ _miller_ _square_miller_ _spline_ _silo_

include("physics/equilibrium.jl")
include("physics/build.jl")
include("physics/currents.jl")
include("physics/fluxsurfaces.jl")
include("physics/misc.jl")
include("physics/neoclassical.jl")
include("physics/profiles.jl")
include("physics/sol.jl")
include("physics/sources.jl")
include("physics/gacode_constants.jl")
include("physics/transport.jl")
include("physics/pedestal.jl")
include("physics/radiation.jl")
include("physics/collisions.jl")
include("physics/boundary.jl")

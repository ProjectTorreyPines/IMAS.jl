# PC vs adaptive-RK4 vs Contour.jl: drift, point count, timing on DIII-D closed surfaces.
# Scratch diagnostic (not a pass/fail test).
import IMAS
using Printf

filename = joinpath(pkgdir(IMAS.IMASdd), "sample", "D3D_eq_ods.json")
dd = IMAS.json2imas(filename; show_warnings=false)
eqt = dd.equilibrium.time_slice[1]
eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
r, z, itp = IMAS.ψ_interpolant(eqt2d)
RA, ZA = eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z
psi_axis = itp(RA, ZA)
eqt1d = eqt.profiles_1d
# a mid-radius closed level
c = psi_axis + 0.5 * (eqt1d.psi[end] - psi_axis)
seed = (maximum(r) - 1e-3, ZA)   # outboard; project finds the surface
seed = IMAS._project_to_level(itp, c, seed)[1]

for method in (:pc, :rk4)
    t = @elapsed (Rs, Zs, closed) = IMAS._trace_surface_cubic(itp, c, seed; method)
    drift = maximum(abs(IMAS.FI.value_gradient(itp,(Rs[k],Zs[k]))[1]-c) for k in eachindex(Rs))
    @printf("%-5s closed=%s  N=%4d  drift=%.2e  t=%.2f ms\n", method, closed, length(Rs), drift, t*1e3)
end

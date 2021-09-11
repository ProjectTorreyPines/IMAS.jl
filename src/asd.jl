using Pkg
Pkg.activate("/Users/meneghini/.julia/dev/IMAS")
# using Revise
using IMAS
# println(VERSION)

time_slice = IMAS.equilibrium__time_slice()
profiles_1d = time_slice.profiles_1d

using ForwardDiff
#using Trapz

n=11
profiles_1d.psi = collect(range(0.0, 1.0, length=n))
profiles_1d.pressure = p0 = collect(1.0.-range(0.0, 1.0, length=n).^2)
profiles_1d.volume = collect(range(0.0, 1.0, length=n))

function test_gradient(mul)
    profiles_1d.pressure = profiles_1d.pressure.*mul.^2
    sum(time_slice.global_quantities.energy_mhd)
end

#test_gradient(3)
println(ForwardDiff.derivative(test_gradient, 2.0))

# using Pkg
# Pkg.activate("/Users/meneghini/.julia/dev/IMAS")
# using Revise
# using IMAS
# println(VERSION)


# time_slice = IMAS.equilibrium__time_slice()
# profiles_1d = time_slice.profiles_1d

# n=11
# profiles_1d.psi = collect(range(0.0, 1.0, length=n))
# profiles_1d.pressure = collect(1.0.-range(0.0, 1.0, length=n).^2)
# profiles_1d.volume = collect(range(0.0, 1.0, length=n))
# time_slice.global_quantities.energy_mhd

# time_slice = IMAS.equilibrium__time_slice()
# profiles_1d = time_slice.profiles_1d

# time_slice = IMAS.equilibrium__time_slice()
# profiles_1d = time_slice.profiles_1d

# #using Zygote
# using ForwardDiff
# #using Trapz

# n=11
# profiles_1d.psi = collect(range(0.0, 1.0, length=n))
# profiles_1d.pressure = p0 = collect(1.0.-range(0.0, 1.0, length=n).^2)
# profiles_1d.volume = collect(range(0.0, 1.0, length=n))

# function test_gradient(mul)
#     profiles_1d.pressure = profiles_1d.pressure.*mul.^2
#     sum(time_slice.global_quantities.energy_mhd)
# end

# #test_gradient(3)
# ForwardDiff.derivative(test_gradient, 2)
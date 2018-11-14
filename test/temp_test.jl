
include("../src/NBIPsMulti.jl")

using JuLIP, NBodyIPs, NBodyIPFitting, NBIPsMulti, FileIO

using NBodyIPs: eval_site_nbody!

include(homedir() * "/Gits/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz() # ; include = ["hess_bcc", "hess_hcp"])
@show length(data)
data2 = data
at2 = data[1].at

r0 = 2.5
BL2 = BondLengthDesc("exp( - 2 * (r/$r0-1))",
                    "(:cos, $(r0-1.5), $(r0))")

basis = [
     nbpolys(2, BL2, 14, [29,29]);
  ]

at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)

# @show E = energy(basis[1], at)
# @show E2 = energy(basis[1],at2)

@show Ei = site_energies(basis[1],at,[29,29])
@show Ei = site_energies(basis[1],at2,[29,29])


V = basis[1]
Z = atomic_numbers(at)
Species = [29,29]
Es = zeros(Float64, length(at))
for (i, j, r, R) in sites(at, cutoff(V))
   Spi = Z[i]
   Spj = Z[j]
   @which eval_site_nbody!(Val(2), R, cutoff(V),
                            ((out, R, J, temp,Spi,Spj,Species) -> out + evaluate(V, descriptor(V), R, J,Spi,Spj,Species)), zero(Float64), nothing, Spi,Spj,Species)
   Es[i] = eval_site_nbody!(Val(2), R, cutoff(V),
                            ((out, R, J, temp,Spi,Spj,Species) -> out + evaluate(V, descriptor(V), R, J,Spi,Spj,Species)), zero(Float64), nothing, Spi,Spj,Species)
end

@which site_energies(basis[1],at,[29,29])
@which site_energies(basis[1],at2,[29,29])

@which energy(basis[1], at)
@which energy(basis[1], at2)

typeof(at)
typeof(at2)

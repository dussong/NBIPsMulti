
include("../src/NBIPsMulti.jl")

using JuLIP, NBodyIPs, FileIO, NBIPsMulti, NBodyIPFitting

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

@which site_energies(basis[1],at,[29,29])
@which site_energies(basis[1],at2,[29,29])

@which energy(basis[1], at)
@which energy(basis[1], at2)

typeof(at)
typeof(at2)

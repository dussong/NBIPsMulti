
include("../src/NBIPsMulti.jl")

using JuLIP, NBodyIPs, FileIO, NBIPsMulti, NBodyIPFitting

include(homedir() * "/Gits/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz() # ; include = ["hess_bcc", "hess_hcp"])
@show length(data)
data2 = data

r0 = 2.5
BL2 = BondLengthDesc("exp( - 2 * (r/$r0-1))",
                    "(:cos, $(r0-1.5), $(r0))")

basis = [
     nbpolys(2, BL2, 14, [29,29]);
  ]

at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)

E = energy(basis[1], at)
E2 = energy(basis[1],data[1].at)


@which energy(basis[1], at)

typeof(at)
typeof(data[1].at)

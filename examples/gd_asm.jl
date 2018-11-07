##
info("Load libraries ...")

using JuLIP, NBodyIPFitting, NBodyIPs, NBIPsMulti

info("Load Butane database ...")

include(homedir() * "/Gits/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz() # ; include = ["hess_bcc", "hess_hcp"])
@show length(data)


r0 = round(rnn(:Ti),2)
E0 = Butane.get_E0()

rcut2 = 2.8 * r0
rcut3 = 2.3 * r0
rcut4 = 1.9 * r0

info("Generate descriptors...")

BL2 = BondLengthDesc("exp( - 2 * (r/$r0-1))",
                    "(:cos, $(rcut2-1.5), $(rcut2))")
BL3 = BondLengthDesc("exp( - 2.5 * (r/$r0-1))",
                   "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))")


##
info("Generate a BL-3B basis ...")
basis = [
      nbpolys(2, BL2a, 14);
      nbpolys(2, BL2b, 8);
      nbpolys(3, BL3, 12);
      nbpolys(4, BL4, 10);
   ]
info("Assemble the LsqDB ...")
@show length(basis)
dbpath = homedir() * "/scratch/nbodyips/Ti_4B_BL_med"
LsqDB(dbpath, basis, data)


##
info("Generate a BA-4B basis ...")
basis = [
      nbpolys(2, BA2a, 14);
      nbpolys(2, BA2b, 8);
      nbpolys(3, BA3, 12);
      nbpolys(4, BA4, 10);
   ]
info("Assemble the LsqDB ...")
@show length(basis)
dbpath = homedir() * "/scratch/nbodyips/Ti_4B_BA_med"
LsqDB(dbpath, basis, data)

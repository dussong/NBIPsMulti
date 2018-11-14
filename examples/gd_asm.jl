##
info("Load libraries ...")

include("../src/NBIPsMulti.jl")

using JuLIP, NBodyIPs, FileIO, NBIPsMulti, NBodyIPFitting

info("Load Butane database ...")

include(homedir() * "/Gits/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz() # ; include = ["hess_bcc", "hess_hcp"])
@show length(data)
data2 = data


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
info("Generate a BL-2B basis ...")
basis = [
      # nbpolys(2, BL2, 10, [1,1]);
      # nbpolys(2, BL2, 10, [1,6]);
      nbpolys(2, BL2, 10, [6,6]);
      # nbpolys(3, BL3, 5, [1,6,6]);
   ]
# basis = [
#       nbpolys(2, BL2, 14);
#    ]

# @show energy(basis,data[1].at)

energy(basis[1],data[1].at)


info("Assemble the LsqDB ...")
@show length(basis)
dbpath = homedir() * "/Dropbox/PIBmat/multispecies/Butane_3B"
# dbpath = homedir() * "/Documents/Butane_3B"


db =  LsqDB(dbpath, basis, data2);



#
# ##
# info("Generate a BA-4B basis ...")
# basis = [
#       nbpolys(2, BA2a, 14);
#       nbpolys(2, BA2b, 8);
#       nbpolys(3, BA3, 12);
#       nbpolys(4, BA4, 10);
#    ]
# info("Assemble the LsqDB ...")
# @show length(basis)
# dbpath = homedir() * "/scratch/nbodyips/Ti_4B_BA_med"
# LsqDB(dbpath, basis, data)

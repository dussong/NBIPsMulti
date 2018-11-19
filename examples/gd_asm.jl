##
info("Load libraries ...")

include("../src/NBIPsMulti.jl")

using JuLIP, NBodyIPs, NBIPsMulti, NBodyIPFitting

info("Load Butane database ...")

include(homedir() * "/Gits/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz() # ; include = ["hess_bcc", "hess_hcp"])
@show length(data)


r0 = 3.*round(rnn(:C),2)
E0 = Butane.get_E0()

rcut2 = 2.8 * r0
rcut3 = 2.3 * r0
rcut4 = 1.9 * r0

info("Generate descriptors...")

BL2 = MultiDesc("exp( - 2 * (r/$r0-1))", "(:cos, $(rcut2-1.5), $(rcut2))",Val(:AA))

# BL3_AAA = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAA))

BL3_AAA = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos, $(rcut2-1.5), $(rcut2))",Val(:AAA))

BL3_AAB = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos, $(rcut2-1.5), $(rcut2))",Val(:AAB))


##
info("Generate a BL-2B basis ...")
basis = [
      nbpolys(BL2, 14, [6,6]);
      nbpolys(BL2, 14, [1,1]);
      nbpolys(BL2, 14, [1,6]);
      nbpolys(BL3_AAA, 3, [1,1,1]);
      nbpolys(BL3_AAA, 3, [6,6,6]);
      nbpolys(BL3_AAB, 3, [1,1,6]);
      nbpolys(BL3_AAB, 3, [1,6,6]);
   ]

info("Assemble the LsqDB ...")
@show length(basis)
dbpath = homedir() * "/Gits/NBIPsMulti/data/Butane_3B"

db =  LsqDB(dbpath, basis, data);

# Check that energy and forces are non zero.
# db1.kron_groups[":14EF"]["F"]
# db1.kron_groups[":14EF"]["E"]

# -------------------------
# To do the fit right away
# -------------------------
info("Fit Butane Database basis...")

dataweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)
configweights = Dict(""  => 1.0)

IP, info = lsqfit( db; E0 = E0,
                       dataweights=dataweights,
                       configweights=configweights,
                       # Ibasis = Ibasis
                       )
info

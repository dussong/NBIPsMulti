##
@info("Load libraries ...")

using JuLIP, NBodyIPs, NBIPsMulti, NBodyIPFitting

@info("Load Butane database ...")

include(homedir() * "/Gits/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz() # ; include = ["hess_bcc", "hess_hcp"])
data
@show length(data)
data[1]

r0 = 3*round(rnn(:C),digits=2)
E0 = Butane.get_E0()

rcut2 = 2.8 * r0
rcut3 = 2.3 * r0
rcut4 = 1.9 * r0

@info("Generate descriptors...")

BL2 = MultiDesc("exp( - 2 * (r/$r0-1))", "(:cos, $(rcut2-1.5), $(rcut2))",Val(:AA))

BL3_AAA = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAA))

BL3_AAB = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAB))

BL4_AAAA = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAAA))

BL4_AAAB = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAAB))

BL4_AABB = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AABB))

##
@info("Generate a 4B basis ...")
basis = [
      nbpolys(BL2, 14, [6,6]);
      nbpolys(BL2, 14, [1,1]);
      nbpolys(BL2, 14, [1,6]);
      nbpolys(BL3_AAA, 5, [1,1,1]);
      nbpolys(BL3_AAA, 5, [6,6,6]);
      nbpolys(BL3_AAB, 5, [1,1,6]);
      nbpolys(BL3_AAB, 5, [1,6,6]);
      nbpolys(BL4_AAAA, 3, [1,1,1,1]);
      nbpolys(BL4_AAAA, 3, [6,6,6,6]);
      nbpolys(BL4_AAAB, 3, [1,1,1,6]);
      nbpolys(BL4_AAAB, 3, [1,6,6,6]);
      nbpolys(BL4_AABB, 3, [1,1,6,6]);
   ]

@info("Assemble the LsqDB ...")
@show length(basis)
dbpath = homedir() * "/Gits/NBIPsMulti/data/Butane_4B"

db =  LsqDB(dbpath, basis, data);
db

# Check that energy and forces are non zero.
# db1.kron_groups[":14EF"]["F"]
# db1.kron_groups[":14EF"]["E"]

# -------------------------
# To do the fit right away
# -------------------------
@info("Fit Butane Database basis...")

p = 1.

obsweights = Dict("E" => 1.0, "F" => 1.0)
configweights = Dict(""  => 1.0)

Y, W = NBodyIPFitting.Lsq.collect_observations(db, configweights, dataweights, OneBody(E0))


IP, lsqinfo = lsqfit( db; E0 = E0,
                       obsweights=obsweights,
                       configweights=configweights,
                       # Ibasis = collect(1:45)
                       Vref = OneBody(E0)
                       )

errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)

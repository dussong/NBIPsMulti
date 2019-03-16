##
@info("Load libraries ...")

using JuLIP, NBodyIPs, NBIPsMulti, NBodyIPFitting

@info("Load Butane database ...")

include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz()
data
@show length(data)
data[1]

r0 = 3 .*round(rnn(:C),digits=2)
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

Vn = ("1.0/(1.0 + exp(1.0 / ($(1.5*r0) - r + 1e-2)))", 1.5*r0)


##
@info("Generate a 4B env basis ...")

basis = [
      envpolysM(BL2, 5, Vn, 2, [6,6]);
      envpolysM(BL2, 5, Vn, 2, [1,1]);
      envpolysM(BL2, 5, Vn, 2, [1,6]);
      envpolysM(BL3_AAA, 3, Vn, 1, [1,1,1]);
      envpolysM(BL3_AAA, 3,Vn, 1, [6,6,6]);
      envpolysM(BL3_AAB, 3,Vn, 1, [1,1,6]);
      envpolysM(BL3_AAB, 3,Vn, 1, [1,6,6]);
      envpolysM(BL4_AAAA, 2,Vn, 1, [1,1,1,1]);
      envpolysM(BL4_AAAA, 2,Vn, 1, [6,6,6,6]);
      envpolysM(BL4_AAAB, 2,Vn, 1, [1,1,1,6]);
      envpolysM(BL4_AAAB, 2,Vn, 1, [1,6,6,6]);
      envpolysM(BL4_AABB, 2,Vn, 1, [1,1,6,6]);
   ]

@info("Assemble the LsqDB ...")
@show length(basis)
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_env"

db =  LsqDB(dbpath, basis, data);
db

# -------------------------
# To do the fit right away
# -------------------------
@info("Fit Butane Database basis...")

obsweights = Dict("E" => 1.0, "F" => 1.0)
configweights = Dict(""  => 1.0)



IP, lsqinfo = lsqfit( db; E0 = E0,
                       obsweights=obsweights,
                       configweights=configweights,
                       # Ibasis = collect(1:45)
                       )

errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)

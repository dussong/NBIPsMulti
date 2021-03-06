##
@info("Load libraries ...")

using JuLIP, NBodyIPs, NBIPsMulti, IPFitting

@info("Load Butane database ...")

# Loading the data
include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz()
@show length(data)


r0 = 3*round(rnn(:C),digits=2) #cutoff parameter
rcut2 = 2.8 * r0
rcut3 = 2.3 * r0
rcut4 = 1.9 * r0


@info("Generate descriptors...")

#Two-body descriptor
BL2 = MultiDesc("exp( - 2 * (r/$r0-1))", CosCut(rcut2-1.5,rcut2),Val(:AA))

#Three-body descriptors
BL3_AAA = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAA))

BL3_AAB = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAB))

#Four-body descriptors
BL4_AAAA = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAAA))

BL4_AAAB = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAAB))

BL4_AABB = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AABB))


##
@info("Generate a 4B basis ...")
#1 is for Hydrogen, 6 is for Carbon
basis = [
      nbpolys(BL2, 14, [6,6]);
      nbpolys(BL2, 14, [1,1]);
      nbpolys(BL2, 14, [1,6]);
      nbpolys(BL3_AAA, 3, [1,1,1]);
      nbpolys(BL3_AAA, 3, [6,6,6]);
      nbpolys(BL3_AAB, 3, [1,1,6]);
      nbpolys(BL3_AAB, 3, [1,6,6]);
      nbpolys(BL4_AAAA, 2, [1,1,1,1]);
      nbpolys(BL4_AAAA, 2, [6,6,6,6]);
      nbpolys(BL4_AAAB, 2, [1,1,1,6]);
      nbpolys(BL4_AAAB, 2, [1,6,6,6]);
      nbpolys(BL4_AABB, 2, [1,1,6,6]);
   ]

@info("Assemble the LsqDB ...")
@show length(basis)
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"

db =  LsqDB(dbpath, basis, data);

# -------------------------
# To do the fit right away
# -------------------------
@info("Fit Butane Database basis...")

for c in db.configs
    c.configtype = "test"
end

obsweights = Dict("E" => 1.0, "F" => 1.0)
configweights = Dict("test"  => 1.0)

# Reference energy
E0 = Butane.get_E0() #reference energy
IP, lsqinfo = lsqfit( db;
                       obsweights=obsweights,
                       configweights=configweights,
                       Ibasis = collect(1:45), #selecting particular basis functions
                       combineIP = NBodyIP,
                       # solver = (:rrqr, 1e-12),
                       Vref = OneBody(E0) #reference potential
                       )


errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)


# Save potential
save_ip(homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_IP.json", IP, lsqinfo); #save potential

# Load the potential directly
load_ip(homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_IP.json")

# For faster evaluation
fast(IP);

@info("Load libraries...")

using NBodyIPs
using IPFitting
using JuLIP
using NBIPsMulti
using StaticArrays

include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"

data = Butane.load_xyz()

r0 = 3*round(rnn(:C),digits=2)
E0 = Butane.get_E0()

rcut2 = 2.8 * r0
# rcut3 = 2.3 * r0
# rcut4 = 1.9 * r0

@info("Generate descriptors...")

BL2 = MultiDesc(ExpTransform(2, r0), CosCut(rcut2-1.5,rcut2),Val(:AA))

# BL3_AAA = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut3-1.5,rcut3),Val(:AAA))
#
# BL3_AAB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut3-1.5,rcut3),Val(:AAB))
#
# BL4_AAAA = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AAAA))
#
# BL4_AAAB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AAAB))
#
# BL4_AABB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AABB))


@info("Generate a 4B basis ...")
basis = [
      nbpolys(BL2, 14, [6,6]);
      nbpolys(BL2, 14, [1,1]);
      nbpolys(BL2, 14, [1,6]);
      # nbpolys(BL3_AAA, 3, [1,1,1]);
      # nbpolys(BL3_AAA, 3, [6,6,6]);
      # nbpolys(BL3_AAB, 3, [1,1,6]);
      # nbpolys(BL3_AAB, 3, [1,6,6]);
      # nbpolys(BL4_AAAA, 2, [1,1,1,1]);
      # nbpolys(BL4_AAAA, 2, [6,6,6,6]);
      # nbpolys(BL4_AAAB, 2, [1,1,1,6]);
      # nbpolys(BL4_AAAB, 2, [1,6,6,6]);
      # nbpolys(BL4_AABB, 2, [1,1,6,6]);
   ]

@info("Assemble the LsqDB ...")
@show length(basis)
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"

db =  LsqDB(dbpath, basis, data);
db

# db = LsqDB(dbpath)
# db
# summary(db)
##
@info("Fit Butane Database basis...")

obsweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)
configweights = Dict("nothing"  => 1.0)

oneB = MOneBody(Dict("E0" => Dict(:H => -13.5203101677, :C => -1027.28368153)))

IP, lsqinfo = lsqfit( db;
                       # E0 = 1.,
                       obsweights=obsweights,
                       configweights=configweights,
                       solver = (:rrqr, 1e-16),
                       # solver = (:qr,),
                       combineIP = NBodyIP,
                       # Ibasis = Ibasis
                       Vref = oneB
                       )

errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)

## 2B repulsion
using NBIPsMulti.RepulsionM: RepulsiveCoreM

V2_1 = IP.components[2]
V2_2 = IP.components[3]
V2_3 = IP.components[4]

V2rep_1 = RepulsiveCoreM(V2_1, 0.0, -1.)
V2rep_2 = RepulsiveCoreM(V2_2, 0.0, -1.)
V2rep_3 = RepulsiveCoreM(V2_3, 0.0, -1.)

IP_rep = deepcopy(IP)

IP_rep.components[2] = V2rep_1
IP_rep.components[3] = V2rep_2
IP_rep.components[4] = V2rep_3


atX = SArray{Tuple{3},Float64,1,3}[[0.0, 0.0, 0.0], [0.6, 0.0, 0.0]]
atZ = [6, 1]
at = Atoms(atZ, atX; cell = [10 0 0;0 10 0;0 0 10])
energy(IP_rep, at)
energy(IP, at)

energy(IP_rep, db.configs[1].at)
energy(IP, db.configs[1].at)
db.configs[1].D["E"][1]

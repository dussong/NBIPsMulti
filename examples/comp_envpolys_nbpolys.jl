@info("Load libraries...")

using NBodyIPs
using IPFitting
using JuLIP
using NBIPsMulti
using StaticArrays


# -------------------------
# For the 2 cases
# -------------------------
include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"

data = Butane.load_xyz()

r0 = 3*round(rnn(:C),digits=2)
E0 = Butane.get_E0()

rcut2 = 2.8 * r0
rcut3 = 2.3 * r0
# rcut4 = 1.9 * r0

@info("Generate descriptors...")

BL2 = MultiDesc(ExpTransform(2, r0), CosCut(rcut2-1.5,rcut2),Val(:AA))
BL3_AAA = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut3-1.5,rcut3),Val(:AAA))

BL3_AAB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut3-1.5,rcut3),Val(:AAB))
#
# BL4_AAAA = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AAAA))
#
# BL4_AAAB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AAAB))
#
# BL4_AABB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AABB))

obsweights = Dict("E" => 1.0, "F" => 1.0)
configweights = Dict("nothing"  => 1.0)

oneB = MOneBody(Dict("E0" => Dict(:H => -13.5203101677, :C => -1027.28368153)))

# -------------------------

# -------------------------
# For the 2 cases
# -------------------------

@info("Generate a 4B basis with nbpolys ...")
basisnb = [
      nbpolys(BL2, 14, [6,6]);
      nbpolys(BL2, 14, [1,1]);
      nbpolys(BL2, 14, [1,6]);
      nbpolys(BL3_AAA, 3, [1,1,1]);
      nbpolys(BL3_AAA, 3, [6,6,6]);
      nbpolys(BL3_AAB, 3, [1,1,6]);
      nbpolys(BL3_AAB, 3, [1,6,6]);
      # nbpolys(BL4_AAAA, 2, [1,1,1,1]);
      # nbpolys(BL4_AAAA, 2, [6,6,6,6]);
      # nbpolys(BL4_AAAB, 2, [1,1,1,6]);
      # nbpolys(BL4_AAAB, 2, [1,6,6,6]);
      # nbpolys(BL4_AABB, 2, [1,1,6,6]);
   ]



@info("Assemble the LsqDB ...")
@show length(basisnb)
dbpathnb = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"

dbnbpol =  LsqDB(dbpathnb, basisnb, data);
dbnbpol

# db = LsqDB(dbpath)
# db
# summary(db)
##
@info("Fit Butane Database basis...")



IPnbpol, lsqinfonbpol = lsqfit( dbnbpol;
                       # E0 = 1.,
                       obsweights=obsweights,
                       configweights=configweights,
                       solver = (:rrqr, 1e-16),
                       # solver = (:qr,),
                       combineIP = NBodyIP,
                       # Ibasis = Ibasis
                       Vref = oneB
                       )

errsnbpol = lsqinfonbpol["errors"]
rmse_table(rmse(errsnbpol)...)

# -------------------------
# Environment-dependence
# -------------------------


Vn = ("1.0/(1.0 + exp(1.0 / ($(1.5*r0) - r + 1e-2)))", 1.5*r0)

weights = Dict( (1,1)=> 1.0,
                (1,6)=> 2.0,
                (6,1)=> 2.0,
                (6,6)=> 5.0
                )

basisenv = [
    envpolysM(BL2, 14, Vn, 0, weights, [6,6]);
    envpolysM(BL2, 14, Vn, 0, weights, [1,1]);
    envpolysM(BL2, 14, Vn, 0, weights, [1,6]);
    envpolysM(BL3_AAA, 3,Vn, 0, weights, [1,1,1]);
    envpolysM(BL3_AAA, 3,Vn, 0, weights, [6,6,6]);
    envpolysM(BL3_AAB, 3,Vn, 0, weights, [1,1,6]);
    envpolysM(BL3_AAB, 3,Vn, 0, weights, [1,6,6]);
    # envpolysM(BL4_AAAA, 4,Vn, 0, weights, [1,1,1,1]);
    # envpolysM(BL4_AAAA, 4,Vn, 0, weights, [6,6,6,6]);
    # envpolysM(BL4_AAAB, 4,Vn, 0, weights, [1,1,1,6]);
    # envpolysM(BL4_AAAB, 4,Vn, 0, weights, [1,6,6,6]);
    # envpolysM(BL4_AABB, 4,Vn, 0, weights, [1,1,6,6]);
 ]

@info("Assemble the environment LsqDB ...")
@show length(basisenv)
dbpathenv = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_env3"

dbenv =  LsqDB(dbpathenv, basisenv, data);
dbenv


@info("Fit Butane Database basis with env dependence...")

IPenv, lsqinfoenv = lsqfit( dbenv;
                    # E0 = 1.,
                    obsweights=obsweights,
                    configweights=configweights,
                    solver = (:rrqr, 1e-16),
                    # solver = (:qr,),
                    combineIP = NBodyIP,
                    # Ibasis = Ibasis
                    Vref = oneB
                    )

errsenv = lsqinfoenv["errors"]
rmse_table(rmse(errsenv)...)


atX = SArray{Tuple{3},Float64,1,3}[[0.0, 0.0, 0.0], [0.6, 0.0, 0.0]]
atZ = [6, 1]
at = Atoms(atZ, atX; cell = [10 0 0;0 10 0;0 0 10])
energy(IPnbpol, at)
energy(IPenv, at)

energy(IPnbpol, data[1].at)
energy(IPenv, data[1].at)
data[1].D["E"][1]

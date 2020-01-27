@info("Load libraries...")

using NBodyIPs
using IPFitting
using JuLIP
using NBIPsMulti
using NBodyIPs.Regularisers: BLRegulariser, EnvBLRegulariser


include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
data = Butane.load_xyz()

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"

@info("Generate descriptors...")
r0 = 3*round(rnn(:C),digits=2) #cutoff parameter
rcut0 = 1.5 * r0
rcut2 = 2.8 * r0
rcut3 = 2.3 * r0
rcut4 = 1.9 * r0
BL2 = MultiDesc(ExpTransform(2, r0), CosCut(rcut2-1.5,rcut2),Val(:AA))
BL3_AAA = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut3-1.5,rcut3),Val(:AAA))

BL3_AAB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut3-1.5,rcut3),Val(:AAB))
#
# BL4_AAAA = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AAAA))
#
# BL4_AAAB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AAAB))
#
# BL4_AABB = MultiDesc(ExpTransform(2.5,r0), CosCut(rcut4-1.5,rcut4),Val(:AABB))


# Environment function
Vn = ("1.0/(1.0 + exp(1.0 / ($(1.5*r0) - r + 1e-2)))", 1.5*r0)
# Environment weights
weights = Dict( (1,1)=> 1.0,
                (1,6)=> 2.0,
                (6,1)=> 2.0,
                (6,6)=> 5.0
                )

# Basis
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

# Regularisers with environment
 regenv = [ EnvBLRegulariser(2, 0, rcut0, rcut3, npoints = 30, creg = 1.0,
                 transform = PolyTransform(1, r0), species = [1,6]),
         EnvBLRegulariser(2, 0, rcut0, rcut3, npoints = 30, creg = 1.0,
                         transform = PolyTransform(1, r0), species = [1,1]),
         EnvBLRegulariser(2, 0, rcut0, rcut3, npoints = 30, creg = 1.0,
                         transform = PolyTransform(1, r0), species = [6,6])                         ]


@info("Assemble the environment LsqDB ...")
@show length(basisenv)
dbpathenv = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_env"

dbenv =  LsqDB(dbpathenv, basisenv, data);
dbenv


@info("Fit Butane Database basis with env dependence...")

#Weights on the observations
obsweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)

# changing configuration names
for c in dbenv.configs
    c.configtype = "test"
end
#configuration weights
configweights = Dict("test"  => 1.0)

E0 = Butane.get_E0() #reference energy

IPenv, lsqinfoenv = lsqfit( dbenv;
                    # E0 = 1.,
                    obsweights=obsweights,
                    configweights=configweights,
                    solver = (:rrqr, 1e-16),
                    # solver = (:qr,),
                    combineIP = NBodyIP,
                    regularisers = regenv,
                    # Ibasis = Ibasis
                    Vref = OneBody(E0)
                    )

errsenv = lsqinfoenv["errors"]
rmse_table(rmse(errsenv)...)

@info("Load libraries...")

using NBodyIPs, IPFitting, JuLIP, NBIPsMulti
using NBodyIPs.Regularisers
include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
E0 = Butane.get_E0()

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_env3"


db = LsqDB(dbpath)
db
# summary(db)
##
@info("Fit Butane Database basis...")

obsweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)
configweights = Dict(""  => 1.0)

rcut0=0.
rcut3=6.
r0 = 2.
reg = [ EnvBAReg(4, 1, rcut0, rcut3, npoints = 6000, creg = 1000., transform = PolyTransform(2, r0), species = [1,6,6,6]) ]

IP, lsqinfo = lsqfit( db; )

IP, lsqinfo = lsqfit( db; E0 = E0,
                       obsweights=obsweights,
                       configweights=configweights,
                       # solver = (:rrqr, 1e-16),
                       solver = (:qr,),
                       combineIP = NBodyIP,
                       # regularisers = reg,
                       # Ibasis = Ibasis
                       )


errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)

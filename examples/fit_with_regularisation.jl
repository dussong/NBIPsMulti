@info("Load libraries...")

using NBodyIPs, IPFitting, JuLIP, NBIPsMulti
using NBodyIPs.Regularisers: BLRegulariser, EnvBLRegulariser
using NBodyIPs: bodyorder
using LinearAlgebra

# Load data
include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
E0 = Butane.get_E0()

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"

# Load database
db = LsqDB(dbpath)
##
@info("Fit Butane Database basis...")

# Change configuration name
for c in db.configs
    c.configtype = "test"
end

obsweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0) #weights on observations
configweights = Dict("test"  => 1.0) #weights on configuration

#Definition of the regularisers
rcut0 = 1.4
rcut3 = 10.
r0 = 2.

reg = [ BLRegulariser(2, rcut0, rcut3, npoints = 30, creg = 1.0,
                transform = PolyTransform(1, r0), species = [1,6]),
        BLRegulariser(2, rcut0, rcut3, npoints = 30, creg = 1.0,
                        transform = PolyTransform(1, r0), species = [1,1]),
        BLRegulariser(2, rcut0, rcut3, npoints = 30, creg = 1.0,
                        transform = PolyTransform(1, r0), species = [6,6]),
        BLRegulariser(3, rcut0, rcut3, npoints = 100, creg = 1.0,
                        transform = PolyTransform(1, r0), species = [1,1,1])                 ]


IP, lsqinfo = lsqfit( db; E0 = E0,
                       obsweights=obsweights,
                       configweights=configweights,
                       solver = (:rrqr, 1e-12),
                       # solver = (:qr,),
                       combineIP = NBodyIP,
                       regularisers = reg,
                       # Ibasis = Ibasis
                       )


errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)

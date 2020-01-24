@info("Load libraries...")

using NBodyIPs, IPFitting, JuLIP, NBIPsMulti
using NBodyIPs.Regularisers
using NBodyIPs: bodyorder
using LinearAlgebra
include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
E0 = Butane.get_E0()

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_env3"

db = LsqDB(dbpath)
db
# summary(db)
##
@info("Fit Butane Database basis...")

for c in db.configs
    c.configtype = "test"
end


obsweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)
configweights = Dict("test"  => 1.0)

rcut0=1.4
rcut3=10.
r0 = 2.
reg = [ EnvBARegulariser(4, 0, rcut0, rcut3, npoints = 30, creg = 1.0,
                 transform = PolyTransform(1, r0), species = [1,6,6,6]) ]
# envdeg(b) = b.t
# I2 = findall( b -> bodyorder(b) == 2, db.basis )
# I3 = findall( b -> bodyorder(b) == 3, db.basis )
# I4 = findall( b -> bodyorder(b) == 4, db.basis )
# Ibasis = 1:370

IP, lsqinfo = lsqfit( db; E0 = E0,
                       obsweights=obsweights,
                       configweights=configweights,
                       solver = (:rrqr, 1e-12),
                       # solver = (:qr,),
                       combineIP = NBodyIP,
                       regularisers = reg,
                       # Ibasis =
                       )


errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)



A = db.Ψ
qrA = qr(A)
qrA.R
qrA.Q

norm(qrA.R[:, 371:end])
444-371
norm(A[:, 371:end])

@info("Load libraries...")

using NBodyIPs, NBodyIPFitting, JuLIP, NBIPsMulti
include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
E0 = Butane.get_E0()

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_env"

db = LsqDB(dbpath)
db
# summary(db)
##
@info("Fit Butane Database basis...")

obsweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)
configweights = Dict(""  => 1.0)

IP, lsqinfo = lsqfit( db; E0 = E0,
                       obsweights=obsweights,
                       configweights=configweights,
                       # Ibasis = Ibasis
                       )


errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)

# NBodyIPs.save_ip("Butane_4B_env.json", IP, info)

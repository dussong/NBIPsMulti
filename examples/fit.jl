@info("Load libraries...")

using NBodyIPs, IPFitting, JuLIP, NBIPsMulti
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

# IP = NBodyIP([IP.components[2]])

errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)

save_ip(homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_env_IP.json", IP, lsqinfo)

load_ip(homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_env_IP.json")

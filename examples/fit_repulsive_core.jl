@info("Load libraries...")

using NBodyIPs, IPFitting, JuLIP, NBIPsMulti

include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
E0 = Butane.get_E0()

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"


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

save_ip(homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_IP.json", IP, lsqinfo)


## 2B repulsion
using NBodyIPs.Repulsion: RepulsiveCore

V2_1 = IP.components[2]
V2_2 = IP.components[3]
V2_3 = IP.components[4]

r0 = 3.
V2rep_1 = RepulsiveCore(V2_1, 0.89*r0, -1.)
V2rep_2 = RepulsiveCore(V2_2, 0.89*r0, -1.)
V2rep_3 = RepulsiveCore(V2_3, 0.89*r0, -1.)


IP.components[2] = V2rep_1
IP.components[3] = V2rep_2
IP.components[4] = V2rep_3
save_ip("pot_with_rep_core.json", IP, lsqinfo)

db.configs[1].at
energy(IP, db.configs[1].at)

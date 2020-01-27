@info("Load libraries...")

using NBodyIPs, IPFitting, JuLIP, NBIPsMulti

# Load data
include(homedir() * "/.julia/dev/NBIPsMulti/src/Butane.jl")
E0 = Butane.get_E0()

@info("Load database...")
dbpath = homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B"

# Loading the database directly
# (this database can be generated with first_example_butane.jl)
db = LsqDB(dbpath)
summary(db)

##
@info("Fit Butane Database basis...")
#Weights on the observations
obsweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)

# changing configuration names
for c in db.configs
    c.configtype = "test"
end
#configuration weights
configweights = Dict("test"  => 1.0)

#Define one-body term
oneB = MOneBody(Dict("E0" => Dict(:H => -13.5203101677, :C => -1027.28368153)))

IP, lsqinfo = lsqfit( db;
                       # E0 = 1.,
                       obsweights=obsweights,
                       configweights=configweights,
                       solver = (:rrqr, 1e-16), # solver = (:qr,), - to choose
                       combineIP = NBodyIP,
                       # Ibasis = Ibasis #for choosing indices
                       Vref = oneB #reference potential
                       )

#Look at errors
errs = lsqinfo["errors"]
rmse_table(rmse(errs)...)

# Save the potential together with the fit info
save_ip(homedir() * "/.julia/dev/NBIPsMulti/data/Butane_4B_IP.json", IP, lsqinfo)


## Adding a 2B repulsion
using NBIPsMulti.RepulsionM: RepulsiveCoreM

using NBodyIPs.Repulsion: RepulsiveCore

# Extract 2-body components
V2_1 = IP.components[2]
V2_2 = IP.components[3]
V2_3 = IP.components[4]

# Adding 2-body repulsion
r0 = 3.
V2rep_1 = RepulsiveCoreM(V2_1, 0.89*r0, -1.)
V2rep_2 = RepulsiveCoreM(V2_2, 0.89*r0, -1.)
V2rep_3 = RepulsiveCoreM(V2_3, 0.89*r0, -1.)

# Redefine the potential with repulsive core
IP_rep = deepcopy(IP)

IP_rep.components[2] = V2rep_1
IP_rep.components[3] = V2rep_2
IP_rep.components[4] = V2rep_3
save_ip("pot_with_rep_core.json", IP_rep, lsqinfo) #save potential

# Show first configuration and its energy
@show db.configs[1].at
@show energy(IP_rep, db.configs[1].at) #with repulsive core
@show energy(IP, db.configs[1].at) #without repulsive core

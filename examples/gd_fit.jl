info("Load libraries...")

include("../src/NBIPsMulti.jl")

using NBodyIPs, NBodyIPFitting, JuLIP, NBIPsMulti
include(homedir() * "/Gits/NBIPsMulti/src/Butane.jl")
E0 = Butane.get_E0()

info("Load database...")
# run(`ls $(homedir() * "/scratch/nbodyips/")`)
dbpath = homedir() * "/Documents/Butane_3B"
# dbpath = homedir() * "/scratch/nbodyips/Si_4B_BA_long"
# dbpath = homedir() * "/scratch/nbodyips/Si_5B_BA"
db = LsqDB(dbpath)
db
# summary(db)
##
info("Fit Butane Database basis...")

dataweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)
configweights = Dict(""  => 1.0)
# Ibasis = [1:12; 21:length(db.basis)]

IP, info = lsqfit( db; E0 = E0,
                       dataweights=dataweights,
                       configweights=configweights,
                       # Ibasis = Ibasis
                       )
info

IP


table_absolute(info["errors"])
table_relative(info["errors"])

at = bulk(:Si) * 3

(IP.components[1].E0)*54

[energy(Vn, at)  for Vn in IP.components]

@show energy(IP, at)
# IPf = fast(IP)
# energy(IPf, at)
NBodyIPs.save_ip("SiPIP_4BBA_short.json", IP, info
   )


IP, info = NBodyIPs.load_ip("SiPIP_4BBA_short.json")
info

## -------------- WITH REGULARISATION -------------

r0 = rnn(:Si)
rcuts = unique(cutoff.(db.basis))
rcut4, rcut3, rcut2 = (sort(rcuts)...)
reg = [ BAReg(2, 0.5*r0, 0.85*r0, creg=3.0),
        BAReg(3, 0.7*r0, 0.85*r0, creg=0.3),
        BAReg(4, 0.7*r0, 0.85*r0, creg=0.3),
        BAReg(2, 0.85*r0, rcut2, creg=0.03),
        BAReg(3, 0.85*r0, rcut3, creg=0.03),
        BAReg(4, 0.85*r0, rcut4, creg=0.03),
     ]

IPreg, errs = lsqfit( db; E0 = E0,
                       dataweights   = dataweights,
                       configweights = configweights,
                       regularisers  = reg,
                       Ibasis = Ibasis
  )

table_absolute(errs["errors"])
table_relative(errs["errors"])

energy(IPreg, at)
NBodyIPs.save_ip("SiPIP_R4BBA_short.json", IPreg)

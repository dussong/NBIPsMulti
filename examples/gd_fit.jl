info("Load libraries...")

include("../src/NBIPsMulti.jl")

using NBodyIPs, NBodyIPFitting, JuLIP, NBIPsMulti
include(homedir() * "/Gits/NBIPsMulti/src/Butane.jl")
E0 = Butane.get_E0()

info("Load database...")
dbpath = homedir() * "/Gits/NBIPsMulti/data/Butane_3B"
db = LsqDB(dbpath)
db
# summary(db)
##
info("Fit Butane Database basis...")

dataweights = Dict("E" => 10.0, "F" => 1.0, "V" => 1.0)
configweights = Dict(""  => 1.0)

IP, info = lsqfit( db; E0 = E0,
                       dataweights=dataweights,
                       configweights=configweights,
                       # Ibasis = Ibasis
                       )
info

IP

table_absolute(info["errors"])
table_relative(info["errors"])

NBodyIPs.save_ip("Butane_3B.json", IP, info)


# IP, info = NBodyIPs.load_ip("SiPIP_4BBA_short.json")
# info

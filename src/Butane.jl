# read and load the Butane data
module Butane

using IPFitting, FileIO
using IPFitting: configtype

export read_Butane

get_E0() = -5.817622899211898

filename() = homedir() * "/.julia/dev/NBIPsMulti/data/1500K_TB_butane"

function loaddb(;  include=nothing, kwargs...)
   # default - read from jld2.
   fname = filename() * ".jld2"
   println("Reading data from $fname")
   data = load(fname, "data")
   if include != nothing
      data = data[ [ configname(d) in include for d in data] ]
   end
   return data
end


function load_xyz(; md=true, hessians=true, perturb=true, kwargs...)
   # the main database file
   data = Dat[]
   if md
      fname = filename() * ".xyz"
      data = IPFitting.Data.read_xyz(fname; kwargs...)
   end
   return data
end

function convert_db(; fname = filename() * ".jld2", kwargs...)
   data = load_xyz(; kwargs...)
   save(fname, "data", data)
end

end

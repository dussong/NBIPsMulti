
module Butane

using NBodyIPFitting, FileIO
using NBodyIPFitting: configtype

export read_Butane

get_E0() = -5.817622899211898

# filename() = homedir() * "/Gits/NBIPsMulti/data/Butane_short"
filename() = homedir() * "/Gits/NBIPsMulti/data/1500K_TB_butane"

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
      data = NBodyIPFitting.Data.read_xyz(fname; kwargs...)
   end
   return data
end

function convert_db(; fname = filename() * ".jld2", kwargs...)
   data = load_xyz(; kwargs...)
   save(fname, "data", data)
end

end

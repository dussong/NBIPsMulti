
module Butane

using NBodyIPFitting, FileIO
using NBodyIPFitting: configtype

export read_Butane

get_E0() = -5.817622899211898

filename() = homedir() * "/Gits/NBIPsMulti/data/Butane_short"

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
   # if perturb
   #    data_h_bcc = NBodyIPFitting.Data.read_xyz("Ti_bcc_unit_cell_disp_atom.xyz"; verbose=false)
   #    for d in data_h_bcc
   #       d.configtype = "perturb_bcc" * d.configtype
   #    end
   #    data_h_hcp = NBodyIPFitting.Data.read_xyz("Ti_hcp_unit_cell_disp_atom.xyz"; verbose=false)
   #    for d in data_h_hcp
   #       d.configtype = "perturb_hcp" * d.configtype
   #    end
   #    data_h_omega = NBodyIPFitting.Data.read_xyz("Ti_omega_unit_cell_disp_atom.xyz"; verbose=false)
   #    for d in data_h_omega
   #       d.configtype = "perturb_omega" * d.configtype
   #    end
   #    data = [data; data_h_bcc; data_h_hcp; data_h_omega]
   # end
   # if hessians
   #    data_h_bcc = NBodyIPFitting.Data.read_xyz("Ti_bcc_hess_train.xyz"; verbose=false)[1:6]
   #    for d in data_h_bcc
   #       d.configtype = "hess_bcc" * d.configtype
   #    end
   #    data_h_hcp = NBodyIPFitting.Data.read_xyz("Ti_hcp_hess_train_v3_relaxed_quippy.xyz"; verbose=false)[1:6]
   #    for d in data_h_hcp
   #       d.configtype = "hess_hcp" * d.configtype
   #    end
   #    data_h_omega = NBodyIPFitting.Data.read_xyz("Ti_omega_hess_train.xyz"; verbose=false)[1:6]
   #    for d in data_h_omega
   #       d.configtype = "hess_omega" * d.configtype
   #    end
   #    data = [data; data_h_bcc; data_h_hcp; data_h_omega]
   # end
   return data
end

function convert_db(; fname = filename() * ".jld2", kwargs...)
   data = load_xyz(; kwargs...)
   save(fname, "data", data)
end

end

if length(ARGS) > 0 && ARGS[1] in ["-convert", "convert"]
   Butane.convert_db()
end

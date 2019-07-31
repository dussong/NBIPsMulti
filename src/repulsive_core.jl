
module RepulsionM

import NBodyIPs: NBodyFunction, descriptor, evaluate_I, evaluate_I_ed,
                 inv_transform, fcut, fast

import JuLIP: AbstractCalculator, energy, forces, virial, decode_dict
import JuLIP.Potentials: @pot, evaluate, evaluate_d, PairPotential, @D, cutoff,
                         @analytic
import Base: Dict, convert

using NBIPsMulti: NBodyFunctionM

struct RepulsiveCoreM{TV1, TV2, DT, SP} <: NBodyFunctionM{2, DT, SP}
   Vout::TV1         # the outer pair potential
   Vin::TV2          # the inner (repulsive) pair potential
   ri::Float64       # the interface between the two
   e0::Float64
   D::DT             # shouldn't be needed, but can't get it to work otherwise...
   Sp::Vector{Int}            #encodes the species
   Sp_type::Val{SP}
end

@pot RepulsiveCoreM

cutoff(V::RepulsiveCoreM) = cutoff(V.Vout)
descriptor(V::RepulsiveCoreM) = V.D

fast(V::RepulsiveCoreM) = RepulsiveCoreM(fast(V.Vout), V.Vin, V.ri, V.e0, V.D, V.Sp, V.Sp_type)

function evaluate_I(V::RepulsiveCoreM, II)
   I1, I2 = II
   @assert length(I1) == 1
   r = inv_transform(V.D.transform, I1[1])
   if r > V.ri
      return evaluate_I(V.Vout, II)
   end
   return V.Vin(r) / fcut(V.Vout.D, r)
end

function evaluate_I_ed(V::RepulsiveCoreM, II)
   I1, I2, dI1, dI2 = II
   @assert length(I1) == 1
   r = inv_transform(V.D.transform, I1[1])
   if r > V.ri
      return evaluate_I_ed(V.Vout, II)
   end
   eV = V.Vin(r)
   dV = @D V.Vin(r)
   return eV, dV
end


function RepulsiveCoreM(Vout::NBodyFunctionM{2}, ri, e0=0.0)
   v = Vout(ri)
   dv = @D Vout(ri)
   if dv >= 0.0
      @error("The slope `Vout'(ri)` is required to be negative")
   end
   if dv > -1.0
      @warn("""The slope `Vout'(ri) = $dv` may not be steep enough to attach a
               repulsive core. Proceed at your own risk.""")
   end
   if v-e0 <= 0.0
      @error("it is required that `Vout(ri) > 0`.")
   end
   if v-e0 <= 1.0
      @warn("""Ideally the repulsive core should not be attached at small
               values of `Vout(ri) = $v`. Proceed at your own risk.""")
   end
   # e0 + B e^{-A (r/ri-1)} * ri/r
   #    => e0 + B = Vout(ri) => = Vout(ri) - e0 = v - e0
   # dv = - A*B/ri e^{-A (r/ri-1)} * ri/r - B*ri*e^{...} / r^2
   #    = - A/ri * (v - 1/ri * (v = - (1+A)/ri * (v-e0)
   #    => -(1+A)/ri * (v-e0) = dv
   #    => 1+A = - ri dv / (v-e0)
   #    => A = -1 - ri dv / (v-e0)
   Vin = let A = -1 - ri * dv / (v-e0), B = v-e0, e0=e0, ri = ri
      @show A
      @show B
      @analytic r -> e0 + B * exp( - A * (r/ri-1) ) * ri/r
   end
   # @show ri
   # @show Vout(ri), (@D Vout(ri))
   # @show Vin(ri), (@D Vin(ri))
   # construct the piecewise potential
   return RepulsiveCoreM(Vout, Vin, ri, e0,
                         descriptor(Vout),Vout.Sp,Vout.Sp_type)
end

# ----------------------------------------------------
#  File IO
# ----------------------------------------------------

Dict(V::RepulsiveCoreM) = Dict("__id__" => "NBodyIPs_RepulsiveCoreM",
                              "Vout" => Dict(V.Vout),
                              "e0" => V.e0,
                              "ri" => V.ri)

function RepulsiveCoreM(D::Dict)
   if haskey(D, "e0")
      return RepulsiveCoreM(decode_dict(D["Vout"]), D["ri"], D["e0"])
   else
      return RepulsiveCoreM(decode_dict(D["Vout"]), D["ri"])
   end
end

convert(::Val{:NBodyIPs_RepulsiveCoreM}, D::Dict) = RepulsiveCoreM(D)

end

# EXAMPLE CODE => TURN THiS INTO A TEST!
# ## try out the repulsive potential
# Vfit = IPall.components[2]
# Vrep = NBodyIPs.Repulsion.RepulsiveCore(Vfit, 2.1)
#
# rp = range(0.3*r0, 9.0, length=500)
# plot(rp, Vfit.(rp), lw=2, label ="fit")
# plot!(rp, Vrep.Vin.(rp), lw=2, label ="inner", ylims = [-2.0, 10.0] )
# plot!(rp, Vrep.(rp), label = "combined")

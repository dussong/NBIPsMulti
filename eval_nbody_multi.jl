using JuLIP, NeighbourLists
using JuLIP: AbstractCalculator
using JuLIP.Potentials: @pot

using NBodyIPs: NBodyFunction, bapolys, eval_site_nbody!, evaluate

import JuLIP: site_energies, energy, forces

# TODO: implement site_energies, energy, forces, virial


# implement site_energies for given species
# Pair potential
function site_energies(V::NBodyFunction{2}, at::Atoms{T},
                                       species::Tuple{Int,Int}) where {T}
   Z = atomic_numbers(at)
   Es = zeros(T, length(at))
   for (i, j, r, R) in sites(at, cutoff(V))
      for k = 1:length(j)
         atnb = sort([Z[i],Z[j][k]])
         if (atnb[1] == species[1])&(atnb[2] == species[2])
            Es[i] += V(r[k])
         end
      end
   end
   return Es
end

function site_energies(V::NBodyFunction{2}, at::Atoms{T},
                                       species::Tuple{Symbol,Symbol}) where {T}
   sp = atomic_number.(species)
   return site_energies(V, at, sp)
end

r0 = 2.5
V = bapolys(2, "($r0/r)^4", "(:cos, 3.6, 4.8)", 2)
Vcucu = V[2]
at = bulk(:Cu, cubic=true)*2
species = (:Cu,:Zn)

atomic_number(:Cu)

Vcucu(3.)

site_energies(Vcucu, at, (29,29))
site_energies(Vcucu, at, (:Cu,:Cu))


Vcucu(x) = x
Vcuzr(x) = 2*x
Vzrzr(x) = 3*x
V2s = Dict(
            (:Cu, :Cu) => Vcucu,
            (:Cu, :Zr) => Vcuzr,
            (:Zr, :Zr) => Vcuzr
         )
at = bulk(:Cu, cubic=true)*2
rcut = 2.7

Z = atomic_numbers(at)
species = (:Cu,:Zn)
for (i, j, r, R) in sites(atu, rcut)
   Zi =
   @show i
   @show j
   @show r
   @show size(R)
end

sites(atu,rcut)
site_energies(V2s, atu,rcut)
energy(V2s,atu,rcut)

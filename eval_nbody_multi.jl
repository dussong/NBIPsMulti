using JuLIP, NeighbourLists
using JuLIP: AbstractCalculator
using JuLIP.Potentials: @pot

using NBodyIPs: NBodyFunction, bapolys, eval_site_nbody!, evaluate, eval_site_nbody!, evaluate_d!

import JuLIP: site_energies, energy, forces

using NeighbourLists: nbodies,
                      maptosites!,
                      maptosites_d!,
                      virial!,
                      max_neigs,
                      sites

# TODO: implement site_energies, energy, forces, virial

struct NBodyFunctionM
   d::Dict{Tuple{Symbol,Symbol},NBodyFunction{2}}
end



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

# using symbols
function site_energies(V::NBodyFunction{2}, at::Atoms{T},
                                       species::Tuple{Symbol,Symbol}) where {T}
   sp = atomic_number.(species)
   return site_energies(V, at, sp)
end

energy(V::NBodyFunction, at::Atoms,species::Tuple{Symbol,Symbol}) = sum_kbn(site_energies(V, at,species))

# Implementation of the forces
function forces(V::NBodyFunction{2}, at::Atoms{T},sp::Tuple{Int,Int}) where { T}
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   F = zeros(JVec{T}, length(at))
   dVsite = zeros(JVec{T}, maxneigs)
   for (i, j, r, R) in sites(nlist)
      fill!(dVsite, zero(JVec{T}))
      eval_site_nbody!(
            Val(2), R, cutoff(V),
            (out, R, J, temp) -> evaluate_d!(out, V, R, J),
            dVsite, nothing )   # dVsite == out, nothing == temp
      # write site energy gradient into forces
      for n = 1:length(j)
         F[j[n]] -= dVsite[n]
         F[i] += dVsite[n]
      end
   end
   return F
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
energy(Vcucu,at,(:Cu,:Cu))
forces(Vcucu,at,(29,29))

# Vcucu(x) = x
# Vcuzr(x) = 2*x
# Vzrzr(x) = 3*x
# V2s = Dict(
#             (:Cu, :Cu) => Vcucu,
#             (:Cu, :Zr) => Vcuzr,
#             (:Zr, :Zr) => Vcuzr
#          )
# at = bulk(:Cu, cubic=true)*2
# rcut = 2.7
#
# Z = atomic_numbers(at)
# species = (:Cu,:Zn)
# for (i, j, r, R) in sites(atu, rcut)
#    Zi =
#    @show i
#    @show j
#    @show r
#    @show size(R)
# end
#
# sites(atu,rcut)
# site_energies(V2s, atu,rcut)
# energy(V2s,atu,rcut)

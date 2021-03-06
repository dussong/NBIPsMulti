using JuLIP, NeighbourLists
using JuLIP: AbstractCalculator
using JuLIP.Potentials: @pot
using StaticArrays
using KahanSummation

using NBodyIPs: NBodyFunction, bapolys, eval_site_nbody!, evaluate, eval_site_nbody!, evaluate_d!, NBSiteDescriptor, _get_loop_ex, _get_Jvec_ex, descriptor, ricoords, skip_simplex, fcut, invariants, evaluate_I, fcut_d, invariants_ed, evaluate_I_ed, gradri2gradR!, evaluate_many!

import JuLIP: site_energies,
              energy,
              forces,
              virial

import NBodyIPs: evaluate,
                 eval_site_nbody!,
                 evaluate_d!,
                 evaluate_many!,
                 evaluate_many_d!

using NeighbourLists: nbodies,
                      maptosites!,
                      maptosites_d!,
                      max_neigs,
                      sites


# skip the simplex if not the right species
function skip_simplex_species!(Spi,Spj,Species,J,tmp)
   tmp[1] = Spi
   @assert(length(tmp) == length(J)+1)
   for i=1:length(J)
      tmp[i+1] = Spj[J[i]]
   end
   sort!(tmp, alg=InsertionSort)
   return tmp != Species
end

# skip the simplex if not the right species
function skip_simplex_species_many!(Spi,Spj,Species,J,tmp)
   tmp[1] = Spi
   @assert(length(tmp) == length(J)+1)
   for i=1:length(J)
      tmp[i+1] = Spj[J[i]]
   end
   sort!(tmp)
   return [tmp == Species[k] for k=1:length(Species)]
end


# 2-Body
function skip_simplex_species_order! end

skip_simplex_species_order!(desc::MultiDesc,
                            Spi,Spj,Species,J,::Val{:AA}) = false

# 3-Body
skip_simplex_species_order!(desc::MultiDesc,
                            Spi,Spj,Species,
                            J,::Val{:AAA}) = false

skip_simplex_species_order!(desc::MultiDesc,
                          Spi,Spj,Species,
                          J,::Val{:AAAba}) = false

function skip_simplex_species_order!(desc::MultiDesc,
                                     Spi,Spj,Species,J,::Val{:AAB})
   return Spj[J[1]] != Spj[J[2]]
end

function skip_simplex_species_order!(desc::MultiDesc,
                                     Spi,Spj,Species,J,::Val{:AABba})
   return Spj[J[1]] != Spj[J[2]]
end

function skip_simplex_species_order!(desc::MultiDesc,
                                     Spi,Spj,Species,J,::Val{:ABC})
   if Spi == Species[1]
      if Spj[J[1]] == Species[2]
         return false
      elseif Spj[J[1]] == Species[3]
         J = [J[2],J[1]]
         return false
      else
         error("I shouldnt be here: skip_simplex_species_order")
      end
   else
      return true
   end
end

skip_simplex_species_order!(desc::MultiDesc,
                           Spi,Spj,Species,J,::Val{:ABCba}) = false

# 4-Body

skip_simplex_species_order!(desc::MultiDesc,
                            Spi,Spj,Species,
                            J,::Val{:AAAA}) = false

skip_simplex_species_order!(desc::MultiDesc,
                          Spi,Spj,Species,
                          J,::Val{:AAAAba}) = false

function skip_simplex_species_order!(desc::MultiDesc,
                                   Spi,Spj,Species,J,::Val{:AAAB})
   return (Spj[J[1]] != Spj[J[2]])||(Spj[J[1]] != Spj[J[3]])||(
         Spj[J[2]] != Spj[J[3]])
end

function skip_simplex_species_order!(desc::MultiDesc,
                                   Spi,Spj,Species,J,::Val{:AAABba})
   return (Spj[J[1]] != Spj[J[2]])||(Spj[J[1]] != Spj[J[3]])||(
         Spj[J[2]] != Spj[J[3]])
end

function skip_simplex_species_order!(desc::MultiDesc,
                                   Spi,Spj,Species,J,::Val{:AABB})
   if Spi == Species[1]
      # ind = sortperm([Spj[J[k]] for k=1:3])
      # J[1],J[2],J[3] = J[ind[1]],J[ind[2]],J[ind[3]]
      # return false
      if (Spj[J[1]] == Species[2])
         return false
      elseif (Spj[J[2]] == Species[2])
         J = [J[2],J[1],J[3]]
         # J[2], J[1] = J[1], J[2]
         return false
      elseif (Spj[J[3]] == Species[2])
         J = [J[3], J[2], J[1]]
         # J[3], J[1] = J[1], J[3]
         return false
      end
   else
      return true
   end
end

function skip_simplex_species_order!(desc::MultiDesc,
                                   Spi,Spj,Species,J,::Val{:AABBba})
   if Spi == Species[1]
      # ind = sortperm([Spj[J[k]] for k=1:3])
      # J[1],J[2],J[3] = J[ind[1]],J[ind[2]],J[ind[3]]
      # return false
      if (Spj[J[1]] == Species[2])
         return false
      elseif (Spj[J[2]] == Species[2])
         J = [J[2],J[1],J[3]]
         # J[2], J[1] = J[1], J[2]
         return false
      elseif (Spj[J[3]] == Species[2])
         J = [J[3], J[2], J[1]]
         # J[3], J[1] = J[1], J[3]
         return false
      end
   else
      return true
   end
end

function skip_simplex_species_order!(desc::MultiDesc,
                                   Spi,Spj,Species,J,::Val{:AABC})
   if Spi == Species[1]
      ind = sortperm([Spj[J[k]] for k=1:3])
      J = [J[ind[1]],J[ind[2]],J[ind[3]]]
      return false
   else
      return true
   end
end

function skip_simplex_species_order!(desc::MultiDesc,
                                   Spi,Spj,Species,J,::Val{:AABCba})
   if Spi == Species[1]
      ind = sortperm([Spj[J[k]] for k=1:3])
      J = [J[ind[1]],J[ind[2]],J[ind[3]]]
      return false
   else
      return true
   end
end



function skip_simplex_species_order!(desc::MultiDesc,
                                     Spi,Spj,Species,J,::Val{:ABCD})
   if Spi == Species[1]
      if (Spj[J[1]] == Species[2])&&(Spj[J[2]] == Species[3])
         return false
      elseif (Spj[J[1]] == Species[2])&&(Spj[J[2]] == Species[4])
         J = [J[1],J[3],J[2]]
         # J[2], J[3] = J[3], J[2]
         return false
      elseif (Spj[J[1]] == Species[3])&&(Spj[J[2]] == Species[2])
         J = [J[2], J[1], J[3]]
         # J[2], J[1], J[3] = J[1], J[2], J[3]
         return false
      elseif (Spj[J[1]] == Species[3])&&(Spj[J[2]] == Species[4])
         J = [J[3], J[1], J[2]]
         # J[2], J[3], J[1] = J[1], J[2], J[3]
         return false
      elseif (Spj[J[1]] == Species[4])&&(Spj[J[2]] == Species[2])
         J = [J[2],J[3],J[1]]
         # J[3], J[1], J[2] = J[1], J[2], J[3]
         return false
      elseif (Spj[J[1]] == Species[4])&&(Spj[J[2]] == Species[3])
         J = [J[3],J[2],J[1]]
         # J[3], J[1] = J[1], J[3]
         return false
      else
         error("I shouldnt be here: skip_simplex_species_order")
      end
   else
      return true
   end
end

skip_simplex_species_order!(desc::MultiDesc,
                           Spi,Spj,Species,J,::Val{:ABCDba}) = false

skip_simplex_species_order!(desc::MultiDesc,
                            Spi,Spj,Species,J) =
                            skip_simplex_species_order!(desc,
                                                        Spi,Spj,Species,J,
                                                        desc.sp_type)





# sort by species for 3B:
# function sort_by_species(Spi,Spj,J::)
#    Sp = [Spi]
#    for i=1:length(J)
#       push!(Sp,Spj[J[i]])
#    end
#    if Sp[1] == Sp[2]
#       return J
#    elseif
# end


@generated function eval_site_nbody!( ::Val{N},
                                      Rs::AbstractVector{JVec{T}},
                                      rcut::T,
                                      reducefun,
                                      out,
                                      temp,
                                      Spi,Spj,Species) where {N, T}
   code = Expr[]
   # initialise the output
   push!(code, :( nR = length(Rs)  ))

   # generate the multi-for-loop
   ex_loop = _get_loop_ex(N)

   # inside the loop
   # ---------------
   code_inner = Expr[]
   # collect the indices into a vector
   push!(code_inner,      _get_Jvec_ex(N) )

   # now call `V` with the simplex-corner vectors and "add" this to the site energy
   push!(code_inner, :(   out = reducefun(out, Rs, J, temp,Spi,Spj,Species) ))

   # put code_inner into the loop expression
   ex_loop.args[2] = Expr(:block, code_inner...)

   # now append the loop to the main code
   push!(code, ex_loop)

   quote
      @inbounds $(Expr(:block, code...))
      return out
   end
end


function site_energies(V::NBodyFunctionM{N},at::Atoms{T}, Species::Vector{Int}) where {N, T}
   # ,nlist = neighbourlist(at, cutoff(V))
   sort!(Species)
   Es = zeros(T, length(at))
   Z = atomic_numbers(at)
   tmp = zeros(Int,N)
   Spi = 0
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   Spj = zeros(Int,maxneigs)
   for (i, j, r, R) in sites(at, cutoff(V))
      Spi = Z[i]
      # if (Spi in Species)
         for n=1:length(j)
            Spj[n] = Z[j[n]]
         end
         Es[i] = eval_site_nbody!(Val(N), R, cutoff(V),
                                  ((out, R, J, temp,Spi,Spj,Species) -> out + evaluate(V, descriptor(V), R, J,Spi,Spj,Species,tmp)), zero(T), nothing, Spi,Spj,Species)
      # else
      #    Es[i] = 0.
      # end
   end
   return Es
end

energy(V::NBodyFunctionM, at::Atoms, Species::Vector{Int}) = sum_kbn(site_energies(V, at, Species))
# ,nlist = neighbourlist(at, cutoff(V))


function energy(V::NBodyFunctionM, at::Atoms)
   # @show V
   # @show V.Sp
   # @show V.Sp_type
   # @show at
   return energy(V,at,V.Sp)
end



function forces(V::NBodyFunctionM{N}, at::Atoms{T},Species::Vector{Int}) where {N, T}
   sort!(Species)
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   F = zeros(JVec{T}, length(at))
   dVsite = zeros(JVec{T}, maxneigs)
   Z = atomic_numbers(at)
   tmp = zeros(Int,N)
   Spi = 0
   Spj = zeros(Int,maxneigs)
   for (i, j, r, R) in sites(nlist)
      Spi = Z[i]
      for n = 1:length(j)
         Spj[n] = Z[j[n]]
      end
      fill!(dVsite, zero(JVec{T}))
      eval_site_nbody!(Val(N), R, cutoff(V),
                               ((out, R, J, temp,Spi,Spj,Species) ->  evaluate_d!(out, V, descriptor(V), R, J,Spi,Spj,Species,tmp)), dVsite, nothing, Spi,Spj,Species)
      # write site energy gradient into forces
      for n = 1:length(j)
         F[j[n]] -= dVsite[n]
         F[i] += dVsite[n]
      end
   end
   return F
end

function forces(V::NBodyFunctionM, at::Atoms)
   return forces(V,at,V.Sp)
end

# ========================= assembly support for LSQ system ====================

# For assembling the LSQ system efficiently we need a way to evaluate all basis
# functions of the same body-order at the same time. Otherwise we would be
# re-computing the invariants many many times, which is very expensive.
# To achieve this we just wrap all basis functions of a body-order into
# a new type `NBBasis` which evaluates to a long vector
#
# at the moment, it seems we need to hard-code this to the Polys
# sub-module, but it would be good if this can be fixed, so we keep this
# "interface" here.


function energy(B::AbstractVector{TB}, at::Atoms{T}
                ) where {TB <: NBodyFunctionM{N}, T} where {N}
   # TODO: assert that all B[j] have the same invariants
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   E = zeros(T, length(B))
   Z = atomic_numbers(at)
   Species = [sort(B[i].Sp) for i=1:length(B)]
   tmp = zeros(Int,N)
   for (i, j, r, R) in sites(nlist)
      Spi = Z[i]
      Spj = Z[j]
      # evaluate all the site energies at the same time
      # for each simplex, write the nB energies into temp
      # then add them to E, which is just passed through all the
      # various loops, so no need to update it here again
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp,Spi,Spj,Species) -> evaluate_many!(out, B, R, J, Spi, Spj,Species, tmp),
                       E, nothing, Spi,Spj,Species)
   end
   return E
end



function forces(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunctionM{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nB = length(B)
   # forces
   F =      [ zeros(JVec{T}, length(at)) for n = 1:nB ]
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:nB ]

   Species = [sort(B[i].Sp) for i=1:length(B)]
   Z = atomic_numbers(at)
   tmp = zeros(Int,N)

   for (i, j, r, R) in sites(nlist)
      Spi = Z[i]
      Spj = Z[j]
      # clear dVsite
      for n = 1:nB; fill!(dVsite[n], zero(JVec{T})); end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp, Spi, Spj, Species) -> evaluate_many_d!(out, B, R, J, Spi, Spj, Species, tmp),
                       dVsite, nothing, Spi,Spj,Species)
      # write it into the force vectors
      for ib = 1:nB, n = 1:length(j)
         F[ib][j[n]] -= dVsite[ib][n]
         F[ib][i] += dVsite[ib][n]
      end
   end
   return F
end


function virial(B::AbstractVector{TB}, at::Atoms{T}
                ) where {TB <: NBodyFunctionM{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nB = length(B)
   # virials (main output)
   S = fill((@SMatrix zeros(3,3)), nB)
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:nB ]

   Species = [sort(B[i].Sp) for i=1:length(B)]
   Z = atomic_numbers(at)
   tmp = zeros(Int,N)

   for (i, j, r, R) in sites(nlist)
      Spi = Z[i]
      Spj = Z[j]
      # clear dVsite
      for n = 1:nB; dVsite[n] .*= 0.0; end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp, Spi, Spj, Species) -> evaluate_many_d!(out, B, R, J, Spi, Spj, Species, tmp),
                       dVsite, nothing,Spi,Spj,Species)
      # update the virials
      for iB = 1:nB
         S[iB] += JuLIP.Potentials.site_virial(dVsite[iB], R)
      end
   end
   return S
end

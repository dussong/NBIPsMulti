using JuLIP, NeighbourLists
using JuLIP: AbstractCalculator
using JuLIP.Potentials: @pot
using StaticArrays
using BenchmarkTools

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
                      virial!,
                      max_neigs,
                      sites


# skip the simplex if not the right species
function skip_simplex_species(Spi,Spj,Species,J)
   Sp = [Spi]
   for i=1:length(J)
      push!(Sp,Spj[J[i]])
   end
   return sort(Sp) != sort(Species)
end

# skip the simplex if not the right species
function skip_simplex_species_many(Spi,Spj,Species,J)
   Sp = [Spi]
   for i=1:length(J)
      push!(Sp,Spj[J[i]])
   end
   return [sort(Sp) == sort(Species[k]) for k=1:length(Species)]
end

# sort by species for 3B:
function sort_by_species(Spi,Spj,J::)
   Sp = [Spi]
   for i=1:length(J)
      push!(Sp,Spj[J[i]])
   end
   if Sp[1] == Sp[2]
      return J
   elseif
end


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


function evaluate(V::NBodyFunction{N},
                  desc::NBSiteDescriptor,
                  Rs::AbstractVector{JVec{T}},
                  J::SVector{K, Int},
                  Spi::Int,Spj::Vector{Int},Species::Vector{Int}) where {N, T, K}
   # check species
   skip_simplex_species(Spi,Spj,Species,J) && return zero(T)
   evaluate(V,desc,Rs,J)
end

evaluate(V::NBodyFunction,Rs::AbstractVector{JVec{T}},J::SVector{K, Int},Spi,Spj,Species) where {T,K} = evaluate(V,descriptor(V),Rs,J,Spi,Spj,Species)

function evaluate_d!(dVsite,
                     V::NBodyFunction{N},
                     desc::NBSiteDescriptor,
                     Rs::AbstractVector{JVec{T}},
                     J,
                     Spi::Int,Spj::Vector{Int},Species::Vector{Int}) where {N,T}
   # check species
   skip_simplex_species(Spi,Spj,Species,J) && return dVsite
   evaluate_d!(dVsite, V, desc, Rs, J)
end








function site_energies(V::NBodyFunctionM{N}, at::Atoms{T},Species::Vector{Int}) where {N, T}
   Es = zeros(T, length(at))
   Z = atomic_numbers(at)
   for (i, j, r, R) in sites(at, cutoff(V))
      Spi = Z[i]
      Spj = Z[j]
      Es[i] = eval_site_nbody!(Val(N), R, cutoff(V),
                               ((out, R, J, temp,Spi,Spj,Species) -> out + evaluate(V, descriptor(V), R, J,Spi,Spj,Species)), zero(T), nothing, Spi,Spj,Species)
   end
   return Es
end

energy(V::NBodyFunctionM, at::Atoms, Species::Vector{Int}) = sum_kbn(site_energies(V, at, Species))


function energy(V::NBodyFunctionM, at::Atoms)
   return energy(V,at,V.Sp)
end



function forces(V::NBodyFunctionM{N}, at::Atoms{T},Species::Vector{Int}) where {N, T}
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   F = zeros(JVec{T}, length(at))
   dVsite = zeros(JVec{T}, maxneigs)
   Z = atomic_numbers(at)
   for (i, j, r, R) in sites(nlist)
      Spi = Z[i]
      Spj = Z[j]
      fill!(dVsite, zero(JVec{T}))
      eval_site_nbody!(Val(N), R, cutoff(V),
                               ((out, R, J, temp,Spi,Spj,Species) ->  evaluate_d!(out, V, descriptor(V), R, J,Spi,Spj,Species)), dVsite, nothing, Spi,Spj,Species)
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


function evaluate_many!(Es,
                        B::AbstractVector{TB},
                        desc::NBSiteDescriptor,
                        Rs, J, Spi,Spj,Species)  where {TB <: NBodyFunctionM{N}} where {N}
   ind = find(skip_simplex_species_many(Spi,Spj,Species,J))
   Es[ind] = evaluate_many!(Es[ind],B[ind],desc,Rs,J)
   return Es
end

evaluate_many!(out, B, Rs, J, Spi, Spj, Species) =
      evaluate_many!(out, B, descriptor(B[1]), Rs, J, Spi, Spj, Species)



function energy(B::AbstractVector{TB}, at::Atoms{T}
                ) where {TB <: NBodyFunctionM{N}, T} where {N}
   # TODO: assert that all B[j] have the same invariants
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   E = zeros(T, length(B))
   Z = atomic_numbers(at)
   Species = [B[i].Sp for i=1:length(B)]
   for (i, j, r, R) in sites(nlist)
      Spi = Z[i]
      Spj = Z[j]
      # evaluate all the site energies at the same time
      # for each simplex, write the nB energies into temp
      # then add them to E, which is just passed through all the
      # various loops, so no need to update it here again
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp,Spi,Spj,Species) -> evaluate_many!(out, B, R, J, Spi, Spj,Species),
                       E, nothing, Spi,Spj,Species)
   end
   return E
end



function evaluate_many_d!(dVsite::AbstractVector,
                          B::AbstractVector{TB},
                          desc::NBSiteDescriptor,
                          Rs,
                          J, Spi,Spj,Species)  where {TB <: NBodyFunctionM{N}} where {N}
   ind = find(skip_simplex_species_many(Spi,Spj,Species,J))
   dVsite[ind] = evaluate_many_d!(dVsite[ind],B[ind],desc,Rs,J)
   return dVsite
end

evaluate_many_d!(out, B, Rs, J, Spi, Spj, Species) =
      evaluate_many_d!(out, B, descriptor(B[1]), Rs, J, Spi,Spj,Species)


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

   Species = [B[i].Sp for i=1:length(B)]
   Z = atomic_numbers(at)

   for (i, j, r, R) in sites(nlist)
      Spi = Z[i]
      Spj = Z[j]
      # clear dVsite
      for n = 1:nB; fill!(dVsite[n], zero(JVec{T})); end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp, Spi, Spj, Species) -> evaluate_many_d!(out, B, R, J, Spi, Spj, Species),
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

   Species = [B[i].Sp for i=1:length(B)]
   Z = atomic_numbers(at)

   for (i, j, r, R) in sites(nlist)
      Spi = Z[i]
      Spj = Z[j]
      # clear dVsite
      for n = 1:nB; dVsite[n] .*= 0.0; end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp, Spi, Spj, Species) -> evaluate_many_d!(out, B, R, J, Spi, Spj, Species),
                       dVsite, nothing,Spi,Spj,Species)
      # update the virials
      for iB = 1:nB
         S[iB] += JuLIP.Potentials.site_virial(dVsite[iB], R)
      end
   end
   return S
end

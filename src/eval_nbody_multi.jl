using JuLIP, NeighbourLists
using JuLIP: AbstractCalculator
using JuLIP.Potentials: @pot
using StaticArrays
using BenchmarkTools

using NBodyIPs: NBodyFunction,
                _get_loop_ex,
                _get_Jvec_ex,
                evaluate_I,
                NBSiteDescriptor,
                descriptor,
                ricoords,
                skip_simplex,
                fcut,
                fcut_d,
                invariants,
                invariants_ed,
                evaluate_I_ed,
                gradri2gradR!

using NeighbourLists: nbodies,
                    maptosites!,
                    maptosites_d!,
                    virial!,
                    max_neigs,
                    sites

import JuLIP: site_energies,
              energy,
              forces,
              virial

import NBodyIPs: evaluate_d!,
                 evaluate_many!,
                 evaluate_many_d!
                 # eval_site_nbody!,
                 # evaluate,



# skip the simplex if not the right species
function skip_simplex_species(Spi::Int,Spj::Vector{Int},Species::Vector{Int},J)
   error("skip simplex")
   Sp = [Spi]
   for i=1:length(J)
      push!(Sp,Spj[J[i]])
   end
   return sort(Sp) != sort(Species)
end

# skip the simplex if not the right species
function skip_simplex_species_many(Spi::Int,Spj::Vector{Int},Species::Vector{Int},J)
   error("skip simplex2")
   Sp = [Spi]
   for i=1:length(J)
      push!(Sp,Spj[J[i]])
   end
   return [sort(Sp) == sort(Species[k]) for k=1:length(Species)]
end


@generated function eval_site_nbodyM!(::Val{N},
                          Rs::AbstractVector{JVec{T}},
                          rcut::T,
                          reducefun,
                          out,
                          temp,
                          Spi::Int,
                          Spj::Vector{Int},
                          Species::Vector{Int}) where {N, T}

   code = Expr[]
   # initialise the output
   push!(code, :( print(".") ))
   push!(code, :( nR = length(Rs)  ))
   push!(code, :( println("nR = ", nR) ))

   # generate the multi-for-loop
   ex_loop = _get_loop_ex(N)

   # inside the loop
   # ---------------
   code_inner = Expr[]
   # collect the indices into a vector
   push!(code_inner,      _get_Jvec_ex(N) )

   # now call `V` with the simplex-corner vectors and "add" this to the site energy
   push!(code_inner, :( print(out) ))
   push!(code_inner, :(   out = reducefun(out, Rs, J, temp,Spi,Spj,Species) ))
   push!(code_inner, :( println(" -> ", out) ))

   # put code_inner into the loop expression
   ex_loop.args[2] = Expr(:block, code_inner...)

   # now append the loop to the main code
   push!(code, ex_loop)

   quote
      @inbounds $(Expr(:block, code...))
      return out
   end
end


function evaluateM(V::NBodyFunctionM{N},
                  desc::NBSiteDescriptor,
                  Rs::AbstractVector{JVec{T}},
                  J::SVector{K, Int},
                  Spi,Spj,Species) where {N, T, K}
   # check species
   error("evaluateM")
   skip_simplex_species(Spi,Spj,Species,J) && return zero(T)
   evaluate(V,desc,Rs,J)
end


evaluateM(V::NBodyFunctionM,Rs::AbstractVector{JVec{T}},J::SVector{K, Int},Spi::Int,Spj::Vector{Int},Species::Vector{Int}) where {T,K} = evaluateM(V,descriptor(V),Rs,J,Spi,Spj,Species)

function evaluate_d!(dVsite,
                     V::NBodyFunctionM{N},
                     desc::NBSiteDescriptor,
                     Rs::AbstractVector{JVec{T}},
                     J,
                     Spi::Int,Spj::Vector{Int},Species::Vector{Int}) where {N,T}
   # check species
   skip_simplex_species(Spi,Spj,Species,J) && return dVsite
   evaluate_d!(dVsite, V, desc, Rs, J)
end








function site_energies(V::NBodyFunctionM{N}, at::Atoms{T},Species::Vector{Int}) where {N, T}
   println("site_energy_multi")
   Es = zeros(T, length(at))
   Z = atomic_numbers(at)
   for (i, j, r, R) in sites(at, cutoff(V))
      Spi = Z[i]
      Spj = Z[j]
      # println("here I am")
      Es[i] = eval_site_nbodyM!(
                 Val(N), R, cutoff(V),
                (out, R, J, temp,Spi,Spj,Species) -> (out +
                                 evaluateM(V, descriptor(V), R,
                                           J,Spi,Spj,Species)),
                zero(T), nothing, Spi,Spj,Species)
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
      eval_site_nbodyM!(Val(N), R, cutoff(V),
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
      eval_site_nbodyM!(Val(N), R, rcut,
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
      eval_site_nbodyM!(Val(N), R, rcut,
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
      eval_site_nbodyM!(Val(N), R, rcut,
                       (out, R, J, temp, Spi, Spj, Species) -> evaluate_many_d!(out, B, R, J, Spi, Spj, Species),
                       dVsite, nothing,Spi,Spj,Species)
      # update the virials
      for iB = 1:nB
         S[iB] += JuLIP.Potentials.site_virial(dVsite[iB], R)
      end
   end
   return S
end

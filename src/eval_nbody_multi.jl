using JuLIP, NeighbourLists
using JuLIP: AbstractCalculator
using JuLIP.Potentials: @pot
using StaticArrays
using BenchmarkTools

using NBodyIPs: NBodyFunction, bapolys, eval_site_nbody!, evaluate, eval_site_nbody!, evaluate_d!, NBSiteDescriptor, _get_loop_ex, _get_Jvec_ex, descriptor, ricoords, skip_simplex, fcut, invariants, evaluate_I, fcut_d, invariants_ed, evaluate_I_ed, gradri2gradR!

import JuLIP: site_energies, energy, forces

import NBodyIPs: evaluate, eval_site_nbody!, evaluate_d!

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
   @show
   return sort(Sp) != sort(Species)
end


@generated function eval_site_nbody!( ::Val{N},
                                      Rs::AbstractVector{JVec{T}},
                                      rcut::T,
                                      reducefun,
                                      out,
                                      temp,
                                      Spi,Spj,Species ) where {N, T}
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


function site_energies(V::NBodyFunction{N}, at::Atoms{T},Species::Vector{Int}) where {N, T}
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

energy(V::NBodyFunction, at::Atoms, Species::Vector{Int}) = sum_kbn(site_energies(V, at, Species))

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


function forces(V::NBodyFunction{N}, at::Atoms{T},Species::Vector{Int}) where {N, T}
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

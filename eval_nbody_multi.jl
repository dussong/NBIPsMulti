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

# TODO: build potential with several species

function skip_simplex_species(Spi,Spj,Species,J)
   Sp = [Spi]
   for i=1:length(J)
      push!(Sp,Spj[J[i]])
   end
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
   # get the physical descriptor: bond-lengths (+ bond-angles)
   rθ = ricoords(desc, Rs, J)
   # check whether to skip this N-body term?
   skip_simplex(desc, rθ) && return zero(T)
   # compute the cutoff (and skip this site if the cutoff is zero)
   fc = fcut(desc, rθ)
   fc == 0 && return zero(T)
   # compute the invariants (this also applies the transform)
   II = invariants(desc, rθ)
   # evaluate the inner potential function (e.g. polynomial)
   return evaluate_I(V, II) * fc
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


function evaluate_d!(dVsite,
                     V::NBodyFunction{N},
                     desc::NBSiteDescriptor,
                     Rs,
                     J,
                     Spi::Int,Spj::Vector{Int},Species::Vector{Int}) where {N}
   # check species
   skip_simplex_species(Spi,Spj,Species,J) && return zero(T)
   evaluate_d!(dVsite, V, Rs, J)
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


r0 = 2.5
V = bapolys(2, "($r0/r)^4", "(:cos, 3.6, 4.8)", 2)
Vcucu = V[2]
at = bulk(:Cu, cubic=true)*2
rattle!(at,0.1)

norm(site_energies(Vcucu, at, [29,29])-site_energies(Vcucu, at))

@btime site_energies(Vcucu, at, [29,29])
@btime site_energies(Vcucu, at)

norm(forces(Vcucu, at, [29,29])-forces(Vcucu, at))
@btime forces(Vcucu,at,[29,29])
@btime forces(Vcucu,at)





# # implement site_energies for given species
# # Pair potential
# function site_energies2(V::NBodyFunction{2}, at::Atoms{T},
#                                        species::Tuple{Int,Int}) where {T}
#    Z = atomic_numbers(at)
#    Es = zeros(T, length(at))
#    for (i, j, r, R) in sites(at, cutoff(V))
#       for k = 1:length(j)
#          atnb = sort([Z[i],Z[j][k]])
#          if (atnb[1] == species[1])&(atnb[2] == species[2])
#             Es[i] += V(r[k])
#          end
#       end
#    end
#    return Es
# end

# using symbols
# function site_energies(V::NBodyFunction{2}, at::Atoms{T},
#                                        species::Tuple{Symbol,Symbol}) where {T}
#    sp = atomic_number.(species)
#    return site_energies(V, at, sp)
# end

# energy(V::NBodyFunction, at::Atoms,species::Tuple{Symbol,Symbol}) = sum_kbn(site_energies(V, at,species))
#
# # Implementation of the forces
# function forces(V::NBodyFunction{2}, at::Atoms{T},sp::Tuple{Int,Int}) where { T}
#    nlist = neighbourlist(at, cutoff(V))
#    maxneigs = max_neigs(nlist)
#    F = zeros(JVec{T}, length(at))
#    dVsite = zeros(JVec{T}, maxneigs)
#    for (i, j, r, R) in sites(nlist)
#       fill!(dVsite, zero(JVec{T}))
#       eval_site_nbody!(
#             Val(2), R, cutoff(V),
#             (out, R, J, temp) -> evaluate_d!(out, V, R, J),
#             dVsite, nothing )   # dVsite == out, nothing == temp
#       # write site energy gradient into forces
#       for n = 1:length(j)
#          F[j[n]] -= dVsite[n]
#          F[i] += dVsite[n]
#       end
#    end
#    return F
# end




# r0 = 2.5
# V = bapolys(2, "($r0/r)^4", "(:cos, 3.6, 4.8)", 2)
# Vcucu = V[2]
# at = bulk(:Cu, cubic=true)
# species = (:Cu,:Zn)
#
# atomic_number(:Cu)
#
# Vcucu(3.)
#
# site_energies2(Vcucu, at, (29,29))
# site_energies(Vcucu, at, (:Cu,:Cu))
# energy(Vcucu,at,(:Cu,:Cu))
# forces(Vcucu,at,(29,29))

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

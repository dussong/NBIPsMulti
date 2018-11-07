
using NBodyIPs: NBodyFunction
using StaticArrays

abstract type NBodyFunctionM{N, DT} <: NBodyFunction{N,DT} end

import Base:              length,
                          Dict,
                          ==
import JuLIP:             cutoff
import JuLIP.Potentials:  @pot
import NBodyIPs:          fast,
                          degree,
                          combinebasis,
                          descriptor,
                          combiscriptor,
                          evaluate_many!,
                          evaluate_many_d!,
                          evaluate_I,
                          evaluate_I_ed,
                          basisname


export NBPolyM

const VecTup{M} = Vector{NTuple{M, Int}}
const Tup{M} = NTuple{M, Int}
# ==================================================================
#           Polynomials of Invariants
# ==================================================================

@pot struct NBPolyM{N, M, T, TD} <: NBodyFunction{N, TD}
   t::VecTup{M}               # tuples M = #edges + 1
   c::Vector{T}               # coefficients
   D::TD                      # Descriptor
   valN::Val{N}               # encodes that this is an N-body function
   Sp::Vector{Int}            #encodes the species

   NBPoly(t::VecTup{M}, c::Vector{T}, D::TD, valN::Val{N}, Sp::Vector{Int}) where {N, M, T, TD} = (
      N <= 1 ? error("""NBPoly must have body-order 2 or larger;
                        use `NBodyIPs.OneBody{T}` for 1-body.""")
             : new{N, M, T, TD}(t, c, D, valN,Sp))
end

"""
`struct NBPolyM`  (N-Body Polynomial with attached species, slow implementation)

A `struct` storing the information for a (pure) N-body potential, i.e.,
containing *only* terms of a specific body-order. Several `NBPoly`s can be
combined into an interatomic potential via `NBodyIP`.

### Fields

* `t::Vector{NTuple{M,TI}}` : list of M-tuples containing basis function information
e.g., if M = 7, α = t[n] is a 7-vector then this corresponds to the basis function
```
I2[α[7]] * ∏_{j=1..6} I1[j]^α[j]
```
where `I1, I2` are the 4-body invariants.

* `c::Vector{T}`: vector of coefficients for the basis functions

* `D`: a descriptor (cf `NBodyIPs.NBodyDescriptor`)

* `Sp`: a vector of Int containing the species
"""
NBPoly

==(V1::NBPolyM, V2::NBPolyM) = ( (V1.t == V2.t) && (V1.c == V2.c) && (V1.D == V2.D) && (V1.Sp == V2.Sp) )

descriptor(V::NBPolyM) = V.D

species(V::NBPolyM) = V.Sp

basisname(::NBPolyM) = "NBPolyM"

combiscriptor(V::NBPolyM) = (NBPolyM, bodyorder(V), combiscriptor(V.D), Sp)

# standard constructor (N can be inferred)
NBPolyM(t::VecTup{K}, c, D,Sp) where {K} = NBPolyM(t, c, D, Val(edges2bo(K-1)),Sp)

# NBPolyM made from a single basis function rather than a collection
NBPolyM(t::Tup, c, D,) = NBPolyM([t], [c], D, Sp)

# collect multiple basis functions represented as NBPoly's into a single NBPolyM
# (for performance reasons)
# TODO: this is not general enough!
NBPolyM(B::Vector{TB}, c, D, Sp) where {TB <: NBPolyM} =
      NBPolyM([b.t[1] for b in B], c .* [b.c[1] for b in B], D, Sp)

# 1-body term (on-site energy)
NBPolyM(c::Float64) = NBPolyM([Tup{0}()], [c], nothing, Val(1), Sp)

# number of basis functions which this term is made from
length(V::NBPolyM) = length(V.t)

cutoff(V::NBPolyM) = cutoff(V.D)

function match_dictionary(V::NBPolyM, V1::NBPolyM)
   if V.D != V1.D
      if V.D.s != V1.D.s
         if V.Sp != V1.Sp
            warn("matching two non-matching dictionaries!")
         end
      end
   end
   return NBPolyM(V.t, V.c, V1.D, V.valN, V.Sp)
end

combinebasis(basis::AbstractVector{TV}, coeffs) where {TV <: NBPolyM} =
      NBPolyM(basis, coeffs, basis[1].D, Sp)


function degree(V::NBPolyM)
   if length(V) == 1
      return NBodyIPs.PolyBasis.tdegree(descriptor(V), V.t[1])
   end
   error("`degree` is only defined for `NBPoly` basis functions, length == 1")
end


function Base.info(B::Vector{T}; indent = 2) where T <: NBPolyM
   ind = repeat(" ", indent)
   println(ind * "body-order = $(bodyorder(B[1]))")
   println(ind * "    length = $(length(B))")
   println(ind * "    Species = $(species(B))")
   if bodyorder(B[1]) > 1
      println(ind * " transform = $(B[1].D.s[1])")
      println(ind * "    cutoff = $(B[1].D.s[2])")
   end
end



# ---------------------------------------------------------------------
#  Construction of basis functions
# ---------------------------------------------------------------------


nbpolys(N::Integer, desc, tdeg; kwargs...) =
   nbpolys(gen_tuples(desc, N, tdeg; kwargs...), desc)

nbpolys(ts::VecTup, desc) = [NBPoly(t, 1.0, desc) for t in ts]

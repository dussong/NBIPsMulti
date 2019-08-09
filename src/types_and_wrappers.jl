

using NBodyIPs: NBodyFunction
using StaticArrays

using NBodyIPs:        BondLengthDesc,
                       BondAngleDesc,
                       edges2bo,
                       bo2edges,
                       _decode_dict,
                       bodyorder,
                       SpaceTransform,
                       NBCutoff,
                       NBSiteDescriptor

using NBodyIPs.Polys:  Tup,
                       VecTup,
                       NBPoly,
                       monomial,
                       monomial_d

using NBodyIPs.PolyBasis: gen_tuples

import Base:              length,
                          Dict,
                          ==
import JuLIP:             cutoff
import JuLIP.Potentials:  @pot
import NBodyIPs:          fast,
                          degree,
                          combinebasis,
                          descriptor,
                          # combiscriptor,
                          evaluate_many!,
                          evaluate_many_d!,
                          evaluate_I,
                          evaluate_I_ed,
                          basisname,
                          nbpolys,
                          NBCutoff,
                          bodyorder,
                          evaluate,
                          evaluate_d
import NBodyIPs.Regularisers: species
import NBodyIPs.Polys: info

import NBodyIPs.PolyBasis: gen_tuples,
                           tdegree

import Base: hash
export NBPolyM, NBodyFunctionM, MultiDesc


const BASIS = Val{:basis}
# tdegrees(desc::MultiDesc, vN) = NBIPsMulti.MultiInvariants.tdegrees(vN)

# ==================================================================
#          Type for species NBodyFunction
# ==================================================================

abstract type NBodyFunctionM{N, DT, SP} <: NBodyFunction{N,DT} end


struct MultiDesc{TT <: SpaceTransform,
                 TC <: NBCutoff, SP, N} <: NBSiteDescriptor
   transform::TT
   cutoff::TC
   sp_type::Val{SP}
   valN::Val{N} #encodes the body-order
end

# ==================================================================
#           Polynomials of Invariants
# ==================================================================

struct NBPolyM{N, M, T, TD, SP} <: NBodyFunctionM{N, TD, SP}
   t::VecTup{M}               # tuples M = #edges + 1
   c::Vector{T}               # coefficients
   D::TD                      # Descriptor
   valN::Val{N}               # encodes that this is an N-body function
   Sp::Vector{Int}            #encodes the species
   Sp_type::Val{SP}

   NBPolyM(t::VecTup{M}, c::Vector{T}, D::TD, valN::Val{N}, Sp::Vector{Int}, Sp_type::Val{SP}) where {N, M, T, TD, SP} = (
      N <= 1 ? error("""NBPoly must have body-order 2 or larger;
                        use `NBodyIPs.OneBody{T}` for 1-body.""")
             : new{N, M, T, TD, SP}(t, c, D, valN, Sp, Sp_type))
end

@pot NBPolyM


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

* `Sp_type`
"""

==(V1::NBPolyM, V2::NBPolyM) = ( (V1.t == V2.t) && (V1.c == V2.c) && (V1.D == V2.D) && (V1.Sp == V2.Sp)  && (V1.Sp_type == V2.Sp_type) )

descriptor(V::NBPolyM) = V.D

species(V::NBPolyM) = V.Sp

species_type(V::NBPolyM) = V.Sp_type

basisname(::NBPolyM) = "NBPolyM"

bodyorder(V::NBPolyM{N, M, T, TD, SP}) where {N, M, T, TD, SP} = N


# combiscriptor(V::NBPolyM) = (NBPolyM, bodyorder(V), combiscriptor(V.D), V.Sp, V.Sp_type)

hash(::BASIS, V::NBPolyM) =
   hash( (hash(NBPolyM), hash(bodyorder(V)), hash(BASIS(), V.D),
        hash(V.Sp), hash(V.Sp_type)))

# standard constructor (N can be inferred)
NBPolyM(t::VecTup{K}, c, D, Sp, Sp_type) where {K} = NBPolyM(t, c, D, Val(edges2bo(K-1)), Sp, Sp_type)

# NBPolyM made from a single basis function rather than a collection
NBPolyM(t::Tup, c, D, Sp, Sp_type) = NBPolyM([t], [c], D, Sp, Sp_type)

# collect multiple basis functions represented as NBPoly's into a single NBPolyM
# (for performance reasons)
# TODO: this is not general enough!
NBPolyM(B::Vector{TB}, c, D, Sp, Sp_type) where {TB <: NBPolyM} =
      NBPolyM([b.t[1] for b in B], c .* [b.c[1] for b in B], D, Sp, Sp_type)

# 1-body term (on-site energy)
NBPolyM(c::Float64) = NBPolyM([Tup{0}()], [c], nothing, Val(1), [], Val{:OneB})

# number of basis functions which this term is made from
length(V::NBPolyM) = length(V.t)

cutoff(V::NBPolyM) = cutoff(V.D)

# function match_dictionary(V::NBPolyM, V1::NBPolyM)
#    if V.D != V1.D
#       if V.D.s != V1.D.s
#          if V.Sp != V1.Sp
#             if V.Sp_type != V1.Sp_type
#                @warn("matching two non-matching dictionaries!")
#             end
#          end
#       end
#    end
#    return NBPolyM(V.t, V.c, V1.D, V.valN, V.Sp, V.Sp_type)
# end



combinebasis(basis::AbstractVector{TV}, coeffs) where {TV <: NBPolyM} =
      NBPolyM(basis, coeffs, basis[1].D, basis[1].Sp, basis[1].Sp_type)


function degree(V::NBPolyM)
   if length(V) == 1
      return NBodyIPs.PolyBasis.tdegree(descriptor(V), V.t[1])
   end
   error("`degree` is only defined for `NBPoly` basis functions, length == 1")
end


function info(B::Vector{T}; indent = 2) where T <: NBPolyM
   ind = repeat(" ", indent)
   println(ind * "body-order = $(bodyorder(B[1]))")
   println(ind * "    length = $(length(B))")
   println(ind * "    Species = $(species(B))")
   println(ind * "    Species type = $(species_type(B))")
   if bodyorder(B[1]) > 1
      println(ind * " transform = $(B[1].D.s[1])")
      println(ind * "    cutoff = $(B[1].D.s[2])")
   end
end


# ---------------  evaluate the n-body terms ------------------

# a tuple α = (α1, …, α6, α7) means the following:
# with f[0] = 1, f[1] = I7, …, f[5] = I11 we construct the basis set
#   f[α7] * g(I1, …, I6)
# this means, that gen_tuples must generate 7-tuples instead of 6-tuples
# with the first 6 entries restricted by degree and the 7th tuple must
# be in the range 0, …, 5

function evaluate_I(V::NBPolyM, II)
   I1, I2 = II
   E = zero(eltype(I1))
   for (α, c) in zip(V.t, V.c)
      E += c * I2[1+α[end]] * monomial(α, I1)
   end
   return E
end

function evaluate_I_ed(V::NBPolyM, II)
   I1, I2, dI1, dI2 = II
   E = zero(eltype(I1))
   dM = zero(typeof(I1))
   dE = zero(typeof(I1))
   #
   for (α, c) in zip(V.t, V.c)
      m, m_d = monomial_d(α, I1)
      E += c * I2[1+α[end]] * m        # just the value of the function itself
      dM += (c * I2[1+α[end]]) * m_d   # the I2 * ∇m term without the chain rule
      dE += (c * m) * dI2[1+α[end]]    # the ∇I2 * m term
   end
   # chain rule
   for i = 1:length(dI1)   # dI1' * dM
      dE += dM[i] * dI1[i]
   end
   return E, dE
end

# # -------------- Body-order from species --------
# bodyorder(::Val{:AA}) = 2
# bodyorder(::Val{:AAA}) = 3
# bodyorder(::Val{:AAB}) = 3
# bodyorder(::Val{:ABC}) = 3

# -------------- Infrastructure to read/write NBPolyM  --------
Dict(V::NBPolyM{N, M, T, TD, SP}) where {N, M, T, TD, SP} = Dict(
                        "__id__" => "NBPolyM",
                        "t" => V.t,
                        "c" => V.c,
                        "D" => Dict(V.D),
                        "N" => N,
                        "Sp" => V.Sp,
                        "Sp_type" => String(SP) )

NBPolyM(D::Dict) = NBPolyM([ tuple(ti...) for ti in D["t"] ],
                           Vector{Float64}(D["c"]),
                           _decode_dict(D["D"]),
                           Val(D["N"]),
                           Vector{Int}(D["Sp"]),
                           Val(Symbol(D["Sp_type"])))

Base.convert(::Val{:NBPolyM}, D::Dict) = NBPolyM(D)

#
# # ==================================================================
# #    StNBPoly
# # ==================================================================
#
# @pot struct StNBPoly{N, TD, TP} <: NBodyFunction{N, TD}
#    D::TD       # Descriptor
#    P::TP       # a static polynomial
#    valN::Val{N}
# end
#
# """
# `struct StNBPoly`  (N-Body Bond-length Polynomial)
#
# fast evaluation of the outer polynomial using `StaticPolynomials`
# """
# StNBPoly
#
# descriptor(V::StNBPoly) = V.D
#
# function StNBPoly(V::NBPoly{N}) where {N}
#    nI1, nI2 = ninvariants(V.D, N)
#    nI = nI1 + nI2
#    nmonomials = length(V.c)
#    # generate the exponents for the StaticPolynomial
#    exps = zeros(Int, nI, nmonomials)
#    for (i, α) in enumerate(V.t)  # i = 1:nmonomials
#       for (j, a) in enumerate(α[1:end-1])   #  ∏ I1[j]^α[j]
#          exps[j, i] = a    # I1[j]^α[j]
#       end
#       exps[nI1+1+α[end], i] = 1   # I2[α[end]] * (...)
#    end
#    # generate the static polynomial
#    return StNBPoly(V.D, StaticPolynomials.Polynomial(V.c, exps), V.valN)
# end
#
# cutoff(V::StNBPoly) = cutoff(V.D)
#
# fast(Vn::StNBPoly)  = Vn
# fast(Vn::NBPoly) =  StNBPoly(Vn)
#
# evaluate_I(V::StNBPoly, II) =
#       StaticPolynomials.evaluate(V.P, vcat(II...))
#
# function evaluate_I_ed(V::StNBPoly, II)
#    V, dV_dI = StaticPolynomials.evaluate_and_gradient(V.P, vcat(II[1], II[2]))
#    if length(dV_dI) != length(II[3]) + length(II[4])
#       @show size.(II)
#       @show size(dV_dI)
#       @show size(vcat(II[1], II[2]))
#    end
#    return V, dot(vcat(II[3], II[4]), dV_dI)  # (dI' * dV_dI)
# end
#



# ---------------------------------------------------------------------
#  Construction of basis functions
# ---------------------------------------------------------------------
function tdegree(desc::MultiDesc, α)
   K = length(α)
   @assert Val(edges2bo(K-1)) == desc.valN
   degs1, degs2 = tdegrees(desc.sp_type)
   # primary invariants
   d = sum(α[j] * degs1[j] for j = 1:K-1)
   # secondary invariants
   d += degs2[1+α[end]]
   return d
end


function gen_tuples(desc::MultiDesc, vN::Val{N}, vK::Val{K}, deg, tuplebound) where {N, K}
   A = Tup{K}[]
   degs1, degs2 = tdegrees(desc, vN)

   α = @MVector zeros(Int, K)
   α[1] = 1
   lastinc = 1

   while true
      admit_tuple = false
      if α[end] <= length(degs2)-1
         if tuplebound(α)
            admit_tuple = true
         end
      end
      if admit_tuple
         push!(A, SVector(α).data)
         α[1] += 1
         lastinc = 1
      else
         if lastinc == K
            return A
         end
         α[1:lastinc] .= 0
         α[lastinc+1] += 1
         lastinc += 1
      end
   end
   error("I shouldn't be here!")
end



function gen_tuples(desc::MultiDesc, deg; tuplebound = (α -> (0 < tdegree(desc, α) <= deg)))
   N = desc.valN
   vN = desc.sp_type
   return gen_tuples(desc, vN, Val(bo2edges(N)+1), deg, tuplebound)
end


function nbpolys(desc::MultiDesc, tdeg, Sp, Sp_type)
   @assert Sp_type == desc.sp_type
   @assert Val(length(Sp)) == desc.valN
   @assert check_species(desc, Sp)
   nbpolys(gen_tuples(desc, tdeg),desc, Sp, Sp_type)
end

nbpolys(ts::VecTup, desc, Sp, Sp_type) = [NBPolyM(t, 1.0, desc, Sp, Sp_type) for t in ts]

nbpolys(desc::MultiDesc, tdeg, Sp) = nbpolys(desc, tdeg, Sp, desc.sp_type)

function check_species(desc::MultiDesc, Sp, ::Val{:AA})
   return length(Sp) == 2
end

function check_species(desc::MultiDesc, Sp, ::Val{:AAAba})
   return (length(Sp) == 3)&&(Sp[1] == Sp[2])&&(Sp[2] == Sp[3])
end

function check_species(desc::MultiDesc, Sp, ::Val{:AAA})
   return (length(Sp) == 3)&&(Sp[1] == Sp[2])&&(Sp[2] == Sp[3])
end

function check_species(desc::MultiDesc, Sp, ::Val{:AAB})
   if length(Sp) == 3
      if (Sp[1] == Sp[2])&&(Sp[2] !== Sp[3])
         return true
      elseif (Sp[2] == Sp[3])&&(Sp[1] !== Sp[3])
         return true
      elseif (Sp[1] == Sp[3])&&(Sp[2] !== Sp[3])
         return true
      else
         return false
      end
   else
      return false
   end
end

function check_species(desc::MultiDesc, Sp, ::Val{:AABba})
   if length(Sp) == 3
      if (Sp[1] == Sp[2])&&(Sp[2] !== Sp[3])
         return true
      elseif (Sp[2] == Sp[3])&&(Sp[1] !== Sp[3])
         return true
      elseif (Sp[1] == Sp[3])&&(Sp[2] !== Sp[3])
         return true
      else
         return false
      end
   else
      return false
   end
end

function check_species(desc::MultiDesc, Sp, ::Val{:ABC})
   return (length(Sp) == 3)&&(length(unique(Sp)) == 3)
            # (Sp[1] != Sp[2])&&
            # (Sp[2] != Sp[3])&&
            # (Sp[1] != Sp[3]))
end

function check_species(desc::MultiDesc, Sp, ::Val{:AAAA})
   return (length(Sp) == 4)&&(length(unique(Sp)) == 1)
end

function check_species(desc::MultiDesc, Sp, ::Val{:AAAAba})
   return (length(Sp) == 4)&&(length(unique(Sp)) == 1)
end

function check_species(desc::MultiDesc, Sp, ::Val{:AAAB})
   if (length(Sp) == 4)&&(length(unique(Sp)) == 2)
      if (Sp[1] == Sp[2])&&(Sp[2] == Sp[3])&&(Sp[3] !== Sp[4])
         return true
      elseif (Sp[1] !== Sp[2])&&(Sp[2] == Sp[3])&&(Sp[3] == Sp[4])
         return true
      elseif (Sp[1] == Sp[3])&&(Sp[3] == Sp[4])&&(Sp[1] !== Sp[2])
         return true
      elseif (Sp[1] == Sp[2])&&(Sp[2] == Sp[4])&&(Sp[1] !== Sp[3])
         return true
      end
   else
      return false
   end
end

function check_species(desc::MultiDesc, Sp, ::Val{:AAABba})
   if (length(Sp) == 4)&&(length(unique(Sp)) == 2)
      if (Sp[1] == Sp[2])&&(Sp[2] == Sp[3])&&(Sp[3] !== Sp[4])
         return true
      elseif (Sp[1] !== Sp[2])&&(Sp[2] == Sp[3])&&(Sp[3] == Sp[4])
         return true
      elseif (Sp[1] == Sp[3])&&(Sp[3] == Sp[4])&&(Sp[1] !== Sp[2])
         return true
      elseif (Sp[1] == Sp[2])&&(Sp[2] == Sp[4])&&(Sp[1] !== Sp[3])
         return true
      end
   else
      return false
   end
end

function check_species(desc::MultiDesc, Sp, ::Val{:AABB})
   if (length(Sp) == 4)&&(length(unique(Sp)) == 2)
      if (Sp[1] == Sp[2])&&(Sp[3] == Sp[4])&&(Sp[2] !== Sp[3])
         return true
      elseif (Sp[1] == Sp[3])&&(Sp[2] == Sp[4])&&(Sp[1] !== Sp[2])
         return true
      elseif (Sp[1] == Sp[4])&&(Sp[2] == Sp[3])&&(Sp[1] !== Sp[2])
         return true
      end
   else
      return false
   end
end

function check_species(desc::MultiDesc, Sp, ::Val{:AABC})
   return (length(Sp) == 4)&&(length(unique(Sp)) == 3)
end

function check_species(desc::MultiDesc, Sp, ::Val{:ABCD})
   return (length(Sp) == 4)&&(length(unique(Sp))==4)
end


function check_species(desc::MultiDesc, Sp)
   return (check_species(desc::MultiDesc, Sp, desc.sp_type))&&(Sp == sort(Sp))
end


# Evaluate functions
evaluate(V::NBodyFunctionM{2}, r::T) where {T <: AbstractFloat} =
      evaluate(V, SVector(r))
evaluate_d(V::NBodyFunctionM{2}, r::T) where {T <: AbstractFloat} =
      evaluate_d(V, SVector(r))[1]

function evaluate(V::NBodyFunctionM{2}, r::SVector{1})
   D = descriptor(V)
   return evaluate_I(V, invariants(D, r)) * fcut(D, r)
end

function evaluate_d(V::NBodyFunction{2}, r::SVector{1})
   D = descriptor(V)
   fc, fc_d = fcut_d(D, r)
   Vn, Vn_d = evaluate_I_ed(V, invariants_ed(D, r))
   return fc * Vn_d + fc_d * Vn
end

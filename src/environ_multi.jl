# TODO:
#   - allow Vn to be an arbitrary pair potential
#   - replace the polynomial with an arbitrary family of pair potentials
#   - (or even nbody?)

module EnvIPsmulti

using NBIPsMulti:         MultiDesc, NBPolyM, species_type
using StaticArrays
using JuLIP:              AbstractCalculator
using JuLIP.Potentials:   Shift,
                          @analytic,
                          @pot
using NBodyIPs:           NBodyDescriptor
using NBodyIPs.PolyBasis: nbpolys
using NBodyIPs.EnvIPs:           analyse_Vn

import Base:              Dict,
                          ==,
                          convert
import NBodyIPs:          NBodyIP,
                          bodyorder,
                          fast,
                          combinebasis,
                          _decode_dict,
                          descriptor,
                          # combiscriptor,
                          degree,
                          basisname,
                          bodyorder
import NBodyIPs.Regularisers: species
import NBIPsMulti:        species_type
import NBodyIPs.Polys: info

import Base: hash

export envpolysM

const BASIS = Val{:basis}

abstract type AbstractEnvIP{N} <: AbstractCalculator end

struct EnvIPM{N, P, TVR, TVN, SP} <: AbstractEnvIP{N}
   t::Int
   Vr::TVR     # N-body potential -> multi
   Vn::TVN     # neighbour counter
   str_Vn::String  # string describing the neighbour counter
   weights::Dict{Tuple{Int64,Int64},Float64}
   valN::Val{N}
   valP::Val{P}
   Sp::Vector{Int}            #encodes the species
   sp_type::Val{SP}
end

@pot EnvIPM


==(V1::EnvIPM, V2::EnvIPM) = ( (V1.t == V2.t) &&
                             (V1.Vr == V2.Vr) &&
                             (V1.str_Vn == V2.str_Vn) &&
                             (V1.weights == V2.weights) &&
                             (V1.valN == V2.valN) &&
                             (V1.valP == V2.valP) &&
                             (V1.sp_type == V2.sp_type))

Dict(V::EnvIPM{N, P, TVR, TVN, SP}) where {N, P, TVR, TVN, SP}  =
                   Dict( "__id__" => "EnvIPM",
                       "t" => V.t,
                       "Vr" => Dict(V.Vr),
                       "str_Vn" => V.str_Vn,
                       "weights" => V.weights,
                       "cutoff_Vn" => cutoff(V.Vn),
                       "N" => N,
                       "P" => P,
                       "Sp" => V.Sp,
                       "Sp_type" => String(SP), )

species(V::EnvIPM) = V.Sp

species_type(V::EnvIPM) = V.sp_type

bodyorder(V::EnvIPM{N, P, TVR, TVN, SP}) where {N, P, TVR, TVN, SP} = N

function _decode_weights(D::Dict)
   Dout = Dict{Tuple{Int64,Int64},Float64}()
   for k in keys(D)
      if typeof(k) == Tuple{Int64,Int64}
         Dout[k] = D[k]
      else
         k1 = parse(Int64,split(split(split(k,"(")[2],")")[1],",")[1])
         k2 = parse(Int64,split(split(split(k,"(")[2],")")[1],",")[2])
         Dout[(k1,k2)] = D[k]
      end
   end
   return Dout
end


EnvIPM(D::Dict) = EnvIPM( D["t"],
                        _decode_dict(D["Vr"]),
                        analyse_Vn(D["str_Vn"], D["cutoff_Vn"]),
                        D["str_Vn"],
                        _decode_weights(D["weights"]),
                        Val(D["N"]),
                        Val(D["P"]),
                        Vector{Int}(D["Sp"]),
                        Val(Symbol(D["Sp_type"])))


convert(::Val{:EnvIPM}, D::Dict) = EnvIPM(D)

function EnvIPM(t, Vr,
                str_Vn::String,
                cutoff_Vn::AbstractFloat,
                weights::Dict{Tuple{Int64,Int64},Float64},
                Sp::Vector{Int},
                sp_type)
   Vn = analyse_Vn(str_Vn, cutoff_Vn)
   return EnvIPM(t, Vr, Vn, str_Vn, weights, Sp, sp_type)
end

EnvIPM(t::Int, Vr::NBPolyM, Vn, str_Vn::String, weights) =
      EnvIPM(t, Vr, Vn, str_Vn, weights, Val(bodyorder(Vr)), Val(t), species(Vr), species_type(Vr))

Vn(V::EnvIPM) = V.Vn
Vr(V::EnvIPM) = V.Vr

descriptor(V::EnvIPM) = descriptor(V.Vr)

# combiscriptor(V::EnvIPM) = (EnvIPM,
#                            combiscriptor(V.Vr),
#                            V.str_Vn,
#                            V.t,
#                            V.Sp, V.sp_type)

# TODO: add hash for weights
hash(::BASIS, V::EnvIPM) = hash((hash(EnvIPM),
                                 hash(V.valN),
                                 hash(V.valP),
                                 hash(BASIS(), V.Vr),
                                 hash(V.str_Vn),
                                 hash(V.t),
                                 hash(V.Sp),
                                 hash(V.sp_type)))


function degree(V::EnvIPM)
   if length(V.Vr) == 1
      return ( degree(V.Vr), V.t )
   end
   error("`degree` is only defined for `EnvIPM` basis functions, length == 1")
end

basisname(::EnvIPM) = "EnvIPM"

# ----------------- generate basis / IP / convert ----------------

function envpolysM(D::MultiDesc, deg_poly::Integer,
                  Vn_descr, deg_n::Integer, weights, Sp; kwargs...)
   B_poly = nbpolys(D, deg_poly, Sp; kwargs...)
   B = EnvIPM[]
   str_Vn = Vn_descr[1]
   Vn = analyse_Vn(Vn_descr...)
   for deg = 0:deg_n
      append!(B, [EnvIPM(deg, Vr, Vn, str_Vn, weights) for Vr in B_poly])
   end
   return [b for b in B]
end

function combinebasis(basis::AbstractVector{TV}, coeffs) where {TV <: EnvIPM}
   # @assert isleaftype(TV)
   @assert all( b.t == basis[1].t for b in basis )
   # combine the Vr components of the basis functions
   # (we get to do this because all t (=P) are the same
   Vr = combinebasis( [b.Vr for b in basis], coeffs )
   return EnvIPM(basis[1].t, Vr, basis[1].Vn, basis[1].str_Vn, basis[1].weights)
end


function info(B::Vector{T}; indent = 2) where T <: EnvIPM
   ind = repeat(" ", indent)
   println(ind * "EnvIPM with P = $(B[1].t)")
   println(ind * "           Vn : $(B[1].str_Vn)")
   println(ind * "           Vr : ...")
   info([ b.Vr for b in B ], indent = indent+5)
end

# fast(V::EnvIPM) = EnvIP(V.t, fast(V.Vr), V.Vn, V.str_Vn)



# # ========================================================
# #    Re-Implement Fast Version of EnvIPs
# # ========================================================
#
# @pot struct EnvPoly{N, TVR, TVN} <: AbstractEnvIP{N}
#    Vr::Vector{TVR}    # N-body potentials multiplied by n^j, j = 0, 1, ...
#    Vn::TVN            # neighbour counter
#    str_Vn::String     # string describing the neighbour counter
#    valN::Val{N}
# end
#
# envdegree(V::EnvPoly) = length(V.P) - 1
#
# ==(V1::EnvPoly, V2::EnvPoly) = ( (V1.Vr == V2.Vr) &&
#                                  (V1.str_Vn == V2.str_Vn) &&
#                                  (V1.valN == V2.valN) )
#
# # Dict(V::EnvPoly) = Dict( "__id__" => "EnvPoly",
# #                          "Vr" => Dict.(V.Vr),
# #                          "str_Vn" => V.str_Vn,
# #                          "cutoff_Vn" => cutoff(V.Vn)  )
# #
# # EnvPoly(D::Dict) = EnvPoly( _decode_dict.(D["Vr"]),
# #                             analyse_Vn(D["str_Vn"], D["cutoff_Vn"]),
# #                             D["str_Vn"] )
#
# convert(::Val{:EnvPoly}, D::Dict) = EnvPoly(D)
#
# EnvPoly(t, Vr, str_Vn::String, cutoff_Vn::AbstractFloat) =
#    EnvPoly(t, Vr, analyse_Vn(str_Vn, cutoff_Vn), str_Vn)
#
# descriptor(V::EnvPoly) = descriptor(V.Vr[1])
#
# combiscriptor(V::EnvPoly) = (EnvPoly,
#                              combiscriptor.(V.Vr),
#                              V.str_Vn)
#
# Vn(V::EnvPoly) = V.Vn
# Vr(V::EnvPoly) = V.Vr
#
# function EnvPoly(Vr::AbstractVector, Vn, str_Vn::String)
#    # for now we put the burden on the user to create a Vr vector where
#    # the types match But we can still collect them here for convenience
#    Vr_coll = [ v for v in Vr ]
#    if !isleaftype(typeof(Vr_coll))
#       error("""For now, `EnvPoly` can only be constructed from an array of
#                potentials that all share the same concrete type.""")
#    end
#    return EnvPoly(Vr_coll, Vn, str_Vn, Val(bodyorder(Vr[1])))
# end
#
#
# # function degree(V::EnvIP)
# #    if length(V.Vr) == 1
# #       return ( degree(V.Vr), V.t )
# #    end
# #    error("`degree` is only defined for `EnvIP` basis functions, length == 1")
# # end
#
# basisname(::EnvPoly) = "EnvPoly"




# all evaluation and assembly in here:
include("eval_env_multi.jl")


end

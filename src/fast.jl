# Extending fast functionalities

using NBodyIPs
using StaticPolynomials

import NBodyIPs: fast, ninvariants, descriptor, evaluate_I, evaluate_I_ed

export StNBPolyM

# ==================================================================
#    StNBPolyM
# ==================================================================


# """
# `struct StNBPolyM`  (N-Body Bond-length Polynomial)
#
# fast evaluation of the outer polynomial using `StaticPolynomials`
# """
struct StNBPolyM{N, TD, SP, TP} <: NBodyFunctionM{N, TD, SP}
   D::TD       # Descriptor
   valN::Val{N}
   Sp::Vector{Int}
   Sp_type::Val{SP}
   P::TP       # a static polynomial
end

@pot StNBPolyM

descriptor(V::StNBPolyM) = V.D

ninvariants(D::MultiDesc) = length.(tdegrees(D.sp_type))

function StNBPolyM(V::NBPolyM{N}) where {N}
   nI1, nI2 = ninvariants(V.D)
   nI = nI1 + nI2
   nmonomials = length(V.c)
   # generate the exponents for the StaticPolynomial
   exps = zeros(Int, nI, nmonomials)
   for (i, α) in enumerate(V.t)  # i = 1:nmonomials
      for (j, a) in enumerate(α[1:end-1])   #  ∏ I1[j]^α[j]
         exps[j, i] = a    # I1[j]^α[j]
      end
      exps[nI1+1+α[end], i] = 1   # I2[α[end]] * (...)
   end
   # generate the static polynomial
   return StNBPolyM(V.D,
                    V.valN, V.Sp, V.Sp_type,
                    StaticPolynomials.Polynomial(V.c, exps))
end

cutoff(V::StNBPolyM) = cutoff(V.D)

fast(Vn::StNBPolyM)  = Vn
fast(Vn::NBPolyM) =  StNBPolyM(Vn)

evaluate_I(V::StNBPolyM, II) =
      StaticPolynomials.evaluate(V.P, vcat(II...))

function evaluate_I_ed(V::StNBPolyM, II)
   V1, dV_dI = StaticPolynomials.evaluate_and_gradient(V.P, vcat(II[1], II[2]))
   # if length(dV_dI) != length(II[3]) + length(II[4])
   #    @show size.(II)
   #    @show size(dV_dI)
   #    @show size(vcat(II[1], II[2]))
   #    @show size(vcat(II[3], II[4]))
   # end

   # TODO: check a few variants how to compute dV
   #       and test performance
   II34 = vcat(II[3], II[4])
   dV = II34[1]*dV_dI[1]
   for i = 2:length(II34)
      dV += II34[i] * dV_dI[i]
   end

   return V1, dV # sum(i * dV_di for (i, dVdi) in zip( .* dV_dI)  # (dI' * dV_dI)
end

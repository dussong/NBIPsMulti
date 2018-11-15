using JuLIP: JVec
using NBodyIPs: SpaceTransform,
                Cutoff,
                edge_lengths,
                fcut,
                fcut_d,
                invariants,
                transform,
                invariants_ed,
                transform_d,
                _sdot,
                gradri2gradR!,
                _grad_len2pos!

const MI = MultiInvariants


import NBodyIPs: ricoords,
                 gradri2gradR!,
                 tdegrees,
                 skip_simplex,
                 fcut,
                 fcut_d,
                 invariants,
                 invariants_ed,
                 gradri2gradR!


# -------------- IO -------------------


MultiDesc(transform::String, cutoff::Union{String, Tuple}, sp_type, valN) =
         MultiDesc(SpaceTransform(transform), Cutoff(cutoff),sp_type, valN)

Dict(D::MultiDesc) = Dict( "__id__"    =>  "MultiDesc",
                                "transform" =>  Dict(D.transform),
                                "cutoff"    =>  Dict(D.cutoff),
                                "sp_type"    =>  Dict(D.sp_type),
                                "valN"    =>  Dict(D.valN) )

MultiDesc(D::Dict) = MultiDesc( SpaceTransform(D["transform"]),
                                          Cutoff(D["cutoff"]),
                                          D["sp_type"],
                                          D["valN"] )

==(D1::MultiDesc, D2::MultiDesc) =
      ( (D1.transform == D2.transform) && (D1.cutoff == D2.cutoff)
        && (D1.sp_type == D2.sp_type) && (D1.valN == D2.valN) )

Base.convert(::Val{:MultiDesc}, D::Dict) = MultiDesc(D)

# ------------- Interface Code ---------------

tdegrees(::MultiDesc, vN::Val{N}) where {N} = MI.tdegrees(vN)

@inline skip_simplex(D::MultiDesc, r) = (maximum(r) > cutoff(D.cutoff))

@inline fcut(D::MultiDesc, r) = fcut(D.cutoff, r)
@inline fcut_d(D::MultiDesc, r) = fcut_d(D.cutoff, r)

@inline invariants(D::MultiDesc, r) = MI.invariants(transform.(D, r),D.sp_type)

@inline function invariants_ed(D::MultiDesc, r)
   x = transform.(D, r)
   I1, I2, DI1, DI2 = MI.invariants_ed(x,D.sp_type)
   x_d = transform_d.(D, r)
   return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
end

# ------------- 3B ricoords and gradri2gradR! ---------------
@inline ricoords(D::MultiDesc, Rs, J) = edge_lengths(Rs, J)

@inline gradri2gradR!(desc::MultiDesc, dVsite, dV_dr, Rs, J, r) =
   _grad_len2pos!(dVsite, dV_dr, Rs, J, r)

# @inline ricoords(D::MultiDesc, Rs, J, Spi, Spj, ::Val{:AAA}) = edge_lengths(Rs, J)
#
# function ricoords(D::MultiDesc, Rs, J, Spi, Spj, ::Val{:AAB})
#    r = edge_lengths(Rs, J)
#    if Spj[J[1]] == Spj[J[2]]
#       return r
#    elseif Spi == Spj[J[1]]
#       return [r[2],r[3],r[1]]
#    elseif Spi == Spj[J[2]]
#       return [r[1],r[3],r[2]]
#    else
#       error("error in ricoords for MultiDesc")
#    end
# end
#
# @inline ricoords(D::MultiDesc, Rs, J, ::Val{:ABC}) = edge_lengths(Rs, J)
# @inline ricoords(D::MultiDesc, Rs, J, ::Val{:AAB}) = lengths_and_angles(Rs, J)



#
# @inline gradri2gradR!(desc::MultiDesc, dVsite, dV_drθ, Rs, J, rθ) =
#    _grad_rθ2pos!(dVsite, dV_drθ, Rs, J, rθ...)
#
# # the cut-off for the bond-angle descriptor depends only on r but not on θ
# @inline fcut(D::MultiDesc, rθ) = fcut(D.cutoff, rθ[1])
#
# @inline function fcut_d(D::BondAngleDesc, rθ)
#    fc, fc_d = fcut_d(D.cutoff, rθ[1])
#    return fc, vcat(fc_d, zero(typeof(rθ[2])))
# end
#
# @inline skip_simplex(D::BondAngleDesc, rθ) = (maximum(rθ[1]) > cutoff(D.cutoff))
#
# @inline _rθ2x(D, r, θ) = vcat(transform.(D, r), θ)
# @inline _rθ2x_d(D, r, θ::SVector{K}) where {K} =
#    vcat(transform_d.(D, r), @SVector ones(K))
#
# @inline invariants(D::BondAngleDesc, rθ) = BAI.invariants(_rθ2x(D, rθ...))
#
# @inline function invariants_ed(D::BondAngleDesc, rθ)
#    x = _rθ2x(D, rθ...)
#    I1, I2, DI1, DI2 = BAI.invariants_ed(x)
#    x_d = _rθ2x_d(D, rθ...)
#    return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
# end


# -------------- Kernel Functions --------------

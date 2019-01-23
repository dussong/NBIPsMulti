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
                _grad_len2pos!,
                lengths_and_angles,
                _rθ2x,
                _rθ2x_d,
                _grad_rθ2pos!

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
to_str(::Val{T}) where {T} = string(T)
get_val(v::Val{T}) where {T} = T

MultiDesc(transform::String, cutoff::Union{String, Tuple}, sp_type, valN) =
         MultiDesc(SpaceTransform(transform), Cutoff(cutoff),sp_type, valN)

MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:AA}) =
         MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                   Val(:AA), Val(2))

MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:AAA}) =
         MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                   Val(:AAA), Val(3))

MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:AAB}) =
         MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                   Val(:AAB), Val(3))

MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:ABC}) =
         MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                   Val(:ABC), Val(3))

# 4-Body

MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:AAAA}) =
                            MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                                      Val(:AAAA), Val(4))


MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:AAAB}) =
                           MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                                     Val(:AAAB), Val(4))


MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:AABB}) =
                          MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                                    Val(:AABB), Val(4))


MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:AABC}) =
                         MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                                   Val(:AABC), Val(4))


MultiDesc(transform::String, cutoff::Union{String, Tuple}, ::Val{:ABCD}) =
                           MultiDesc(SpaceTransform(transform), Cutoff(cutoff),
                                     Val(:ABCD), Val(4))



Dict(D::MultiDesc) = Dict( "__id__"    =>  "MultiDesc",
                            "transform" =>  Dict(D.transform),
                            "cutoff"    =>  Dict(D.cutoff),
                            "sp_type"    =>  to_str(D.sp_type),
                            "valN"    =>  get_val(D.valN)
                             )

MultiDesc(D::Dict) = MultiDesc( SpaceTransform(D["transform"]),
                                          Cutoff(D["cutoff"]),
                                          Val(Symbol(D["sp_type"])),
                                          Val(D["valN"])
                                          )

==(D1::MultiDesc, D2::MultiDesc) =
      ( (D1.transform == D2.transform) && (D1.cutoff == D2.cutoff)
        && (D1.sp_type == D2.sp_type) && (D1.valN == D2.valN) )

Base.convert(::Val{:MultiDesc}, D::Dict) = MultiDesc(D)

# ------------- Interface Code ---------------

tdegrees(::MultiDesc, vN::Val{N}) where {N} = MI.tdegrees(vN)

# ------------- fcut function ---------------

# 2-body
@inline fcut(D::MultiDesc, ::Val{:AA}, r) = fcut(D.cutoff, r)

# 3-body
@inline fcut(D::MultiDesc, ::Val{:AAA}, r) = fcut(D.cutoff, r)
@inline fcut(D::MultiDesc, ::Val{:AAB}, rθ) = fcut(D.cutoff, rθ[1])
@inline fcut(D::MultiDesc, ::Val{:ABC}, r) = fcut(D.cutoff, r)

# 4-body
@inline fcut(D::MultiDesc, ::Val{:AAAA}, r) = fcut(D.cutoff, r)
@inline fcut(D::MultiDesc, ::Val{:AAAB}, rθ) = fcut(D.cutoff, rθ[1])
@inline fcut(D::MultiDesc, ::Val{:AABB}, r) = fcut(D.cutoff, r)
@inline fcut(D::MultiDesc, ::Val{:AABC}, r) = fcut(D.cutoff, r)
@inline fcut(D::MultiDesc, ::Val{:ABCD}, r) = fcut(D.cutoff, r)

# wrap-up
@inline fcut(D::MultiDesc, r) = fcut(D::MultiDesc, D.sp_type, r)


# ------------- fcut_d function ---------------

# 2-body
@inline fcut_d(D::MultiDesc, ::Val{:AA}, r) = fcut_d(D.cutoff, r)

# 3-body
@inline fcut_d(D::MultiDesc, ::Val{:AAA}, r) = fcut_d(D.cutoff, r)
@inline function fcut_d(D::MultiDesc, ::Val{:AAB}, rθ)
   fc, fc_d = fcut_d(D.cutoff, rθ[1])
   return fc, vcat(fc_d, zero(typeof(rθ[2])))
end
@inline fcut_d(D::MultiDesc, ::Val{:ABC}, r) = fcut_d(D.cutoff, r)

# 4-body
@inline fcut_d(D::MultiDesc, ::Val{:AAAA}, r) = fcut_d(D.cutoff, r)
@inline function fcut_d(D::MultiDesc, ::Val{:AAAB}, rθ)
   fc, fc_d = fcut_d(D.cutoff, rθ[1])
   return fc, vcat(fc_d, zero(typeof(rθ[2])))
end
@inline fcut_d(D::MultiDesc, ::Val{:AABB}, r) = fcut_d(D.cutoff, r)
@inline fcut_d(D::MultiDesc, ::Val{:AABC}, r) = fcut_d(D.cutoff, r)
@inline fcut_d(D::MultiDesc, ::Val{:ABCD}, r) = fcut_d(D.cutoff, r)

# wrap-up
@inline fcut_d(D::MultiDesc, r) = fcut_d(D, D.sp_type, r)


# ------------- skip_simplex function ---------------
# 2-body
@inline skip_simplex(D::MultiDesc, ::Val{:AA}, r) = (maximum(r) > cutoff(D.cutoff))

# 3-body
@inline skip_simplex(D::MultiDesc, ::Val{:AAA}, r) = (maximum(r) > cutoff(D.cutoff))
@inline skip_simplex(D::MultiDesc, ::Val{:AAB}, rθ) = (maximum(rθ[1]) > cutoff(D.cutoff))
@inline skip_simplex(D::MultiDesc, ::Val{:ABC}, r) = (maximum(r) > cutoff(D.cutoff))

# wrap-up
@inline skip_simplex(D::MultiDesc, r) = skip_simplex(D, D.sp_type, r)


# ------------- invariants function ---------------
# 2-body
@inline invariants(D::MultiDesc, ::Val{:AA}, r) = MI.invariants(transform.(D, r),D.sp_type)

# 3-body
@inline invariants(D::MultiDesc, ::Val{:AAA}, r) = MI.invariants(transform.(D, r),D.sp_type)
@inline invariants(D::MultiDesc, ::Val{:AAB}, rθ) = MI.invariants(_rθ2x(D, rθ...),D.sp_type)
@inline invariants(D::MultiDesc, ::Val{:ABC}, r) = MI.invariants(transform.(D, r),D.sp_type)

# wrap-up
@inline invariants(D::MultiDesc, r) = invariants(D, D.sp_type, r)



# ------------- invariants_ed function ---------------
# 2-body
@inline function invariants_ed(D::MultiDesc, ::Val{:AA}, r)
   x = transform.(D, r)
   I1, I2, DI1, DI2 = MI.invariants_ed(x,D.sp_type)
   x_d = transform_d.(D, r)
   return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
end

# 3-body
@inline function invariants_ed(D::MultiDesc, ::Val{:AAA}, r)
   x = transform.(D, r)
   I1, I2, DI1, DI2 = MI.invariants_ed(x,D.sp_type)
   x_d = transform_d.(D, r)
   return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
end

@inline function invariants_ed(D::MultiDesc, ::Val{:AAB}, rθ)
   x = _rθ2x(D, rθ...)
   I1, I2, DI1, DI2 = MI.invariants_ed(x,D.sp_type)
   x_d = _rθ2x_d(D, rθ...)
   return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
end

@inline function invariants_ed(D::MultiDesc, ::Val{:ABC}, r)
   x = transform.(D, r)
   I1, I2, DI1, DI2 = MI.invariants_ed(x,D.sp_type)
   x_d = transform_d.(D, r)
   return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
end

# wrap-up
@inline invariants_ed(D::MultiDesc, r) = invariants_ed(D, D.sp_type, r)


# ------------- ricoords function ---------------
# 2-body
@inline ricoords(::Val{:AA}, Rs, J) = edge_lengths(Rs, J)

# 3-body
@inline ricoords(::Val{:AAA}, Rs, J) = edge_lengths(Rs, J)
@inline ricoords(::Val{:AAB}, Rs, J) = lengths_and_angles(Rs, J)
@inline ricoords(::Val{:ABC}, Rs, J) = edge_lengths(Rs, J)

# 4-body
@inline ricoords(::Val{:AAAA}, Rs, J) = edge_lengths(Rs, J)
@inline ricoords(::Val{:AAAB}, Rs, J) = lengths_and_angles(Rs, J)
@inline ricoords(::Val{:AABB}, Rs, J) = edge_lengths(Rs, J)
@inline ricoords(::Val{:AABC}, Rs, J) = edge_lengths(Rs, J)
@inline ricoords(::Val{:ABCD}, Rs, J) = edge_lengths(Rs, J)

# wrap-up
@inline ricoords(D::MultiDesc, Rs, J) = ricoords(D.sp_type, Rs, J)


# ------------- gradri2gradR! function ---------------
# 2-body
@inline gradri2gradR!(::Val{:AA}, dVsite, dV_dr, Rs, J, r) =
   _grad_len2pos!(dVsite, dV_dr, Rs, J, r)

# 3-body
@inline gradri2gradR!(::Val{:AAA}, dVsite, dV_dr, Rs, J, r) =
_grad_len2pos!(dVsite, dV_dr, Rs, J, r)

@inline gradri2gradR!(::Val{:AAB}, dVsite, dV_drθ, Rs, J, rθ) =
_grad_rθ2pos!(dVsite, dV_drθ, Rs, J, rθ...)

@inline gradri2gradR!(::Val{:ABC}, dVsite, dV_dr, Rs, J, r) =
_grad_len2pos!(dVsite, dV_dr, Rs, J, r)

# wrap-up
@inline gradri2gradR!(D::MultiDesc, dVsite, dV_dr, Rs, J, r) = gradri2gradR!(D.sp_type, dVsite, dV_dr, Rs, J, r)


# -------------- Kernel Functions --------------

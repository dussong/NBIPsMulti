using JuLIP: JVec

const MI = MultiInvariants


# -------------- IO -------------------


MultiDesc(transform::String, cutoff::Union{String, Tuple}, sp_type) =
         MultiDesc(SpaceTransform(transform), Cutoff(cutoff),sp_type)

Dict(D::MultiDesc) = Dict( "__id__"    =>  "MultiDesc",
                                "transform" =>  Dict(D.transform),
                                "cutoff"    =>  Dict(D.cutoff),
                                "sp_type"    =>  Dict(D.sp_type) )

MultiDesc(D::Dict) = MultiDesc( SpaceTransform(D["transform"]),
                                          Cutoff(D["cutoff"]),
                                          D["sp_type"] )

==(D1::MultiDesc, D2::MultiDesc) =
      ( (D1.transform == D2.transform) && (D1.cutoff == D2.cutoff)
        && (D1.sp_type == D2.sp_type) )

Base.convert(::Val{:MultiDesc}, D::Dict) = MultiDesc(D)

# ------------- Interface Code ---------------

tdegrees(::MultiDesc, vN::Val{N}) where {N} = MI.tdegrees(vN)


# ------------- 3B ricoords and gradri2gradR! ---------------
@inline ricoords(D::MultiDesc, Rs, J, ::Val{:AAA}) = edge_lengths(Rs, J)
@inline ricoords(D::MultiDesc, Rs, J, ::Val{:ABC}) = edge_lengths(Rs, J)
@inline ricoords(D::MultiDesc, Rs, J, ::Val{:AAB}) = lengths_and_angles(Rs, J)

@inline gradri2gradR!(desc::MultiDesc, dVsite, dV_drθ, Rs, J, rθ) =
   _grad_rθ2pos!(dVsite, dV_drθ, Rs, J, rθ...)

# the cut-off for the bond-angle descriptor depends only on r but not on θ
@inline fcut(D::MultiDesc, rθ) = fcut(D.cutoff, rθ[1])

@inline function fcut_d(D::BondAngleDesc, rθ)
   fc, fc_d = fcut_d(D.cutoff, rθ[1])
   return fc, vcat(fc_d, zero(typeof(rθ[2])))
end

@inline skip_simplex(D::BondAngleDesc, rθ) = (maximum(rθ[1]) > cutoff(D.cutoff))

@inline _rθ2x(D, r, θ) = vcat(transform.(D, r), θ)
@inline _rθ2x_d(D, r, θ::SVector{K}) where {K} =
   vcat(transform_d.(D, r), @SVector ones(K))

@inline invariants(D::BondAngleDesc, rθ) = BAI.invariants(_rθ2x(D, rθ...))

@inline function invariants_ed(D::BondAngleDesc, rθ)
   x = _rθ2x(D, rθ...)
   I1, I2, DI1, DI2 = BAI.invariants_ed(x)
   x_d = _rθ2x_d(D, rθ...)
   return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
end


# -------------- Kernel Functions --------------

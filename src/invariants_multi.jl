module MultiInvariants

using StaticArrays

using NBodyIPs
using NBodyIPs: invariants, invariants_d, invariants_ed, tdegrees

import NBodyIPs: invariants, invariants_d, invariants_ed, tdegrees




# ------------------------------------------------------------------------
#             2-BODY Invariants
#             fully equivalent to BL/BA invariants
# ------------------------------------------------------------------------

# x = (r12,)

invariants(x::SVector{1, T},::Val{:AA}) where {T} =
      copy(x),
      SVector{1, T}(1.0)

invariants_d(x::SVector{1, T},::Val{:AA}) where {T} =
      (@SVector [ SVector(one(T))  ]),
      (@SVector [ SVector(zero(T)) ])

invariants_ed(x::SVector{1,T},::Val{:AA}) where {T} =
      copy(x),
      SVector{1, T}(1.0),
      (@SVector [ SVector(one(T))  ]),
      (@SVector [ SVector(zero(T)) ])

tdegrees(::Val{:AA}) = (1,), (0,)


# ------------------------------------------------------------------------
#             3-BODY Invariants
# ------------------------------------------------------------------------

# Case :AAA (three identical species), using bond-lengths invariants

# r = (r12, r13, r23)

# the 1.0 is a "secondary invariant"
invariants(r::SVector{3, T},::Val{:AAA}) where {T} =
      (@SVector T[ r[1]+r[2]+r[3],
                   r[1]*r[2] + r[1]*r[3] + r[2]*r[3],
                   r[1]*r[2]*r[3] ]),
      (@SVector T[ 1.0 ])


invariants_d(r::SVector{3, T},::Val{:AAA}) where {T} =
      (@SVector [ (@SVector [1.0, 1.0, 1.0]),
                  (@SVector [r[2]+r[3], r[1]+r[3], r[1]+r[2]]),
                  (@SVector [r[2] * r[3], r[1] * r[3], r[1] * r[2]]) ]),
      (@SVector [ (@SVector [0.0, 0.0, 0.0]) ])

invariants_ed(r::SVector{3, T},::Val{:AAA}) where {T} =
      (@SVector T[ r[1]+r[2]+r[3],
                   r[1]*r[2] + r[1]*r[3] + r[2]*r[3],
                   r[1]*r[2]*r[3] ]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector [1.0, 1.0, 1.0]),
                  (@SVector [r[2]+r[3], r[1]+r[3], r[1]+r[2]]),
                  (@SVector [r[2] * r[3], r[1] * r[3], r[1] * r[2]]) ]),
      (@SVector [ (@SVector [0.0, 0.0, 0.0]) ])

tdegrees(::Val{:AAA}) = (1, 2, 3), (0,)

corners(::Val{:AAA}) = ( SVector(1,2), SVector(1,3), SVector(2,3) )


# Case :AAA (one species), using bond-angle invariants
invariants(x::SVector{3, T},::Val{:AAAba}) where {T} =
      (@SVector T[ x[2] + x[3], x[2] * x[3], x[1] ]),
      (@SVector T[ 1.0 ])


invariants_d(x::SVector{3, T},::Val{:AAAba}) where {T} =
      (@SVector [ (@SVector T[0.0, 1.0, 1.0]),
                  (@SVector T[0.0, x[3], x[2]]),
                  (@SVector T[1.0, 0.0, 0.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])

invariants_ed(x::SVector{3, T},::Val{:AAAba}) where {T} =
      (@SVector T[ x[2] + x[3], x[2] * x[3], x[1] ]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector T[0.0, 1.0, 1.0]),
                  (@SVector T[0.0, x[3], x[2]]),
                  (@SVector T[1.0, 0.0, 0.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])


tdegrees(::Val{:AAAba}) = (1, 2, 1), (0,)


# Case :AAB (two different species), using bond-angle invariants

# x = (r1, r2, θ12) where 0 is the B species, 1 and 2 are A.

# the 1.0 is a "secondary invariant"
invariants(x::SVector{3, T},::Val{:AABba}) where {T} =
      (@SVector T[ x[2] + x[3], x[2] * x[3], x[1] ]),
      (@SVector T[ 1.0 ])


invariants_d(x::SVector{3, T},::Val{:AABba}) where {T} =
      (@SVector [ (@SVector T[0.0, 1.0, 1.0]),
                  (@SVector T[0.0, x[3], x[2]]),
                  (@SVector T[1.0, 0.0, 0.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])

invariants_ed(x::SVector{3, T},::Val{:AABba}) where {T} =
      (@SVector T[ x[2] + x[3], x[2] * x[3], x[1] ]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector T[0.0, 1.0, 1.0]),
                  (@SVector T[0.0, x[3], x[2]]),
                  (@SVector T[1.0, 0.0, 0.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])


tdegrees(::Val{:AABba}) = (1, 2, 1), (0,)


# Case :AAB (two different species), using bond-length variables

# x = (r12, r13, r23) where 1 is the B species, 2 and 3 are A.

# the 1.0 is a "secondary invariant"
invariants(x::SVector{3, T},::Val{:AAB}) where {T} =
      (@SVector T[ x[2] + x[3], x[2] * x[3], x[1] ]),
      (@SVector T[ 1.0 ])


invariants_d(x::SVector{3, T},::Val{:AAB}) where {T} =
      (@SVector [ (@SVector T[0.0, 1.0, 1.0]),
                  (@SVector T[0.0, x[3], x[2]]),
                  (@SVector T[1.0, 0.0, 0.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])

invariants_ed(x::SVector{3, T},::Val{:AAB}) where {T} =
      (@SVector T[ x[2] + x[3], x[2] * x[3], x[1] ]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector T[0.0, 1.0, 1.0]),
                  (@SVector T[0.0, x[3], x[2]]),
                  (@SVector T[1.0, 0.0, 0.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])


tdegrees(::Val{:AAB}) = (1, 2, 1), (0,)


   # Case :ABC (two different species), no real invariants since there is no symmetry

   # x = (r1, r2, r3)

   # the 1.0 is a "secondary invariant"
   invariants(x::SVector{3, T},::Val{:ABC}) where {T} =
         (@SVector T[ x[1], x[2], x[3] ]),
         (@SVector T[ 1.0 ])


   invariants_d(x::SVector{3, T},::Val{:ABC}) where {T} =
         (@SVector [ (@SVector T[1.0, 0.0, 0.0]),
                     (@SVector T[0.0, 1.0, 0.0]),
                     (@SVector T[0.0, 0.0, 1.0]) ]),
         (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])

   invariants_ed(x::SVector{3, T},::Val{:ABC}) where {T} =
         (@SVector T[ x[1], x[2], x[3] ]),
         (@SVector T[ 1.0 ]),
         (@SVector [ (@SVector T[1.0, 0.0, 0.0]),
                     (@SVector T[0.0, 1.0, 0.0]),
                     (@SVector T[0.0, 0.0, 1.0]) ]),
         (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])


   tdegrees(::Val{:ABC}) = (1, 1, 1), (0,)


   # Case :ABC (three different species), no real invariants since there is no symmetry

   # the 1.0 is a "secondary invariant"
   invariants(x::SVector{3, T},::Val{:ABCba}) where {T} =
         (@SVector T[ x[1], x[2], x[3] ]),
         (@SVector T[ 1.0 ])


   invariants_d(x::SVector{3, T},::Val{:ABCba}) where {T} =
         (@SVector [ (@SVector T[1.0, 0.0, 0.0]),
                     (@SVector T[0.0, 1.0, 0.0]),
                     (@SVector T[0.0, 0.0, 1.0]) ]),
         (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])

   invariants_ed(x::SVector{3, T},::Val{:ABCba}) where {T} =
         (@SVector T[ x[1], x[2], x[3] ]),
         (@SVector T[ 1.0 ]),
         (@SVector [ (@SVector T[1.0, 0.0, 0.0]),
                     (@SVector T[0.0, 1.0, 0.0]),
                     (@SVector T[0.0, 0.0, 1.0]) ]),
         (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])


   tdegrees(::Val{:ABCba}) = (1, 1, 1), (0,)



# ------------------------------------------------------------------------
#             4-BODY Invariants
# ------------------------------------------------------------------------

# Case :AAAA (four identical species), using bond-lengths invariants

invariants(x::SVector{6, T},::Val{:AAAA}) where {T} = NBodyIPs.BLInvariants.invariants(x)

invariants_d(x::SVector{6, T},::Val{:AAAA}) where {T} =  NBodyIPs.BLInvariants.invariants_d(x)

invariants_ed(x::SVector{6, T},::Val{:AAAA}) where {T} = NBodyIPs.BLInvariants.invariants_ed(x)

tdegrees(::Val{:AAAA}) = NBodyIPs.BLInvariants.tdegrees(Val(4))


# Case :AAAA (3+1 atoms), using bond-angle invariants

invariants(x::SVector{6, T},::Val{:AAAAba}) where {T} = NBodyIPs.BAInvariants.invariants(x)

invariants_d(x::SVector{6, T},::Val{:AAAAba}) where {T} =  NBodyIPs.BAInvariants.invariants_d(x)

invariants_ed(x::SVector{6, T},::Val{:AAAAba}) where {T} = NBodyIPs.BAInvariants.invariants_ed(x)

tdegrees(::Val{:AAAAba}) = NBodyIPs.BAInvariants.tdegrees(Val(4))




# Case :AAAB (3+1 atoms), using bond-angle invariants

invariants(x::SVector{6, T},::Val{:AAABba}) where {T} = NBodyIPs.BAInvariants.invariants(x)

invariants_d(x::SVector{6, T},::Val{:AAABba}) where {T} =  NBodyIPs.BAInvariants.invariants_d(x)

invariants_ed(x::SVector{6, T},::Val{:AAABba}) where {T} = NBodyIPs.BAInvariants.invariants_ed(x)

tdegrees(::Val{:AAABba}) = NBodyIPs.BAInvariants.tdegrees(Val(4))


# Case AAAB with bond-length

include("NB_4B_AAAB_BL_invariants.jl")
@inline invariants(x::SVector{6},::Val{:AAAB}) = NB_4B_AAAB_BL.invariants_gen(x)
@inline invariants_d(x::SVector{6},::Val{:AAAB}) = NB_4B_AAAB_BL.invariants_d_gen(x)
@inline invariants_ed(x::SVector{6},::Val{:AAAB}) = NB_4B_AAAB_BL.invariants_ed_gen(x)

tdegrees(::Val{:AAAB}) = (1, 1, 2, 2, 3, 3,), (0, 2, 3, 3, 4, 6,)



# Case AABB with bond-lengths

include("NB_4B_AABB_invariants.jl")
@inline invariants(x::SVector{6},::Val{:AABB}) = NB_4B_AABB.invariants_gen(x)
@inline invariants_d(x::SVector{6},::Val{:AABB}) = NB_4B_AABB.invariants_d_gen(x)
@inline invariants_ed(x::SVector{6},::Val{:AABB}) = NB_4B_AABB.invariants_ed_gen(x)

tdegrees(::Val{:AABB}) = (1, 1, 1, 2, 2, 2,), (0, 3,)

# Case AABB with bond-angles

include("NB_4B_AABBba_invariants.jl")
@inline invariants(x::SVector{6},::Val{:AABBba}) = NB_4B_AABBba.invariants_gen(x)
@inline invariants_d(x::SVector{6},::Val{:AABBba}) = NB_4B_AABBba.invariants_d_gen(x)
@inline invariants_ed(x::SVector{6},::Val{:AABBba}) = NB_4B_AABBba.invariants_ed_gen(x)

tdegrees(::Val{:AABBba}) = (1, 1, 1, 2, 2, 2,), (0, 3,)




# Case AABC

include("NB_4B_AABC_invariants.jl")
@inline invariants(x::SVector{6},::Val{:AABC}) = NB_4B_AABC.invariants_gen(x)
@inline invariants_d(x::SVector{6},::Val{:AABC}) = NB_4B_AABC.invariants_d_gen(x)
@inline invariants_ed(x::SVector{6},::Val{:AABC}) = NB_4B_AABC.invariants_ed_gen(x)

tdegrees(::Val{:AABC}) = (1, 1, 1, 1, 2, 2,), (0, 2,)

# Case AABC with bond-angles, only trivial invariants

# the 1.0 is a "secondary invariant"
invariants(x::SVector{6, T},::Val{:AABCba}) where {T} =
      (@SVector T[ x[1], x[2], x[3], x[4], x[5], x[6]]),
      (@SVector T[ 1.0 ])


invariants_d(x::SVector{6, T},::Val{:AABCba}) where {T} =
      (@SVector [ (@SVector T[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
                   ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ])

invariants_ed(x::SVector{6, T},::Val{:AABCba}) where {T} =
      (@SVector T[ x[1], x[2], x[3], x[4], x[5], x[6]]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector T[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
                   ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ])

tdegrees(::Val{:AABCba}) = (1, 1, 1, 1, 1, 1,), (0,)


# Case :ABCD (4 different atoms), only trivial invariants

# x = (r1, r2, r3)

# the 1.0 is a "secondary invariant"
invariants(x::SVector{6, T},::Val{:ABCD}) where {T} =
      (@SVector T[ x[1], x[2], x[3], x[4], x[5], x[6]]),
      (@SVector T[ 1.0 ])


invariants_d(x::SVector{6, T},::Val{:ABCD}) where {T} =
      (@SVector [ (@SVector T[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
                   ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ])

invariants_ed(x::SVector{6, T},::Val{:ABCD}) where {T} =
      (@SVector T[ x[1], x[2], x[3], x[4], x[5], x[6]]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector T[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
                   ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ])

tdegrees(::Val{:ABCD}) = (1, 1, 1, 1, 1, 1,), (0,)



# Case :ABCDba (4 different atoms), only trivial invariants


# the 1.0 is a "secondary invariant"
invariants(x::SVector{6, T},::Val{:ABCDba}) where {T} =
      (@SVector T[ x[1], x[2], x[3], x[4], x[5], x[6]]),
      (@SVector T[ 1.0 ])


invariants_d(x::SVector{6, T},::Val{:ABCDba}) where {T} =
      (@SVector [ (@SVector T[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
                   ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ])

invariants_ed(x::SVector{6, T},::Val{:ABCDba}) where {T} =
      (@SVector T[ x[1], x[2], x[3], x[4], x[5], x[6]]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector T[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
                  (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
                   ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ])

tdegrees(::Val{:ABCDba}) = (1, 1, 1, 1, 1, 1,), (0,)

end

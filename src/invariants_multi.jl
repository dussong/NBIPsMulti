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


# Case :AAB (two different species), using bond-angle invariants

# x = (r1, r2, θ12) where 0 is the B species, 1 and 2 are A.

# the 1.0 is a "secondary invariant"
invariants(x::SVector{3, T},::Val{:AAB}) where {T} =
      (@SVector T[ x[1] + x[2], x[1] * x[2], x[3] ]),
      (@SVector T[ 1.0 ])


invariants_d(x::SVector{3, T},::Val{:AAB}) where {T} =
      (@SVector [ (@SVector T[1.0, 1.0, 0.0]),
                  (@SVector T[x[2], x[1], 0.0]),
                  (@SVector T[0.0, 0.0, 1.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])

invariants_ed(x::SVector{3, T},::Val{:AAB}) where {T} =
      (@SVector T[ x[1] + x[2], x[1] * x[2], x[3] ]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector T[1.0, 1.0, 0.0]),
                  (@SVector T[x[2], x[1], 0.0]),
                  (@SVector T[0.0, 0.0, 1.0]) ]),
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


# ------------------------------------------------------------------------
#             4-BODY Invariants
# ------------------------------------------------------------------------

# Case :AAAA (four identical species), using bond-lengths invariants

invariants(x::SVector{6, T},::Val{:AAAA}) where {T} = NBodyIPs.BLInvariants.invariants(x)

invariants_d(x::SVector{6, T},::Val{:AAAA}) where {T} =  NBodyIPs.BLInvariants.invariants_d(x)

invariants_ed(x::SVector{6, T},::Val{:AAAA}) where {T} = NBodyIPs.BLInvariants.invariants_ed(x)

tdegrees(::Val{:AAAA}) = NBodyIPs.BLInvariants.tdegrees(Val(4))

# Case :AAAB (3+1 atoms), using bond-angle invariants

invariants(x::SVector{6, T},::Val{:AAAB}) where {T} = NBodyIPs.BAInvariants.invariants(x)

invariants_d(x::SVector{6, T},::Val{:AAAB}) where {T} =  NBodyIPs.BAInvariants.invariants_d(x)

invariants_ed(x::SVector{6, T},::Val{:AAAB}) where {T} = NBodyIPs.BAInvariants.invariants_ed(x)

tdegrees(::Val{:AAAB}) = NBodyIPs.BAInvariants.tdegrees(Val(4))

# Case AABB

include("NB_4B_AABB_invariants.jl")
@inline invariants(x::SVector{6},::Val{:AABB}) = NB_4B_AABB.invariants_gen(x)
@inline invariants_d(x::SVector{6},::Val{:AABB}) = NB_4B_AABB.invariants_d_gen(x)
@inline invariants_ed(x::SVector{6},::Val{:AABB}) = NB_4B_AABB.invariants_ed_gen(x)

# Case AABC

include("NB_4B_AABC_invariants.jl")
@inline invariants(x::SVector{6},::Val{:AABC}) = NB_4B_AABC.invariants_gen(x)
@inline invariants_d(x::SVector{6},::Val{:AABC}) = NB_4B_AABC.invariants_d_gen(x)
@inline invariants_ed(x::SVector{6},::Val{:AABC}) = NB_4B_AABC.invariants_ed_gen(x)


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



end

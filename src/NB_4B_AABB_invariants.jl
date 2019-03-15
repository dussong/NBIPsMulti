module NB_4B_AABB

using NBodyIPs.FastPolys
using StaticArrays

import NBodyIPs.tdegrees

const G_NB_4B_AABB = [
[ 1, 2, 3, 4, 5, 6 ]
,[ 1, 4, 5, 2, 3, 6 ]
,[ 1, 3, 2, 5, 4, 6 ]
,[ 1, 5, 4, 3, 2, 6 ]
,]
simplex_permutations(x::SVector{6}) = [x[G_NB_4B_AABB[i]] for i=1:4]
# Primary invariants for NB_4B_AABB
 # : definitions at the beginning of the file
const P1_1 = (1,)

const P2_1 = (2,4,3,5,)

const P3_1 = (6,)

const P4_1 = (2,4,)
const P4_2 = (3,5,)

const P5_1 = (2,3,)
const P5_2 = (4,5,)

const P6_1 = (2,4,3,5,)

# Irreducible secondaries for group NB_4B_AABB
 # : definitions at the beginning of the file
const IS1_1 = (2,4,3,5,)


# Primary invariants for NB_4B_AABB
 # : definitions of the types at the beginning of the file
const P1 = Val((P1_1,))
const P2 = Val((P2_1,))
const P3 = Val((P3_1,))
const P4 = Val((P4_1,P4_2,))
const P5 = Val((P5_1,P5_2,))
const P6 = Val((P6_1,))
# Irreducible secondaries for group NB_4B_AABB
 # : definitions of the types at the beginning of the file
const IS1 = Val((IS1_1,))


function invariants_gen(x1::SVector{6, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NB_4B_AABB
 # : what goes in the function for the evaluation
P1 = fpoly((x1,) , NB_4B_AABB.P1)
P2 = fpoly((x1,) , NB_4B_AABB.P2)
P3 = fpoly((x1,) , NB_4B_AABB.P3)
P4 = fpoly((x1,x1,) , NB_4B_AABB.P4)
P5 = fpoly((x1,x1,) , NB_4B_AABB.P5)
P6 = fpoly((x2,) , NB_4B_AABB.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B_AABB
 # : what goes in the function for the evaluation
IS1 = fpoly((x3,) , NB_4B_AABB.IS1)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1
SEC2  = IS1


return (@SVector [P1,P2,P3,P4,P5,P6,]), (@SVector [SEC1,SEC2,])
 end



function invariants_d_gen(x1::SVector{6, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1

   dx1 = @SVector ones(6)
   dx2 = 2 * x1
   dx3 = 3 * x2
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NB_4B_AABB
 # : what goes in the function for the derivatives
dP1 = fpoly_d((x1,),(dx1,) , NB_4B_AABB.P1)
dP2 = fpoly_d((x1,),(dx1,) , NB_4B_AABB.P2)
dP3 = fpoly_d((x1,),(dx1,) , NB_4B_AABB.P3)
dP4 = fpoly_d((x1,x1,),(dx1,dx1,) , NB_4B_AABB.P4)
dP5 = fpoly_d((x1,x1,),(dx1,dx1,) , NB_4B_AABB.P5)
dP6 = fpoly_d((x2,),(dx2,) , NB_4B_AABB.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B_AABB
 # : what goes in the function for the evaluation
IS1 = fpoly((x3,) , NB_4B_AABB.IS1)


# Irreducible secondaries for group NB_4B_AABB
 # : what goes in the function for the derivatives
dIS1 = fpoly_d((x3,),(dx3,) , NB_4B_AABB.IS1)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


dSEC1   = @SVector zeros(6)
dSEC2  = dIS1


return (dP1,dP2,dP3,dP4,dP5,dP6,), (dSEC1,dSEC2,)
 end



function invariants_ed_gen(x1::SVector{6, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1

   dx1 = @SVector ones(6)
   dx2 = 2 * x1
   dx3 = 3 * x2
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NB_4B_AABB
 # : what goes in the function for the evaluation and derivatives
P1, dP1 = fpoly_ed((x1,),(dx1,) , NB_4B_AABB.P1)
P2, dP2 = fpoly_ed((x1,),(dx1,) , NB_4B_AABB.P2)
P3, dP3 = fpoly_ed((x1,),(dx1,) , NB_4B_AABB.P3)
P4, dP4 = fpoly_ed((x1,x1,),(dx1,dx1,) , NB_4B_AABB.P4)
P5, dP5 = fpoly_ed((x1,x1,),(dx1,dx1,) , NB_4B_AABB.P5)
P6, dP6 = fpoly_ed((x2,),(dx2,) , NB_4B_AABB.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B_AABB
 # : what goes in the function for the evaluation and derivatives
IS1, dIS1 = fpoly_ed((x3,),(dx3,) , NB_4B_AABB.IS1)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1
SEC2  = IS1


dSEC1   = @SVector zeros(6)
dSEC2  = dIS1


return (@SVector [P1,P2,P3,P4,P5,P6,]), (@SVector [SEC1,SEC2,]), (@SVector [dP1,dP2,dP3,dP4,dP5,dP6,]), (@SVector [dSEC1,dSEC2,])
 end

end

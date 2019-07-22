module NB_4B_AAAB_BL

using NBodyIPs.FastPolys
using StaticArrays

import NBodyIPs.tdegrees

const G_NB_4B_AAAB_BL = [
[ 1, 2, 3, 4, 5, 6 ]
,[ 2, 3, 1, 6, 4, 5 ]
,[ 3, 1, 2, 5, 6, 4 ]
,[ 2, 1, 3, 4, 6, 5 ]
,[ 3, 2, 1, 6, 5, 4 ]
,[ 1, 3, 2, 5, 4, 6 ]
,]
simplex_permutations(x::SVector{6}) = [x[G_NB_4B_AAAB_BL[i]] for i=1:6]
# Primary invariants for NB_4B_AAAB_BL
 # : definitions at the beginning of the file
const P1_1 = (1,3,2,)

const P2_1 = (4,5,6,)

const P3_1 = (1,3,2,)

const P4_1 = (4,5,6,)

const P5_1 = (1,3,2,)

const P6_1 = (4,5,6,)

# Irreducible secondaries for group NB_4B_AAAB_BL
 # : definitions at the beginning of the file
const IS1_1 = (1,3,2,)
const IS1_2 = (6,4,5,)

const IS2_1 = (1,1,2,1,2,1,)
const IS2_2 = (2,3,3,2,3,3,)
const IS2_3 = (5,6,4,6,5,4,)

const IS3_1 = (1,3,2,2,3,1,)
const IS3_2 = (4,4,5,4,4,5,)
const IS3_3 = (6,5,6,5,6,6,)


# Primary invariants for NB_4B_AAAB_BL
 # : definitions of the types at the beginning of the file
const P1 = Val((P1_1,))
const P2 = Val((P2_1,))
const P3 = Val((P3_1,))
const P4 = Val((P4_1,))
const P5 = Val((P5_1,))
const P6 = Val((P6_1,))
# Irreducible secondaries for group NB_4B_AAAB_BL
 # : definitions of the types at the beginning of the file
const IS1 = Val((IS1_1,IS1_2,))
const IS2 = Val((IS2_1,IS2_2,IS2_3,))
const IS3 = Val((IS3_1,IS3_2,IS3_3,))


function invariants_gen(x1::SVector{6, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NB_4B_AAAB_BL
 # : what goes in the function for the evaluation
P1 = fpoly((x1,) , NB_4B_AAAB_BL.P1)
P2 = fpoly((x1,) , NB_4B_AAAB_BL.P2)
P3 = fpoly((x2,) , NB_4B_AAAB_BL.P3)
P4 = fpoly((x2,) , NB_4B_AAAB_BL.P4)
P5 = fpoly((x3,) , NB_4B_AAAB_BL.P5)
P6 = fpoly((x3,) , NB_4B_AAAB_BL.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B_AAAB_BL
 # : what goes in the function for the evaluation
IS1 = fpoly((x1,x1,) , NB_4B_AAAB_BL.IS1)
IS2 = fpoly((x1,x1,x1,) , NB_4B_AAAB_BL.IS2)
IS3 = fpoly((x1,x1,x1,) , NB_4B_AAAB_BL.IS3)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1
SEC2  = IS1
SEC3  = IS2
SEC4  = IS3
SEC5  = IS1^2
SEC6  = IS1^3


return (@SVector [P1,P2,P3,P4,P5,P6,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,])
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

# Primary invariants for NB_4B_AAAB_BL
 # : what goes in the function for the derivatives
dP1 = fpoly_d((x1,),(dx1,) , NB_4B_AAAB_BL.P1)
dP2 = fpoly_d((x1,),(dx1,) , NB_4B_AAAB_BL.P2)
dP3 = fpoly_d((x2,),(dx2,) , NB_4B_AAAB_BL.P3)
dP4 = fpoly_d((x2,),(dx2,) , NB_4B_AAAB_BL.P4)
dP5 = fpoly_d((x3,),(dx3,) , NB_4B_AAAB_BL.P5)
dP6 = fpoly_d((x3,),(dx3,) , NB_4B_AAAB_BL.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B_AAAB_BL
 # : what goes in the function for the evaluation
IS1 = fpoly((x1,x1,) , NB_4B_AAAB_BL.IS1)
IS2 = fpoly((x1,x1,x1,) , NB_4B_AAAB_BL.IS2)
IS3 = fpoly((x1,x1,x1,) , NB_4B_AAAB_BL.IS3)


# Irreducible secondaries for group NB_4B_AAAB_BL
 # : what goes in the function for the derivatives
dIS1 = fpoly_d((x1,x1,),(dx1,dx1,) , NB_4B_AAAB_BL.IS1)
dIS2 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , NB_4B_AAAB_BL.IS2)
dIS3 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , NB_4B_AAAB_BL.IS3)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


dSEC1   = @SVector zeros(6)
dSEC2  = dIS1
dSEC3  = dIS2
dSEC4  = dIS3
dSEC5  =  + dIS1*2IS1
dSEC6  =  + dIS1*3 * IS1 ^ 2


return (dP1,dP2,dP3,dP4,dP5,dP6,), (dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,)
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

# Primary invariants for NB_4B_AAAB_BL
 # : what goes in the function for the evaluation and derivatives
P1, dP1 = fpoly_ed((x1,),(dx1,) , NB_4B_AAAB_BL.P1)
P2, dP2 = fpoly_ed((x1,),(dx1,) , NB_4B_AAAB_BL.P2)
P3, dP3 = fpoly_ed((x2,),(dx2,) , NB_4B_AAAB_BL.P3)
P4, dP4 = fpoly_ed((x2,),(dx2,) , NB_4B_AAAB_BL.P4)
P5, dP5 = fpoly_ed((x3,),(dx3,) , NB_4B_AAAB_BL.P5)
P6, dP6 = fpoly_ed((x3,),(dx3,) , NB_4B_AAAB_BL.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B_AAAB_BL
 # : what goes in the function for the evaluation and derivatives
IS1, dIS1 = fpoly_ed((x1,x1,),(dx1,dx1,) , NB_4B_AAAB_BL.IS1)
IS2, dIS2 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , NB_4B_AAAB_BL.IS2)
IS3, dIS3 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , NB_4B_AAAB_BL.IS3)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1
SEC2  = IS1
SEC3  = IS2
SEC4  = IS3
SEC5  = IS1^2
SEC6  = IS1^3


dSEC1   = @SVector zeros(6)
dSEC2  = dIS1
dSEC3  = dIS2
dSEC4  = dIS3
dSEC5  =  + dIS1*2IS1
dSEC6  =  + dIS1*3 * IS1 ^ 2


return (@SVector [P1,P2,P3,P4,P5,P6,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,]), (@SVector [dP1,dP2,dP3,dP4,dP5,dP6,]), (@SVector [dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,])
 end

end

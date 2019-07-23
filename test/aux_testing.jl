using Combinatorics, StaticArrays


# 3-Body

const G_3B_AAA = [
[ 1, 2, 3 ]
,[ 1, 3, 2 ]
,[ 2, 1, 3 ]
,[ 2, 3, 1 ]
,[ 3, 1, 2 ]
,[ 3, 2, 1 ]
,]
simplex_permutations(x::SVector{3},::Val{:AAA}) = [x[G_3B_AAA[i]] for i=1:length(G_3B_AAA)]


const G_3B_AAB_BA = [
[ 1, 2, 3 ]
,[ 1, 3, 2 ]
,]
simplex_permutations(x::SVector{3},::Val{:AABba}) = [x[G_3B_AAB_BA[i]] for i=1:length(G_3B_AAB_BA)]

const G_3B_AAB_BL = [
[ 1, 2, 3 ]
,[ 1, 3, 2 ]
,]
simplex_permutations(x::SVector{3},::Val{:AAB}) = [x[G_3B_AAB_BL[i]] for i=1:length(G_3B_AAB_BL)]

const G_3B_ABC = [
[ 1, 2, 3 ]
,]
simplex_permutations(x::SVector{3},::Val{:ABC}) = [x[G_3B_ABC[i]] for i=1:length(G_3B_ABC)]


# 4-Body
const b4_e_inds = [0 1 2 3
                   1 0 4 5
                   2 4 0 6
                   3 5 6 0]

S4_to_S6(π::Vector{Int}) = Int[
   b4_e_inds[π[1], π[2]], b4_e_inds[π[1], π[3]], b4_e_inds[π[1], π[4]],
   b4_e_inds[π[2], π[3]], b4_e_inds[π[2], π[4]], b4_e_inds[π[3], π[4]] ]

simplex_permutations(x::SVector{6},::Val{:AAAA}) =
   (  [ x[S4_to_S6(πX)]
              for πX in permutations(1:4) ]  )


const G_4B_AAAB_BA = [
[ 1, 2, 3, 4, 5, 6 ]
,[ 2, 3, 1, 6, 4, 5 ]
,[ 3, 1, 2, 5, 6, 4 ]
,]
simplex_permutations(x::SVector{6},::Val{:AAABba}) = [x[G_4B_AAAB_BA[i]] for i=1:3]

const G_NB_4B_AAAB_BL = [
[ 1, 2, 3, 4, 5, 6 ]
,[ 2, 3, 1, 6, 4, 5 ]
,[ 3, 1, 2, 5, 6, 4 ]
,[ 2, 1, 3, 4, 6, 5 ]
,[ 3, 2, 1, 6, 5, 4 ]
,[ 1, 3, 2, 5, 4, 6 ]
,]
simplex_permutations(x::SVector{6},::Val{:AAAB}) = [x[G_NB_4B_AAAB_BL[i]] for i=1:6]

const G_NB_4B_AABB = [
[ 1, 2, 3, 4, 5, 6 ]
,[ 1, 4, 5, 2, 3, 6 ]
,[ 1, 3, 2, 5, 4, 6 ]
,[ 1, 5, 4, 3, 2, 6 ]
,]
simplex_permutations(x::SVector{6},::Val{:AABB}) = [x[G_NB_4B_AABB[i]] for i=1:4]



const G_NB_4B_AABC = [
[ 1, 2, 3, 4, 5, 6 ]
,[ 1, 4, 5, 2, 3, 6 ]
,]
simplex_permutations(x::SVector{6},::Val{:AABC}) = [x[G_NB_4B_AABC[i]] for i=1:2]


const G_4B_ABCD = [
[ 1, 2, 3, 4, 5, 6 ],
]
simplex_permutations(x::SVector{6},::Val{:ABCD}) = [x[G_4B_ABCD[i]] for i=1:length(G_4B_ABCD)]

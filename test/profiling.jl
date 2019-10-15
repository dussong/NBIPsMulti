using Profile, ProfileView
using NBodyIPs, NBIPsMulti
# using NBodyIPs: tdegrees
using JuLIP

# println("Setting up the test systems ...")
# r0 = rnn(:Cu)
# TRANSFORM = "r -> exp( - 3 * ((r/$r0) - 1))"
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.05)
# rcut3 = 3.1 * r0
# D3 = BondLengthDesc(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
# rcut4 = 2.1 * r0
# D4 = BondLengthDesc(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )
# rcut5 = 1.5 * r0
# D5 = BondLengthDesc(TRANSFORM, (:cos, 0.66*rcut5, rcut5) )
# DD = [nothing, D3, D3, D4, D5]
#
# # B2 = nbpolys(2, D2, 23)
# # c2 = rand(length(B2))
# B3 = nbpolys(3, D3, 12)
# c3 = rand(length(B3))
# B4 = nbpolys(4, D4, 8)
# c4 = rand(length(B4))
# B = [B3; B4]
# c = [c3; c4]
#
# IP = NBodyIP(B,c)
#
# @time energy(IP,at)
#
# IPf = fast(IP)
# @time energy(IPf,at)



r0 = 1.2*rnn(:Cu)
rcut = 5.0
transform = "exp( - 2 * (r/$r0-1))"
cutoff = "(:cos, $(rcut-1.5), $(rcut))"

# 2B
D2 = MultiDesc(transform,cutoff,Val(:AA))
B2 = nbpolys(D2, 5, [29,30])

# # 3B
D3AAA = MultiDesc(transform,
    cutoff,Val(:AAA))
D3AAB = MultiDesc(transform,
    cutoff,Val(:AAB))
D3ABC = MultiDesc(transform,
    cutoff,Val(:ABC))
B3AAA = nbpolys(D3AAA, 5, [29,29,29])
B3AAB = nbpolys(D3AAB, 5, [29,29,30])
B3ABC = nbpolys(D3ABC, 5, [29,30,31])

D3AAAba = MultiDesc(transform,
    cutoff,Val(:AAAba))
D3AABba = MultiDesc(transform,
    cutoff,Val(:AABba))
D3ABCba = MultiDesc(transform,
    cutoff,Val(:ABCba))
B3AAAba = nbpolys(D3AAAba, 5, [29,29,29])
B3AABba = nbpolys(D3AABba, 5, [29,29,30])
B3ABCba = nbpolys(D3ABCba, 5, [29,30,31])


# # 4B
D4AAAA = MultiDesc(transform,
    cutoff,Val(:AAAA))
D4AAAB = MultiDesc(transform,
    cutoff,Val(:AAAB))
D4AABB = MultiDesc(transform,
    cutoff,Val(:AABB))
D4AABC = MultiDesc(transform,
    cutoff,Val(:AABC))
D4ABCD = MultiDesc(transform,
    cutoff,Val(:ABCD))
B4AAAA = nbpolys(D4AAAA, 5, [29,29,29,29])
B4AAAB = nbpolys(D4AAAB, 5, [29,29,29,30])
B4AABB = nbpolys(D4AABB, 5, [29,29,30,30])
B4AABC = nbpolys(D4AABC, 5, [29,29,30,31])
B4ABCD = nbpolys(D4ABCD, 5, [29,30,31,32])


D4AAAAba = MultiDesc(transform,
    cutoff,Val(:AAAAba))
D4AAABba = MultiDesc(transform,
    cutoff,Val(:AAABba))
D4ABCDba = MultiDesc(transform,
    cutoff,Val(:ABCDba))
B4AAAAba = nbpolys(D4AAAAba, 5, [29,29,29,29])
B4AAABba = nbpolys(D4AAABba, 5, [29,29,29,30])
B4ABCDba = nbpolys(D4ABCDba, 5, [29,30,31,32])


B = [B2;
     B3AAA;
     B3AAB; B3ABC; B3AAAba; B3AABba; B3ABCba;
     B4AAAA;
     B4AAAB; B4AABB; B4AABC; B4ABCD;
     B4AAAAba; B4AAABba; B4ABCDba
     ]
c = rand(length(B))
IPM = NBodyIP(B, c)

@time energy(IPM,at)
@time energy(IPM,at)

IPf = fast(IPM)
@time energy(IPf,at)
@time energy(IPf,at)

22

# function skip_simplex_species!(Sp,Spi,Spj,Species,J)
#    Sp[1] = [Spi]
#    for i=1:length(J)
#       push!(Sp,Spj[J[i]])
#    end
#    return sort(Sp) != sort(Species)
# end
#
# @time skip_simplex_species(3,[4,5],[3,4],SVector(1,2))
#
Profile.clear()
Profile.init(delay = 1e-8)
@profile E = forces(IPM,at)
# @code_warntype energy(IPM,at)
# @
# # @profile dEfast = -forces(IPM, at) |> mat
ProfileView.view()
#
# a = [1]
#
# push!(a,2)

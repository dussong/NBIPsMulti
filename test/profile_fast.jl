# Test energy and forces: 4-body
using JuLIP, NBodyIPs, NBIPsMulti, StaticArrays
using Test
using LinearAlgebra: norm
using Printf
using BenchmarkTools
using Profile, ProfileView

using NBodyIPs: tdegrees, invariants, invariants_d, invariants_ed

at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
Z1 = atomic_numbers(at)
Z1[2] = 30
Z1[3] = 30
Z1[4] = 30
Z1[10] = 30
Z1[11] = 31
Z1[12] = 31
Z1[22] = 31
Z1[14] = 32
Z1[15] = 32
Z1[16] = 32
Z1[17] = 32
at.Z = Z1

at_positions = copy(positions(at)) |> mat

println("generate some basis functions")
B1 = [OneBody(1.0)]

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


# 4B
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
     B3AAA; B3AAB; B3ABC; B3AAAba; B3AABba; B3ABCba;
     B4AAAA; B4AAAB; B4AABB; B4AABC; B4ABCD;
     B4AAAAba; B4AAABba; B4ABCDba
     ]
c = rand(length(B))
IPM = NBodyIP(B, c)
IPMf = fast(IPM)



@time E = energy(IPM, at)
@time Efast = energy(IPMf, at)

@test(abs.(E1-E) < 1e-10)

@time dE = -forces(IPM, at) |> mat
@time dEfast = -forces(IPMf, at) |> mat

@test(norm(dE - dEfast, Inf) < 1e-10)

using Test
using NBodyIPs
using NBodyIPs.Polys
using NBodyIPs: BondLengthDesc, AnalyticTransform, OneBody, BASIS
using JuLIP: save_json, load_json, decode_dict
using NBIPsMulti
##
println("Check (De-)Dictionisation of `MultiDesc`")
D = MultiDesc("exp( - 2 * (r/2.5-1))", "(:cos, 6., 9.)",Val(:AA))
Ds = Dict(D)
D1 = MultiDesc(Ds)
println(@test D1 == D)
println(@test hash(BASIS(), D) == hash(BASIS(), D1))
##
println("generate some basis functions")
rcuts = [9.2, 6.2, 4.5]
TRANSFORM = PolyTransform(3, 2.9)
CUTOFF = [CosCut(0.66*rcut, rcut) for rcut in rcuts]

B1 = [OneBody(1.0)]

r0 = 2.3
transform = "exp( - 2 * (r/$r0-1))"
cutoff = "(:cos, $(r0-1.5), $(r0))"

# 2B
D2 = MultiDesc(transform,cutoff,Val(:AA))
B2 = nbpolys(D2, 5, [1,2])

# 3B
D3AAA = MultiDesc(transform,
    cutoff,Val(:AAA))
D3AAB = MultiDesc(transform,
    cutoff,Val(:AAB))
D3ABC = MultiDesc(transform,
    cutoff,Val(:ABC))
B3AAA = nbpolys(D3AAA, 5, [1,1,1])
B3AAB = nbpolys(D3AAB, 5, [1,1,2])
B3ABC = nbpolys(D3ABC, 5, [1,2,3])

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
B4AAAA = nbpolys(D4AAAA, 5, [1,1,1,1])
B4AAAB = nbpolys(D4AAAB, 5, [1,1,1,2])
B4AABB = nbpolys(D4AABB, 5, [1,1,2,2])
B4AABC = nbpolys(D4AABC, 5, [1,1,2,3])
B4ABCD = nbpolys(D4ABCD, 5, [1,2,3,4])

B = [B2;
     B3AAA; B3AAB; B3ABC
     B4AAAA; B4AAAB; B4AABB; B4AABC; B4ABCD]
c = rand(length(B))
IP = NBodyIP(B, c)
println("Check (De-)Dictionisation of `NBodyIP`")
D_IP = Dict(IP)
IP1 = NBodyIP(D_IP)
println(@test Dict(IP1) == D_IP)
println(@test IP1 == IP)

##
# store it as a JSON file and load it again
fname = tempname() * ".json"
save_ip(fname, IP)
IP2, _ = load_ip(fname)
println(@test IP2 == IP)
run(`rm $fname`)


# Same test with environment
println("same test with environment")
Vn = ("1.0/(1.0 + exp(1.0 / ($(1.5*r0) - r + 1e-2)))", 1.5*r0)

weights = Dict( (1,1)=> 1.0,
                (1,2)=> 2.0,
                (2,1)=> 2.0,
                (1,3)=> 5.0,
                (3,1)=> 5.0,
                (1,4)=> 5.0,
                (4,1)=> 5.0,
                (2,3)=> 5.0,
                (3,2)=> 5.0,
                (2,4)=> 5.0,
                (4,2)=> 5.0,
                (3,4)=> 5.0,
                (4,3)=> 5.0,
                (2,2)=> 5.0,
                (3,3)=> 5.0,
                (4,4)=> 5.0)

# 2B
r0 = 2.3
B2 = envpolysM(D2, 5, Vn, 2, weights, [1,1]);

# 3B
B3AAA = envpolysM(D3AAA, 3, Vn, 1, weights, [1,1,1])
B3AAB = envpolysM(D3AAB, 3, Vn, 1, weights, [1,1,2])
B3ABC = envpolysM(D3ABC, 3, Vn, 1, weights, [1,2,3])

# 4B
B4AAAA = envpolysM(D4AAAA, 3, Vn, 1, weights, [1,1,1,1])
B4AAAB = envpolysM(D4AAAB, 3, Vn, 1, weights, [1,1,1,2])
B4AABB = envpolysM(D4AABB, 3, Vn, 1, weights, [1,1,2,2])
B4AABC = envpolysM(D4AABC, 3, Vn, 1, weights, [1,1,2,3])
B4ABCD = envpolysM(D4ABCD, 3, Vn, 1, weights, [1,2,3,4])

B = [B2;
     B3AAA; B3AAB; B3ABC
     B4AAAA; B4AAAB; B4AABB; B4AABC; B4ABCD]
c = rand(length(B))
IP = NBodyIP(B, c)
println("Check (De-)Dictionisation of `NBodyIP`")
D_IP = Dict(IP)
IP1 = NBodyIP(D_IP)
println(@test Dict(IP1) == D_IP)
println(@test IP1 == IP)

##
# store it as a JSON file and load it again
fname = tempname() * ".json"
save_ip(fname, IP)
IP2, _ = load_ip(fname)
println(@test IP2 == IP)
run(`rm $fname`)

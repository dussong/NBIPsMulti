using NBodyIPs, NBIPsMulti
using JuLIP, Test, StaticArrays, ForwardDiff, Combinatorics
using Printf
using LinearAlgebra: norm

using JuLIP.Potentials: evaluate, evaluate_d

const minvariants = NBIPsMulti.MultiInvariants.invariants
const minvariants_d = NBIPsMulti.MultiInvariants.invariants_d
const minvariants_ed = NBIPsMulti.MultiInvariants.invariants_ed

include("aux_testing.jl")

all_minvariants(r,sp_type) = vcat(minvariants(r,sp_type)...)  # [I1; I2]
ad_minvariants(r,sp_type) = ForwardDiff.jacobian(all_minvariants, r)

println("-------------------------------------------")
println("   Testing implementation of `minvariants`")
println("-------------------------------------------")

# TODO: test correctness of the minvariants implementation
#       against the MAGMA output

println("3-Body ")
println("-------------------------------------------")

Sp_type = [Val(:AAA), Val(:AAB), Val(:ABC),
           Val(:AAAA), Val(:AAAB), Val(:AABB), Val(:AABC), Val(:ABCD)]
dim = [3,3,3,6,6,6,6,6]

for (i,sp_type) in enumerate(Sp_type)
   println(sp_type)
   println("[1] Correctness of gradients")
   r = 1.0 + SVector(rand(dim[i])...)

   I = all_minvariants(r,sp_type)
   dI1, dI2 = minvariants_d(r,sp_type)
   dI = [hcat(dI1...)'; hcat(dI2...)']
   dIh = zeros(size(dI))
   r0 = Vector(r)
   errs = []
   for p = 2:9
      h = .1^p
      dIh = zeros(size(dI))
      for j = 1:length(r)
         r0[j] += h
         Ih = all_minvariants(SVector(r0...),sp_type)
         dIh[:, j] = (Ih - I) / h
         r0[j] -= h
      end
      push!(errs, norm(dIh - dI, Inf))
      @printf(" %d | %.2e \n", p, errs[end])
   end
   println("---------------")
   @test minimum(errs) <= 1e-3 * maximum(errs)
   println()

   println("[2] Symmetry")
   for n = 1:3
      r = 1.0 + SVector(rand(dim[i])...)
      I = all_minvariants(r,sp_type)
      # @show I
      for rπ in simplex_permutations(r,sp_type)
         # @show all_minvariants(SVector(rπ...),sp_type)
         @test I ≈ all_minvariants(SVector(rπ...),sp_type)
      end
      print(".")
   end
   println()

   println("[3] minvariants_ed")
   for n = 1:3
      r = 1.0 + SVector(rand(dim[i])...)
      I1, I2 = minvariants(r,sp_type)
      dI1, dI2 = minvariants_d(r,sp_type)
      J1, J2, dJ1, dJ2 = minvariants_ed(r,sp_type)
      @test all(i ≈ j for (i,j) in zip(I1, J1))
      @test all(i ≈ j for (i,j) in zip(I2, J2))
      @test all(di ≈ dj for (di,dj) in zip(dI1, dJ1))
      @test all(di ≈ dj for (di,dj) in zip(dI2, dJ2))
      print(".")
   end
   println()
end

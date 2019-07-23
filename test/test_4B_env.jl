# Test energy and forces: 4-body+env
using JuLIP, NBodyIPs, NBIPsMulti, StaticArrays
using Test
using LinearAlgebra: norm
using Printf

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


r0 = 1.2*rnn(:Cu)

valSp = [
         Val(:AAAA),
         Val(:AAAB),Val(:AAAB),
         Val(:AAABba),Val(:AAABba),
         Val(:AABB),
         Val(:AABC),
         Val(:ABCD)
         ]
Sp = [
      [29,29,29,29],
      [29,29,29,30],[29,30,30,30],
      [29,29,29,30],[29,30,30,30],
      [29,29,30,30],
      [29,29,30,31],
      [29,30,31,32]
      ]
weights = Dict( (29,29)=> 0.5,
                (29,30)=> 0.1,
                (30,29)=> 0.1,
                (29,31)=> 0.1,
                (31,29)=> 0.1,
                (29,32)=> 0.1,
                (32,29)=> 0.1,
                (30,30) => 0.4,
                (30,31)=> 0.1,
                (31,30)=> 0.1,
                (30,32)=> 0.1,
                (32,30)=> 0.1,
                (31,31)=> 0.1,
                (31,32)=> 0.2,
                (32,31)=> 0.2,
                (32,32)=> 0.2
                )


for i in 1:length(valSp)
   println("-------------------------------------------------")
   println(" Tests $(valSp[i]) ")
   println("-------------------------------------------------")

   BL = MultiDesc("exp( - 2 * (r/$r0-1))",
                    "(:cos, $(r0-1.5), $(r0))",valSp[i])

   Vn = ("exp(- 3 * ((r/$r0)-1))", r0)
   basis = [ envpolysM(BL, 3, Vn, 1, weights,  Sp[i]); ]
   basis[1].Sp
   energy(basis[1],at,basis[1].Sp)
   # site_energies(basis[1], at, [29,29])






   println("-------------------------------------------------")
   println(" Test finite difference energy vs forces - implementation with evaluate ")
   println("-------------------------------------------------")
   E = energy(basis[1], at)
   @show E

   dE = -forces(basis[1], at) |> mat

   errs = []
   # loop through finite-difference step-lengths
   @printf("---------|----------- \n")
   @printf("    p    | error \n")
   @printf("---------|----------- \n")
   for p = 2:9
      h = .1^p
      dEh = zeros(size(at_positions))
      for j = 1:length(at_positions)
         at_positions[j] += h
         set_positions!(at,at_positions)
         Eh = energy(basis[1],at)
         dEh[j] =  (Eh - E) / h
         at_positions[j] -= h
         set_positions!(at,at_positions)
      end
      push!(errs, norm(dEh - dE, Inf))
      @printf(" %d | %.2e \n", p, errs[end])
   end
   println("---------------")
   @test minimum(errs) <= 1e-3 * maximum(errs)



   println("-------------------------------------------------")
   println(" Test finite difference energy vs forces - implementation with evaluate_many ")
   println("-------------------------------------------------")
   E = energy(basis, at)
   F = -forces(basis, at)
   dE = [mat(F[i]) for i = 1:length(F)]

   errs = []
   # loop through finite-difference step-lengths
   @printf("---------|----------- \n")
   @printf("    p    | error \n")
   @printf("---------|----------- \n")
   for p = 2:9
      h = .1^p
      dEh = Matrix{Float64}[]
      for k=1:length(basis)
         push!(dEh,zeros(size(at_positions)))
      end
      for j = 1:length(at_positions)
         at_positions[j] += h
         set_positions!(at,at_positions)
         Eh = energy(basis,at)
         for k = 1:length(basis)
            dEh[k][j] =  (Eh[k] - E[k]) / h
         end
         at_positions[j] -= h
         set_positions!(at,at_positions)
      end
      push!(errs, maximum(norm(dEh[k] - dE[k], Inf) for k=1:length(basis)))
      @printf(" %d | %.2e \n", p, errs[end])
   end
   println("---------------")
   @test minimum(errs) <= 1e-3 * maximum(errs)


   println("-------------------------------------------------")
   println(" Test comparison between 2 implementations - evaluate and evaluate_many ")
   println("-------------------------------------------------")

   println(" Energy difference between 2 implementations - ")
   E1 = energy(basis, at)
   E2 = [energy(basis[k], at) for k=1:length(basis)]
   print(norm(E1-E2,Inf))
   @test norm(E1-E2,Inf) <= 1e-9

   println(" Forces difference between 2 implementations -")
   F = forces(basis,at)
   F1 = [mat(F[i]) for i = 1:length(F)]
   F2 = [mat(forces(basis[k], at)) for k=1:length(basis)]
   print(norm([norm(F1[k]-F2[k],Inf) for k=1:length(basis)], Inf))
   @test norm([norm(F1[k]-F2[k],Inf) for k=1:length(basis)], Inf) <= 1e-10

end

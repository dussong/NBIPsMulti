# Test energy and forces
using JuLIP, NBodyIPs, NBIPsMulti, StaticArrays
using BenchmarkTools, Base.Test


r0 = 2.5
BL2 = BondLengthDesc("exp( - 2 * (r/$r0-1))",
                    "(:cos, $(r0-1.5), $(r0))")

basis = [
     nbpolys(2, BL2, 14, [29,29], Val{:AA});
     nbpolys(2, BL2, 14, [29,30], Val{:AA});
     nbpolys(2, BL2, 14, [30,30], Val{:AA});
  ]

basis[1].Sp
energy(basis[1],at,basis[1].Sp)
site_energies(basis[1], at, [29,29])




at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
Z1 = atomic_numbers(at)
Z1[2] = 30
Z1[10] = 30
at.Z = Z1

at_positions = copy(positions(at)) |> mat

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
   push!(errs, vecnorm(dEh - dE, Inf))
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
   push!(errs, maximum(vecnorm(dEh[k] - dE[k], Inf) for k=1:length(basis)))
   @printf(" %d | %.2e \n", p, errs[end])
end
println("---------------")
@test minimum(errs) <= 1e-3 * maximum(errs)


println("-------------------------------------------------")
println(" Test comparison between 2 implementations - evaluate and evaluate_many ")
println("-------------------------------------------------")

println(" Energy ")
E1 = energy(basis, at)
E2 = [energy(basis[k], at) for k=1:length(basis)]
@test vecnorm(E1-E2,Inf) <= 1e-10

println(" Forces ")
F = forces(basis,at)
F1 = [mat(F[i]) for i = 1:length(F)]
F2 = [mat(forces(basis[k], at)) for k=1:length(basis)]
@test vecnorm([vecnorm(F1[k]-F2[k],Inf) for k=1:length(basis)], Inf) <= 1e-10



# r0 = 2.5
# V = bapolys(2, "($r0/r)^4", "(:cos, 3.6, 4.8)", 4)
# Vcucu = V[3]
#
# Species = [30,30]
#
# # :Cu = 29, :Zn = 30
# at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
# Z1 = atomic_numbers(at)
# Z1[2] = 30
# Z1[10] = 30
# at.Z = Z1
#
# at_positions = copy(positions(at)) |> mat
#
# E = energy(Vcucu, at, Species)
# dE = -forces(Vcucu, at, Species) |> mat
#
# println("-------------------------------------------------")
# println(" Test finite difference energy vs forces ")
# println("-------------------------------------------------")
#
# errs = []
# # loop through finite-difference step-lengths
# @printf("---------|----------- \n")
# @printf("    p    | error \n")
# @printf("---------|----------- \n")
# for p = 2:9
#    h = .1^p
#    dEh = zeros(size(at_positions))
#    for j = 1:length(at_positions)
#       at_positions[j] += h
#       set_positions!(at,at_positions)
#       Eh = energy(Vcucu,at,Species)
#       # Eh = energy(Vcucu,at)
#       dEh[j] =  (Eh - E) / h
#       at_positions[j] -= h
#       set_positions!(at,at_positions)
#    end
#    push!(errs, vecnorm(dEh - dE, Inf))
#    @printf(" %d | %.2e \n", p, errs[end])
# end
# println("---------------")
# @test minimum(errs) <= 1e-3 * maximum(errs)

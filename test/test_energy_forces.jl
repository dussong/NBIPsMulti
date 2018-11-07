# Test energy and forces
using JuLIP, NBodyIPs, StaticArrays
# NBIPsMulti,
using BenchmarkTools, Base.Test

r0 = 2.5
V = blpolys(2, "($r0/r)^4", "(:cos, 3.6, 4.8)", 4)
Vcucu = V[3]
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
at_positions = copy(positions(at)) |> mat

# E = energy(Vcucu, at, [29,29])
E = energy(Vcucu, at)
# dE = forces(Vcucu, at, [29,29]) |> mat
dE = -forces(Vcucu, at) |> mat

println("`NBPoly` finite-difference test on configurations")
println("------------------------------------------------")
at1 = rattle!(bulk(:Cu, cubic=true) * (1,2,2), 0.02)
# at2 = bulk(:Cu)
# set_constraint!(at2, VariableCell(at2, free = []))
VN = Vcucu
(@test JuLIP.Testing.fdtest(VN, at1)) |> println


println("-------------------------------------------------")
println(" Test finite difference energy vs forces ")
println("-------------------------------------------------")

errs = []
# loop through finite-difference step-lengths
@printf("---------|----------- \n")
@printf("    h    | error \n")
@printf("---------|----------- \n")
for p = 2:9
   h = .1^p
   dEh = zeros(size(at_positions))
   for j = 1:length(at_positions)
      at_positions[j] += h
      set_positions!(at,at_positions)
      # Eh = energy(Vcucu,at,[29,29])
      Eh = energy(Vcucu,at)
      dEh[j] =  (Eh - E) / h
      at_positions[j] -= h
      set_positions!(at,at_positions)
   end
   push!(errs, vecnorm(dEh - dE, Inf))
   @printf(" %d | %.2e \n", p, errs[end])
end
println("---------------")
@test minimum(errs) <= 1e-3 * maximum(errs)

using NBIPsMulti
using JuLIP, NBodyIPs, Test, StaticArrays

@testset "NBodyIPsMulti" begin
   @testset "invariants" begin include("test_invariants.jl") end
   @testset "Energy_forces" begin include("test_energy_forces.jl") end
   @testset "2-body" begin include("test_2B.jl") end
   @testset "3-body" begin include("test_3B.jl") end
   @testset "4-body" begin include("test_4B.jl") end

   @testset "2-body-env" begin include("test_2B_env.jl") end
   @testset "3-body-env" begin include("test_3B_env.jl") end
   @testset "4-body-env" begin include("test_4B_env.jl") end
end

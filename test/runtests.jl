include("../src/NBIPsMulti.jl")

using NBIPsMulti
using JuLIP, NBodyIPs, Base.Test, StaticArrays

@testset "NBodyIPsMulti" begin
   @testset "Energy_forces" begin include("test_energy_forces.jl") end
   @testset "2-body" begin include("test_2B.jl") end
   @testset "3-body" begin include("test_3B.jl") end
end

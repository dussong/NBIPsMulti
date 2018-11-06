using NBodyIPs
using JuLIP, Base.Test, StaticArrays

@testset "NBodyIPsMulti" begin
   @testset "Energy_forces" begin include("test_energy_forces.jl") end
end

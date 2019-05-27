using SubstitutionModels
using BioSymbols
using Test
using LinearAlgebra
using StaticArrays


import SubstitutionModels._π,
       SubstitutionModels.P_generic


function Qtest(mod::NASM)
  x = Q(mod)
  return x[DNA_A, DNA_A] ≈ -(x[DNA_A, DNA_C] +
                              x[DNA_A, DNA_G] +
                              x[DNA_A, DNA_T]) &&
         x[DNA_C, DNA_C] ≈ -(x[DNA_C, DNA_A] +
                              x[DNA_C, DNA_G] +
                              x[DNA_C, DNA_T]) &&
         x[DNA_G, DNA_G] ≈ -(x[DNA_G, DNA_A] +
                              x[DNA_G, DNA_C] +
                              x[DNA_G, DNA_T]) &&
         x[DNA_T, DNA_T] ≈ -(x[DNA_T, DNA_A] +
                              x[DNA_T, DNA_C] +
                              x[DNA_T, DNA_G])
end


function test_scaling(mod::NASM)
  q = Q(mod, true)
  return _π(mod) ⋅ -diag(q) ≈ 1.0
end


"""
Test that P(t) = exp(Q × t) and Q × t = log(P(t))
"""
function test_p_q_relationship(mod::NASM, scaled::Bool)
  q = Q(mod, scaled)
  p_1 = P(mod, 1.0, scaled)
  p_2 = P(mod, 2.0, scaled)

  return p_1            ≈ exp(q) &&
         p_2            ≈ exp(q * 2.0) &&
         q              ≈ real(log(Array(p_1))) &&
         q * 2.0        ≈ real(log(Array(p_2)))
end

@testset "JC69" begin
  testmod1 = JC69()
  testmod2 = JC69(1.0)
  testmod3 = JC69(2.0)
  @test Qtest(testmod1)
  @test Q(testmod1) == Q(testmod2)
  @test test_scaling(testmod1)
  @test test_scaling(testmod2)
  @test test_scaling(testmod3)
  @test test_p_q_relationship(testmod1, true)
  @test test_p_q_relationship(testmod2, true)
  @test test_p_q_relationship(testmod3, true)
  @test test_p_q_relationship(testmod1, false)
  @test test_p_q_relationship(testmod2, false)
  @test test_p_q_relationship(testmod3, false)
  @test P_generic(testmod1, [1.5 3.0])[1] ≈ P(testmod1, 1.5)
  @test P_generic(testmod1, [1.5 3.0])[2] ≈ P(testmod1, 3.0)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
  @test sum(P(testmod1, 2.5), dims=2) ≈ [1 1 1 1]'
  @test sum(P(testmod2, 2.5), dims=2) ≈ [1 1 1 1]'
end


@testset "K80" begin
  testmod1 = K80(0.5)
  testmod2 = K80(0.5, 1.0)
  testmod3 = K80(0.5, 2.0)
  @test Qtest(testmod1)
  @test Q(testmod1) == Q(testmod2)
  @test test_scaling(testmod1)
  @test test_scaling(testmod2)
  @test test_scaling(testmod3)
  @test test_p_q_relationship(testmod1, true)
  @test test_p_q_relationship(testmod2, true)
  @test test_p_q_relationship(testmod3, true)
  @test test_p_q_relationship(testmod1, false)
  @test test_p_q_relationship(testmod2, false)
  @test test_p_q_relationship(testmod3, false)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
  @test sum(P(testmod1, 2.5), dims=2) ≈ [1 1 1 1]'
  @test sum(P(testmod2, 2.5), dims=2) ≈ [1 1 1 1]'
end


@testset "F81" begin
  testmod1 = F81(0.1, 0.2, 0.3, 0.4)
  testmod2 = F81(1.0, 0.1, 0.2, 0.3, 0.4)
  testmod3 = F81(2.0, 0.1, 0.2, 0.3, 0.4)
  @test Qtest(testmod1)
  @test Q(testmod1) == Q(testmod2)
  @test test_scaling(testmod1)
  @test test_scaling(testmod2)
  @test test_scaling(testmod3)
  @test test_p_q_relationship(testmod1, true)
  @test test_p_q_relationship(testmod2, true)
  @test test_p_q_relationship(testmod3, true)
  @test test_p_q_relationship(testmod1, false)
  @test test_p_q_relationship(testmod2, false)
  @test test_p_q_relationship(testmod3, false)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
  @test sum(P(testmod1, 2.5), dims=2) ≈ [1 1 1 1]'
  @test sum(P(testmod2, 2.5), dims=2) ≈ [1 1 1 1]'
end


@testset "F84" begin
  testmod1 = F84(0.75, 0.1, 0.2, 0.3, 0.4)
  testmod2 = F84(0.75, 1.0, 0.1, 0.2, 0.3, 0.4)
  testmod3 = F84(0.75, 2.0, 0.1, 0.2, 0.3, 0.4)
  @test Qtest(testmod1)
  @test Q(testmod1) == Q(testmod2)
  @test test_scaling(testmod1)
  @test test_scaling(testmod2)
  @test test_scaling(testmod3)
  @test test_p_q_relationship(testmod1, true)
  @test test_p_q_relationship(testmod2, true)
  @test test_p_q_relationship(testmod3, true)
  @test test_p_q_relationship(testmod1, false)
  @test test_p_q_relationship(testmod2, false)
  @test test_p_q_relationship(testmod3, false)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
  @test sum(P(testmod1, 2.5), dims=2) ≈ [1 1 1 1]'
  @test sum(P(testmod2, 2.5), dims=2) ≈ [1 1 1 1]'
end


@testset "HKY85" begin
  testmod1 = HKY85(0.75, 0.1, 0.2, 0.3, 0.4)
  testmod2 = HKY85(0.75, 1.0, 0.1, 0.2, 0.3, 0.4)
  testmod3 = HKY85(0.75, 2.0, 0.1, 0.2, 0.3, 0.4)
  @test Qtest(testmod1)
  @test Q(testmod1) == Q(testmod2)
  @test test_scaling(testmod1)
  @test test_scaling(testmod2)
  @test test_scaling(testmod3)
  @test test_p_q_relationship(testmod1, true)
  @test test_p_q_relationship(testmod2, true)
  @test test_p_q_relationship(testmod3, true)
  @test test_p_q_relationship(testmod1, false)
  @test test_p_q_relationship(testmod2, false)
  @test test_p_q_relationship(testmod3, false)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
  @test sum(P(testmod1, 2.5), dims=2) ≈ [1 1 1 1]'
  @test sum(P(testmod2, 2.5), dims=2) ≈ [1 1 1 1]'
end


@testset "TN93" begin
  testmod1 = TN93(0.6, 0.7, 0.1, 0.2, 0.3, 0.4)
  testmod2 = TN93(0.6, 0.7, 1.0, 0.1, 0.2, 0.3, 0.4)
  testmod3 = TN93(0.6, 0.7, 2.0, 0.1, 0.2, 0.3, 0.4)
  @test Qtest(testmod1)
  @test Q(testmod1) == Q(testmod2)
  @test test_scaling(testmod1)
  @test test_scaling(testmod2)
  @test test_scaling(testmod3)
  @test test_p_q_relationship(testmod1, true)
  @test test_p_q_relationship(testmod2, true)
  @test test_p_q_relationship(testmod3, true)
  @test test_p_q_relationship(testmod1, false)
  @test test_p_q_relationship(testmod2, false)
  @test test_p_q_relationship(testmod3, false)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
  @test sum(P(testmod1, 2.5), dims=2) ≈ [1 1 1 1]'
  @test sum(P(testmod2, 2.5), dims=2) ≈ [1 1 1 1]'
end


@testset "GTR" begin
  testmod1 = GTR(1.0, 2.0, 3.0, 4.0, 5.0, 0.1, 0.2, 0.3, 0.4)
  testmod2 = GTR(1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 0.1, 0.2, 0.3, 0.4)
  testmod3 = GTR(1.0, 2.0, 3.0, 4.0, 5.0, 2.0, 0.1, 0.2, 0.3, 0.4)
  @test Qtest(testmod1)
  @test Q(testmod1) == Q(testmod2)
  @test test_scaling(testmod1)
  @test test_scaling(testmod2)
  @test test_scaling(testmod3)
  @test test_p_q_relationship(testmod1, true)
  @test test_p_q_relationship(testmod2, true)
  @test test_p_q_relationship(testmod3, true)
  @test test_p_q_relationship(testmod1, false)
  @test test_p_q_relationship(testmod2, false)
  @test test_p_q_relationship(testmod3, false)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test P(testmod1, [1.0 2.0])[1] ≈ P_generic(testmod1, 1.0)
  @test P(testmod2, [1.0 2.0])[1] ≈ P_generic(testmod2, 1.0)
  @test P(testmod1, [1.0 2.0])[2] ≈ P_generic(testmod1, 2.0)
  @test P(testmod2, [1.0 2.0])[2] ≈ P_generic(testmod2, 2.0)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
  @test sum(P(testmod1, 2.5), dims=2) ≈ [1 1 1 1]'
  @test sum(P(testmod2, 2.5), dims=2) ≈ [1 1 1 1]'
end

@testset "Indexing" begin
  testmat = MMatrix{4, 4,Int64}(1:16)
  @test testmat[DNA_A] == testmat[DNA_A, DNA_A]
  testmat[DNA_G] = 100
  @test testmat[DNA_G] == 100
end

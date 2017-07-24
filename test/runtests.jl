using SubstitutionModels
using Base.Test

@testset "JC69" begin
  testmod1 = JC69()
  testmod2 = JC69(1.0)
  @test Q(testmod1) == Q(testmod2)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
end

@testset "K80" begin
  testmod1 = K80(0.5)
  testmod2 = K80(0.5, 1.0)
  @test Q(testmod1) == Q(testmod2)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
end

@testset "F81" begin
  testmod1 = F81(0.3, 0.2, 0.2, 0.3)
  testmod2 = F81(1.0, 0.3, 0.2, 0.2, 0.3)
  @test Q(testmod1) == Q(testmod2)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
end

@testset "F84" begin
  testmod1 = F84(0.75, 0.3, 0.2, 0.2, 0.3)
  testmod2 = F84(0.75, 1.0, 0.3, 0.2, 0.2, 0.3)
  @test Q(testmod1) == Q(testmod2)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
end

@testset "HKY85" begin
  testmod1 = HKY85(0.75, 0.3, 0.2, 0.2, 0.3)
  testmod2 = HKY85(0.75, 1.0, 0.3, 0.2, 0.2, 0.3)
  @test Q(testmod1) == Q(testmod2)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
end

@testset "TN93" begin
  testmod1 = TN93(0.6, 0.7, 0.3, 0.2)
  testmod2 = TN93(0.6, 0.7, 1.0, 0.3, 0.2)
  @test Q(testmod1) == Q(testmod2)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
end

@testset "GTR" begin
  testmod1 = GTR(1.0, 2.0, 3.0, 4.0, 5.0, 0.3, 0.2, 0.2, 0.3)
  testmod2 = GTR(1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 0.3, 0.2, 0.2, 0.3)
  @test Q(testmod1) == Q(testmod2)
  @test P(testmod1, 1.0e2) ≈ P_generic(testmod1, 1.0e2)
  @test P(testmod2, 1.0e2) ≈ P_generic(testmod2, 1.0e2)
  @test isapprox(diag(P(testmod1, 1.0e9)), _π(testmod1), atol = 1.0e-5)
  @test isapprox(diag(P(testmod2, 1.0e9)), _π(testmod2), atol = 1.0e-5)
  @test sum(_π(testmod1)) == 1.0
  @test sum(_π(testmod2)) == 1.0
end

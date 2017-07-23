using SubstitutionModels2
using Base.Test

@testset "JC69" begin
  testmod = JC69()
  @test Q(testmod) == Q(JC69(1.0))
  @test isapprox(diag(P(testmod, 1.0e10)), _π(testmod), atol=1.0e-5)
  @test sum(_π(testmod)) == 1.0
  @test _πR(testmod) + _πY(testmod) == 1.0
end

@testset "K80" begin
  testmod = K80(0.5)
  @test Q(testmod) == Q(K80(0.5, 1.0))
  @test isapprox(diag(P(testmod, 1.0e10)), _π(testmod), atol=1.0e-5)
  @test sum(_π(testmod)) == 1.0
  @test _πR(testmod) + _πY(testmod) == 1.0
end

@testset "F81" begin
  testmod = F81(0.3, 0.2, 0.2, 0.3)
  @test Q(testmod) == Q(F81(1.0, 0.3, 0.2, 0.2, 0.3))
  @test isapprox(diag(P(testmod, 1.0e10)), _π(testmod), atol=1.0e-5)
  @test sum(_π(testmod)) == 1.0
  @test _πR(testmod) + _πY(testmod) == 1.0
end

@testset "F84" begin
  testmod = F84(0.75, 0.3, 0.2, 0.2, 0.3)
  @test Q(testmod) == Q(F84(0.75, 1.0, 0.3, 0.2, 0.2, 0.3))
  @test isapprox(diag(P(testmod, 1.0e10)), _π(testmod), atol=1.0e-5)
  @test sum(_π(testmod)) == 1.0
  @test _πR(testmod) + _πY(testmod) == 1.0
end

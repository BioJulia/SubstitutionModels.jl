using SubstitutionModels
using BioSymbols
using Test
using LinearAlgebra
using StaticArrays


import SubstitutionModels._π,
       SubstitutionModels._scale,
       SubstitutionModels.scale_generic,
       SubstitutionModels.P_generic


function test_mod_fun(mod::Type{T}, n_params::Int64, equal_base_freqs::Bool, closed_form_p::Bool) where T <: NASM
  _dummy_params = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
  _dummy_freqs = [0.21, 0.29, 0.23, 0.27]
  if equal_base_freqs
    @test_throws ErrorException convert(mod, _dummy_params[1:n_params+1], safe=true)
    if n_params > 0
      @test_throws BoundsError mod(_dummy_params[1:n_params-1], false) # safe = false
    end
    @test_nowarn mod(_dummy_params[1:n_params+1], false) # safe = false
    for i in 1:n_params
      flip = fill(1.0, n_params)
      flip[i] *= -1
      @test_throws ErrorException convert(mod, _dummy_params[1:n_params] .* flip, safe=true)
      @test_nowarn mod(_dummy_params[1:n_params] .* flip, false) # safe = false
    end
    @test_throws MethodError convert(mod, _dummy_params[1:n_params], _dummy_freqs)
    @test_nowarn convert(mod, _dummy_params[1:n_params], safe=true)
    x = mod(_dummy_params[1:n_params])
    @test x == supertype(mod)(_dummy_params[1:n_params]) # Convenience constructor
  else
    @test_throws ErrorException convert(mod, _dummy_params[1:n_params+1], _dummy_freqs)
    if n_params > 0
      @test_throws BoundsError mod(_dummy_params[1:n_params-1], _dummy_freqs, false) # safe = false
    end
    for i in 1:n_params
      flip = fill(1.0, n_params)
      flip[i] *= -1
      @test_throws ErrorException convert(mod, _dummy_params[1:n_params] .* flip, _dummy_freqs)
      @test_nowarn mod(_dummy_params[1:n_params] .* flip, _dummy_freqs, false) # safe = false
    end
    @test_throws ErrorException convert(mod, _dummy_params[1:n_params], _dummy_freqs .+ 0.1)
    @test_throws ErrorException convert(mod, _dummy_params[1:n_params], _dummy_freqs[1:3])
    @test_nowarn mod(_dummy_params[1:n_params], [_dummy_freqs; 0.1], false) # safe = false
    @test_throws BoundsError mod(_dummy_params[1:n_params], _dummy_freqs[1:3], false) # safe = false
    @test_throws MethodError convert(mod, _dummy_params[1:n_params])
    @test_nowarn convert(mod, _dummy_params[1:n_params], _dummy_freqs, safe=true)
    x = mod(_dummy_params[1:n_params], _dummy_freqs)
    @test x == supertype(mod)(_dummy_params[1:n_params], _dummy_freqs) # Convenience constructor
  end
  @test_nowarn Q(x)
  @test_nowarn Q(x, true) # Scaled q matrix
  q1 = Q(x)
  q2 = Q(x, true) # Scaled q matrix
  @test all(.≈(sum(q1, dims=2), 0.0, atol=1e-13)) # Q matrix col sums
  @test all(.≈(sum(q2, dims=2), 0.0, atol=1e-13)) # Scaled Q matrix col sums
  @test _scale(x) ≈ scale_generic(x) # Test specific vs. generic scale method
  @test _π(x) ⋅ -diag(q2) ≈ 1.0 # Consistency of π with scaled Q matrix
  @test_throws ErrorException P(x, -1e3)
  @test P(x, [1e3]) ≈ P_generic(x, [1e3]) # Test P generic function
  @test P(x, [1e3], true) ≈ P_generic(x, [1e3], true) # Test P generic function with scaling
  if closed_form_p # If closed form solution to P matrix calculation
    @test P(x, Inf) ≈ _π(x)' .* [1, 1, 1, 1] # P matrix asymptotics
  end
end


for testmod in [(JC69abs, 1, true, true)
                (JC69rel, 0, true, true)
                (K80abs, 2, true, true)
                (K80rel, 1, true, true)
                (F81abs, 1, false, true)
                (F81rel, 0, false, true)
                (F84abs, 2, false, true)
                (F84rel, 1, false, true)
                (HKY85abs, 2, false, true)
                (HKY85rel, 1, false, true)
                (TN93abs, 3, false, true)
                (TN93rel, 2, false, true)
                (GTRabs, 6, false, false)
                (GTRrel, 5, false, false)]
  @testset "$(testmod[1])" begin
    test_mod_fun(testmod...)
  end
end


@testset "Indexing" begin
  testmat = MMatrix{4, 4,Int64}(1:16)
  @test testmat[DNA_A] == testmat[DNA_A, DNA_A]
  testmat[DNA_G] = 100
  @test testmat[DNA_G] == 100
end

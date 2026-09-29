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
      @test_throws BoundsError mod(_dummy_params[1:n_params-1], safe=false)
    end
    @test_nowarn mod(_dummy_params[1:n_params+1], safe=false)
    for i in 1:n_params
      flip = fill(1.0, n_params)
      flip[i] *= -1
      @test_throws ErrorException convert(mod, _dummy_params[1:n_params] .* flip, safe=true)
      @test_nowarn mod(_dummy_params[1:n_params] .* flip, safe=false)
    end
    @test_throws MethodError convert(mod, _dummy_params[1:n_params], _dummy_freqs)
    @test_nowarn convert(mod, _dummy_params[1:n_params], safe=true)
    x = mod(_dummy_params[1:n_params])
    @test x == supertype(mod)(_dummy_params[1:n_params]) # Convenience constructor
  else
    @test_throws ErrorException convert(mod, _dummy_params[1:n_params+1], _dummy_freqs)
    if n_params > 0
      @test_throws BoundsError mod(_dummy_params[1:n_params-1], _dummy_freqs, safe=false)
    end
    for i in 1:n_params
      flip = fill(1.0, n_params)
      flip[i] *= -1
      @test_throws ErrorException convert(mod, _dummy_params[1:n_params] .* flip, _dummy_freqs)
      @test_nowarn mod(_dummy_params[1:n_params] .* flip, _dummy_freqs, safe=false)
    end
    @test_throws ErrorException convert(mod, _dummy_params[1:n_params], _dummy_freqs .+ 0.1)
    @test_throws ErrorException convert(mod, _dummy_params[1:n_params], _dummy_freqs[1:3])
    @test_nowarn mod(_dummy_params[1:n_params], [_dummy_freqs; 0.1], safe=false)
    @test_throws BoundsError mod(_dummy_params[1:n_params], _dummy_freqs[1:3], safe=false)
    @test_throws MethodError convert(mod, _dummy_params[1:n_params])
    @test_nowarn convert(mod, _dummy_params[1:n_params], _dummy_freqs, safe=true)
    x = mod(_dummy_params[1:n_params], _dummy_freqs)
    @test mod(_dummy_params[1:n_params], view(_dummy_freqs, :)) == x
    @test mod(view(_dummy_params, 1:n_params), _dummy_freqs) == x
    @test convert(mod, _dummy_params[1:n_params], view(_dummy_freqs, :)) == x
    @test x == supertype(mod)(_dummy_params[1:n_params], _dummy_freqs) # Convenience constructor
  end
  @test_nowarn Q(x)
  @test_nowarn Q(x, true) # Scaled q matrix
  @test Q(x) isa SMatrix{4,4,Float64}
  @test P(x, 0.1) isa SMatrix{4,4,Float64}
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


@testset "Nucleotide indexing" begin
  f = [0.21, 0.29, 0.23, 0.27]
  @test f[DNA_A] == 0.21
  @test_throws ArgumentError f[DNA_N]
  f2 = [0.21, 0.29, 0.23]
  @test_throws BoundsError f2[DNA_T]
  sm = reshape(1:16, 4, 4)
  @test [sm[i,j] for i in (RNA_A, RNA_C, RNA_G, RNA_U), j in (RNA_A, RNA_C, RNA_G, RNA_U)] == sm
  @test_throws ArgumentError sm[DNA_G, DNA_N]
  sm2 = reshape(1:9, 3, 3)
  @test_throws BoundsError sm2[DNA_A, DNA_T]
end

const unambiguous_states = ((DNA_A, DNA_C, DNA_G, DNA_T),
                          (RNA_A, RNA_C, RNA_G, RNA_U))
const invalid_states = (DNA_N, RNA_N, DNA_R, RNA_Y, DNA_B, RNA_V,
                        DNA_Gap, RNA_Gap)

@testset "Nucleotide index conversion" begin
  for states in unambiguous_states, (i, nt) in enumerate(states)
    @test nucleotide_index(nt) == i
  end
  for nt in invalid_states
    @test_throws ArgumentError nucleotide_index(nt)
  end
end

@testset "Legacy array indexing" begin
  for a in ([1,2,3,4], SVector(1,2,3,4), MVector(1,2,3,4), view(1:4, :))
    @test [a[nt] for nt in (DNA_A, DNA_C, DNA_G, DNA_T)] == collect(a)
    @test checkbounds(a, RNA_U) === nothing
  end
  for a in ([1,2,3], MVector(1,2,3), view([1,2,3], :))
    @test_throws BoundsError checkbounds(a, DNA_T)
    @test_throws BoundsError setindex!(a, 9, DNA_T)
    @test_throws ArgumentError setindex!(a, 9, DNA_N)
    @test a == [1,2,3]
  end
  a = MMatrix{3,3}(1:9)
  @test_throws BoundsError setindex!(a, 9, DNA_A, DNA_T)
  @test_throws ArgumentError setindex!(a, 9, DNA_N, DNA_A)
  @test a == reshape(1:9, 3,3)
  @test setindex!(a, 42, DNA_C, DNA_G)[2,3] == 42
end

@testset "Legacy deprecation migration" begin
  mktemp() do path, io
    close(io)
    script = "using SubstitutionModels, BioSymbols; a = collect(1:4); a[DNA_A]"
    cmd = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) --startup-file=no --depwarn=yes -e $script`
    run(pipeline(cmd, stderr=path))
    warning = read(path, String)
    @test occursin("deprecated", warning)
    @test occursin("0.6.0", warning)
    @test occursin("nucleotide_index", warning)
  end
  mktemp() do path, io
    close(io)
    script = "using SubstitutionModels; F81([1.0], [0.21, 0.29, 0.23, 0.27])"
    cmd = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) --startup-file=no --depwarn=yes -e $script`
    run(pipeline(cmd, stderr=path))
    @test !occursin("deprecated", read(path, String))
  end
end

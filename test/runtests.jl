using SubstitutionModels
using BioSymbols
using Test
using LinearAlgebra
using StaticArrays


import SubstitutionModels._π,
       SubstitutionModels._scale,
       SubstitutionModels.scale_generic,
       SubstitutionModels.P_generic


struct BrokenFrequencies{E} <: NASM
  failure::E
end
SubstitutionModels.Q(::BrokenFrequencies) = Q(JC69())
SubstitutionModels._π(mod::BrokenFrequencies) = throw(mod.failure)

@testset "Array and model interfaces" begin
  @test isdefined(Main, :NucleotideView)
  for wrap in (identity, a -> view(a, :))
    a = collect(1:4)
    v = SubstitutionModels.NucleotideView(wrap(a))
    @test Base.mightalias(a, v)
    a .= view(v, 4:-1:1)
    @test a == [4,3,2,1]
  end
  v = SubstitutionModels.NucleotideView(reshape(collect(1:16), 4,4))
  @test v[DNA_A, :] == [1,5,9,13]
  @test collect(view(v, :, RNA_G)) == [9,10,11,12]
  @test !occursin(r"[\e\r\n]", repr(JC69()))
  @test !occursin(r"[\e\r]", sprint(show, MIME"text/plain"(), JC69()))
  for err in (ErrorException("frequency bug"), InterruptException(), LinearAlgebra.LAPACKException(1))
    @test_throws typeof(err) P_generic(BrokenFrequencies(err), [0.1])
  end
  mod = F81rel(0.0, 0.3, 0.3, 0.4; safe=false)
  @test P_generic(mod, [0.0, 0.1]) ≈ [exp(Q(mod) * t) for t in [0.0, 0.1]]
end


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
  @test eval(Meta.parse(repr(x))) == x
  @test !occursin(r"[\e\r\n]", repr(x))
  @test !occursin(r"[\e\r]", sprint(show, MIME"text/plain"(), x))
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


@testset "No global nucleotide indexing" begin
  for A in (Vector{Int}, Matrix{Int}), N in (DNA, RNA)
    @test which(getindex, (A, N)).module !== SubstitutionModels
    @test which(getindex, (A, N, N)).module !== SubstitutionModels
    @test which(setindex!, (A, Int, N)).module !== SubstitutionModels
    @test which(setindex!, (A, Int, N, N)).module !== SubstitutionModels
    @test which(checkbounds, (A, N)).module !== SubstitutionModels
    @test which(checkbounds, (A, N, N)).module !== SubstitutionModels
  end
  @test_throws ArgumentError [1, 2, 3, 4][DNA_A]
end

@testset "Nucleotide views" begin
  for a in ([11, 23, 37, 41], SVector{4}(11, 23, 37, 41),
            MVector{4}(11, 23, 37, 41), view([11, 23, 37, 41], :),
            reshape(collect(11:26), 4, 4), SMatrix{4,4}(11:26),
            MMatrix{4,4}(11:26), view(reshape(collect(11:26), 4, 4), :, :))
    v = NucleotideView(a)
    @test parent(v) === a
    @test size(v) == size(a)
    @test axes(v) == axes(a)
    @test length(v) == length(a)
    @test eltype(v) == eltype(a)
    @test collect(v) == a
    @test copy(v) == v
    @test collect(eachindex(v)) == collect(eachindex(a))
    @test v[1] == a[1]
    @test v[:] == a[:]
    for states in ((DNA_A, DNA_C, DNA_G, DNA_T), (RNA_A, RNA_C, RNA_G, RNA_U))
      for (i, nt) in enumerate(states)
        @test nucleotide_index(nt) == i
        @test v[nt] == a[i]
      end
      if ndims(a) == 2
        for (i, ni) in enumerate(states), (j, nj) in enumerate(states)
          @test v[ni, nj] == a[i, j]
        end
        @test v[DNA_A, 3] == a[1, 3]
        @test v[2, RNA_G] == a[2, 3]
        @test_throws BoundsError v[DNA_A, 5]
      end
    end
    for nt in (DNA_N, RNA_N, DNA_R, RNA_Y, DNA_B, RNA_V, DNA_Gap, RNA_Gap)
      @test_throws ArgumentError nucleotide_index(nt)
      @test_throws ArgumentError v[nt]
      @test_throws ArgumentError setindex!(v, 0, nt)
    end
    @test_throws BoundsError v[length(v)+1]
  end
  for a in (ones(3), ones(3,4), ones(4,3), ones(4,4,4))
    @test_throws ArgumentError NucleotideView(a)
  end
  a = [1, 2, 3, 4]
  v = NucleotideView(a)
  v[DNA_G] = 9
  @test a[3] == 9
  c = copy(v)
  c[DNA_A] = 10
  @test a[1] == 1
  resize!(a, 3)
  @test_throws BoundsError v[DNA_T]
  @test_throws BoundsError setindex!(v, 10, DNA_T)
  @test a == [1, 2, 9]
  @test_throws ErrorException setindex!(NucleotideView(SVector(1,2,3,4)), 0, DNA_A)
end

@testset "Explicit nucleotide indices" begin
  f = [0.21, 0.29, 0.23, 0.27]
  @test f[nucleotide_index(DNA_A)] == 0.21
  @test_throws ArgumentError f[nucleotide_index(DNA_N)]
  f2 = [0.21, 0.29, 0.23]
  @test_throws BoundsError f2[nucleotide_index(DNA_T)]
  sm = reshape(1:16, 4, 4)
  @test_throws ArgumentError sm[nucleotide_index(DNA_G), nucleotide_index(DNA_N)]
  sm2 = reshape(1:9, 3, 3)
  @test_throws BoundsError sm2[nucleotide_index(DNA_A), nucleotide_index(DNA_T)]
end

@testset "Nucleotide view writes and bounds" begin
  for a in (reshape(collect(1:16), 4,4), MMatrix{4,4}(1:16),
            view(reshape(collect(1:16), 4,4), :, :))
    v = NucleotideView(a)
    for states in ((DNA_A, DNA_C, DNA_G, DNA_T), (RNA_A, RNA_C, RNA_G, RNA_U))
      for (i, ni) in enumerate(states), (j, nj) in enumerate(states)
        @test setindex!(v, 10i+j, ni, nj) === v
        @test a[i,j] == 10i+j
        @test checkbounds(Bool, v, ni, nj)
      end
    end
    v[DNA_A, 2] = 71
    v[3, RNA_U] = 81
    v[2] = 91
    @test a[1,2] == 71
    @test a[3,4] == 81
    @test a[2] == 91
    @test !checkbounds(Bool, v, DNA_A, 5)
    @test !checkbounds(Bool, v, 5, RNA_A)
    for nt in (DNA_N, RNA_N, DNA_R, RNA_Y, DNA_B, RNA_V, DNA_Gap, RNA_Gap)
      before = copy(a)
      @test_throws ArgumentError v[nt, DNA_A]
      @test_throws ArgumentError v[DNA_A, nt]
      @test_throws ArgumentError setindex!(v, 0, nt, DNA_A)
      @test_throws ArgumentError setindex!(v, 0, DNA_A, nt)
      @test a == before
    end
    before = copy(a)
    @test_throws BoundsError setindex!(v, 0, DNA_A, 5)
    @test a == before
    @test collect(view(v, :, 2)) == a[:,2]
    @test v .+ 1 == a .+ 1
  end
  v = NucleotideView([1,2,3,4])
  @test checkbounds(v, DNA_T) === nothing
  resize!(parent(v), 3)
  @test !checkbounds(Bool, v, DNA_T)
  @test_throws BoundsError checkbounds(v, DNA_T)
end

@testset "Nucleotide vector writes" begin
  for a in ([1,2,3,4], MVector{4}(1,2,3,4), view([1,2,3,4], :))
    v = NucleotideView(a)
    for (i, nt) in enumerate((RNA_A, RNA_C, RNA_G, RNA_U))
      @test setindex!(v, i+10, nt) === v
      @test a[i] == i+10
    end
    before = copy(a)
    @test_throws ArgumentError setindex!(v, 0, RNA_N)
    @test_throws BoundsError setindex!(v, 0, 5)
    @test a == before
  end
  a = view(reshape(collect(1:32), 4,8), :, 1:2:8)
  v = NucleotideView(a)
  @test collect(v) == a
  @test v[DNA_T, DNA_G] == a[4,3]
  @test v[DNA_G] == a[3]
  @test copy(v) == v
  @test v[:,2] == a[:,2]
end

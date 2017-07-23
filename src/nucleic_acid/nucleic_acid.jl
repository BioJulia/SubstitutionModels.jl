_πAT(mod) = _πACGT(mod)
_πGC(mod) = _πACGT(mod)
_πG(mod) = _πGC(mod)
_πC(mod) = _πGC(mod)
_πA(mod) = _πAT(mod)
_πT(mod) = _πAT(mod)
_πACGU(mod) = _πACGT(mod)
_πAU(mod) = _πAT(mod)
_πU(mod) = _πT(mod)

"π = [πA, πC, πG, πT/πU]"
@inline function _π(mod::NASM)
  return SVector(_πA(mod), _πC(mod), _πG(mod), _πT(mod))
end

"πR = πA + πG"
@inline function _πR(mod::NASM)
  return _πA(mod) + _πG(mod)
end

"πY = πT + πC"
@inline function _πY(mod::NASM)
  return _πT(mod) + _πC(mod)
end

"Generate a rate matrix for a `NucleicAcidSubstitutionModel`, of the form:
[[A→A, A→C, A→G, A→T]
 [C→A, C→C, C→G, C→T]
 [G→A, G→C, G→G, G→T]
 [T→A, T→C, T→G, T→T]]

 or

 [[A→A, A→C, A→G, A→U]
  [C→A, C→C, C→G, C→U]
  [G→A, G→C, G→G, G→U]
  [U→A, U→C, U→G, U→U]]"
@inline function _R(mod::NASM)
    α = _α(mod) # α = r(T/U → C) = r(C → T/U)
    β = _β(mod) # β = r(T/U → A) = r(A → T/U)
    γ = _γ(mod) # γ = r(T/U → G) = r(G → T/U)
    δ = _δ(mod) # δ = r(C → A) = r(A → C)
    ϵ = _ϵ(mod) # ϵ = r(C → G) = r(G → C)
    η = _η(mod) # η = r(A → G) = r(G → A)

    return SMatrix{4,4, Float64}(-(δ+η+β), δ, η, β,
                                 δ, -(δ+ϵ+α), ϵ, α,
                                 η, ϵ, -(η+ϵ+γ), γ,
                                 β, α, γ, -(β+α+γ))
end

"Generate a Q matrix for a `NucleicAcidSubstitutionModel`, of the form:
[[A→A, A→C, A→G, A→T]
 [C→A, C→C, C→G, C→T]
 [G→A, G→C, G→G, G→T]
 [T→A, T→C, T→G, T→T]]

 or

 [[A→A, A→C, A→G, A→U]
  [C→A, C→C, C→G, C→U]
  [G→A, G→C, G→G, G→U]
  [U→A, U→C, U→G, U→U]]"
@inline function Q(mod::NASM)
  return SDiagonal(_π(mod)) * _R(mod) * _μ(mod)
end

"Generate a P matrix for a `NucleicAcidSubstitutionModel`, of the form:
[[A→A, A→C, A→G, A→T]
 [C→A, C→C, C→G, C→T]
 [G→A, G→C, G→G, G→T]
 [T→A, T→C, T→G, T→T]]

 or

 [[A→A, A→C, A→G, A→U]
  [C→A, C→C, C→G, C→U]
  [G→A, G→C, G→G, G→U]
  [U→A, U→C, U→G, U→U]]

  for a specified time"
@inline function P(mod::NASM, t::Float64)
  return expm(Q(mod) * t)
end

"Generate an array of P matrices for a specified array of times"
function P(mod::NASM, t::Array{Float64})
  eig_vals, eig_vecs = eig(Q(mod))
  return [expm(eig_vecs * (diagm(eig_vals)*i) * eig_vecs') for i in t]
end

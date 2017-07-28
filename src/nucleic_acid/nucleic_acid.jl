const _πAT(mod) = _πACGT(mod)
const _πCG(mod) = _πACGT(mod)
const _πG(mod) = _πCG(mod)
const _πC(mod) = _πCG(mod)
const _πA(mod) = _πAT(mod)
const _πT(mod) = _πAT(mod)
const _πACGU(mod) = _πACGT(mod)
const _πAU(mod) = _πAT(mod)
const _πU(mod) = _πT(mod)


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


"Generate a Q matrix for a `NucleicAcidSubstitutionModel`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]"
@inline function Q(mod::NASM)
    α = _α(mod)
    β = _β(mod)
    γ = _γ(mod)
    δ = _δ(mod)
    ϵ = _ϵ(mod)
    η = _η(mod)
    πA = _πA(mod)
    πC = _πC(mod)
    πG = _πG(mod)
    πT = _πT(mod)

    return Qmatrix(-(δ*πC+η*πG+β*πT), δ*πA, η*πA, β*πA,
                   δ*πC, -(δ*πA+ϵ*πG+α*πT), ϵ*πC, α*πC,
                   η*πG, ϵ*πG, -(η*πA+ϵ*πC+γ*πT), γ*πG,
                   β*πT, α*πT, γ*πT, -(β*πA+α*πC+γ*πG))
end


"Generate a P matrix for a `NucleicAcidSubstitutionModel`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
@inline function P_generic(mod::NASM, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  return expm(Q(mod) * t)
end


"Generate an array of P matrices for a specified array of times"
function P_generic(mod::NASM, t::Array{Float64})
  if any(t .< 0.0)
    error("t must be positive")
  end
  try
    eig_vals, eig_vecs = eig(Q(mod))
    return [eig_vecs * expm(diagm(eig_vals)*i) * eig_vecs' for i in t]
  catch
    eig_vals, eig_vecs = eig(Array(Q(mod)))
    return [SMatrix(eig_vecs * expm(diagm(eig_vals)*i) * eig_vecs') for i in t]
  end
end


P(mod, t) = P_generic(mod, t)

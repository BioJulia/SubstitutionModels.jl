abstract type F81 <: NASM end

struct F81abs <: F81
  β::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function F81abs(β::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if β <= 0.
      error("F81 parameter β must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("F81 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<=0.0)
      error("F81 frequencies must be positive")
    end
    new(β, πA, πC, πG, πT)
  end
end

struct F81rel <: F81
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function F81rel(πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if sum([πA,πC,πG,πT]) != 1.0
      error("F81 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<=0.0)
      error("F81 frequencies must be positive")
    end
    new(πA, πC, πG, πT)
  end
end

F81(β, πA, πC, πG, πT) = F81abs(β, πA, πC, πG, πT)
F81(πA, πC, πG, πT) = F81rel(πA, πC, πG, πT)

@inline function _μ(mod::F81abs)
  return mod.β
end

@inline function _μ(mod::F81rel)
  return 1.0
end

@inline function _πA(mod::F81)
  return mod.πA
end

@inline function _πC(mod::F81)
  return mod.πC
end

@inline function _πG(mod::F81)
  return mod.πG
end

@inline function _πT(mod::F81)
  return mod.πT
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::F81)
  return 1.0
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::F81)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::F81)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::F81)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::F81)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::F81)
  return 1.0
end

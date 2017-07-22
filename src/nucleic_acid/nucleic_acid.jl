@inline function R(mod::NASM)
    α = _α(mod)
    β = _β(mod)
    γ = _γ(mod)
    δ = _δ(mod)
    ϵ = _ϵ(mod)
    ζ = _ζ(mod)

    return SMatrix{4,4}(-(α+β+γ), α, β, γ,
                        α, -(α+δ+ϵ), δ, ϵ,
                        β, δ, -(β+δ+ζ), ζ,
                        γ, ϵ, ζ, -(γ+ϵ+ζ))
end

@inline function Q(mod::NASM)
  return SDiagonal(π(mod)) * R(mod) * μ(mod)
end

@inline function P(mod::NASM, t::Float64)
  return expm(Q(mod) * t)
end

function P(mod::NASM, t::Array{Float64})
  eig_vals, eig_vecs = eig(Q(mod))
  return [expm(eig_vecs * (diagm(eig_vals)*i) * eig_vecs') for i in t]
end

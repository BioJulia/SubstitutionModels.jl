_πACGT(mod::NASM) = 0.25
_πACGU(mod) = _πACGT(mod)
_πA(mod) = _πACGT(mod)
_πC(mod) = _πACGT(mod)
_πG(mod) = _πACGT(mod)
_πT(mod) = _πACGT(mod)
_π(mod::NASM) = SVector(_πA(mod), _πC(mod), _πG(mod), _πT(mod))
_πR(mod::NASM) = _πA(mod) + _πG(mod)
_πY(mod::NASM) = _πT(mod) + _πC(mod)


"""
Generate a Q matrix for a `NucleicAcidSubstitutionModel`, of the form:

\$Q = \\begin{bmatrix}
Q_{A, A} & Q_{A, C} & Q_{A, G} & Q_{A, T} \\\\
Q_{C, A} & Q_{C, C} & Q_{C, G} & Q_{C, T} \\\\
Q_{G, A} & Q_{G, C} & Q_{G, G} & Q_{G, T} \\\\
Q_{T, A} & Q_{T, C} & Q_{T, G} & Q_{T, T} \\end{bmatrix}\$

Call as either
1) `Q(model)`, or
2) `Q(model, bool)`
Form (2) scales the matrix so that ``_π(model) ⋅ -diag(Q(model)) = 1``
when bool=true. `Q(model, false)` is equivalent to `Q(model)`.
"""
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

    return Qmatrix(-(δ * πC + η * πG + β * πT), δ * πA, η * πA, β * πA,
                   δ * πC, -(δ * πA + ϵ * πG + α * πT), ϵ * πC, α * πC,
                   η * πG, ϵ * πG, -(η * πA + ϵ * πC + γ * πT), γ * πG,
                   β * πT, α * πT, γ * πT, -(β * πA + α * πC + γ * πG))
end


function Q(mod::NASM, scale::Bool)
  if scale
    q = Q(mod)
    scale = 1/(-diag(q)' * _π(mod))
    return q * scale
  else
    Q(mod)
  end
end


scale_generic(mod::NASM) = 1/(-diag(Q(mod))' * _π(mod))


_scale(mod::NASM) = scale_generic(mod)


@inline function P_generic(mod::NASM, t::Float64)
    if t < 0.0
        error("t must be positive")
    end
    return exp(Q(mod) * t)
end


function P_generic(mod::NASM, t, scale::Bool)
    if scale
        P(mod, t * _scale(mod))
    else
        P(mod, t)
    end
end


function P_generic(mod::NASM, t::Array{Float64})
    if any(t .< 0.0)
        error("t must be positive")
    end

    try
        q = Q(mod)

        # Use symmetrical similar matrix if possible:
        # B is similar (has same eigenvalues) to Q, but is symmetrical if the model is reversible.
        # Eigenvalue decomposition is easier with symmetrical matrix
        # NB: computation of B described in Inferring Phylogenies, Felsenstein, p.206

        rootπ = sqrt.(_π(mod))
        b = diagm(0 => rootπ) * q * diagm(0 => 1.0 ./ rootπ)

        # If B is symmetrical, then do eigenvalue decomposition on B. The resulting
        # eigenvalues are the same as for Q, then Q's eigenvectors can be obtained by
        # a simple conversion.
        if ishermitian(round.(b, digits=12))
            eig_vals, r = eigen(b)
            eig_vecs = diagm(0 => 1.0 ./ rootπ) * r # eig_vecs are the eigenvectors of Q
            eig_vecs_inv = r' * diagm(0 => rootπ)

        # If B is not symmetrical, fall back to directly obtaining eigenvalues from Q
        else
            eig_vals, eig_vecs = eigen(Array(q))
            eig_vecs_inv = inv(eig_vecs)
        end

        return [SMatrix{size(q)...}(eig_vecs * diagm(0 => exp.(eig_vals .* i)) * eig_vecs_inv) for i in t]

    catch
        # Any errors, fall back to direct use of matrix exponential
        return [P_generic(mod, i) for i in t]

    end
end


"""
Generate a P matrix for a `NucleicAcidSubstitutionModel`, of the form:

\$P = \\begin{bmatrix}
P_{A, A} & P_{A, C} & P_{A, G} & P_{A, T} \\\\
P_{C, A} & P_{C, C} & P_{C, G} & P_{C, T} \\\\
P_{G, A} & P_{G, C} & P_{G, G} & P_{G, T} \\\\
P_{T, A} & P_{T, C} & P_{T, G} & P_{T, T} \\end{bmatrix}\$

for specified time

Call as either
1) `P(model, t)`, or
2) `P(model, t, bool)`
Form (2) obtains its probabilities from the scaled Q matrix if bool=true.
Branch lengths estimated from a scaled P matrix are in units of expected
number of substitutions per site. `P(model, t, false)` is equivalent to
`P(model, t)`.
"""
P(mod::NASM, t) = P_generic(mod, t)

P(mod::NASM, t, scale) = P_generic(mod, t, scale)

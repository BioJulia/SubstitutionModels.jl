GTR(α::F, β::F, γ::F, δ::F, ϵ::F, η::F,
    πA::F, πC::F, πG::F, πT::F,
    safe::Bool=true) where F <: Float64 =
  GTRabs(α, β, γ, δ, ϵ, η, πA, πC, πG, πT, safe)


GTR(α::F, β::F, γ::F, δ::F, ϵ::F,
    πA::F, πC::F, πG::F, πT::F,
    safe::Bool=true) where F <: Float64 =
  GTRrel(α, β, γ, δ, ϵ, πA, πC, πG, πT, safe)


function GTR(θ_vec::A,
               π_vec::A,
               safe::Bool=true) where A <: AbstractArray
  if safe && length(π_vec) != 4
    error("Incorrect base frequency vector length")
  end
  if length(θ_vec) == 5
    return GTRrel(θ_vec[1], θ_vec[2], θ_vec[3], θ_vec[4], θ_vec[5],
                  π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
  elseif length(θ_vec) == 6
    return GTRabs(θ_vec[1], θ_vec[2], θ_vec[3], θ_vec[4], θ_vec[5], θ_vec[6],
                  π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
  else
    error("Parameter vector length incompatiable with absolute or relative rate form of substitution model")
  end
end


function GTRrel(θ_vec::A,
                π_vec::A,
                safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 5
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return GTRrel(θ_vec[1], θ_vec[2], θ_vec[3], θ_vec[4], θ_vec[5],
                π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                safe)
end


function GTRabs(θ_vec::A,
                π_vec::A,
                safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 6
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return GTRabs(θ_vec[1], θ_vec[2], θ_vec[3], θ_vec[4], θ_vec[5], θ_vec[6],
                π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                safe)
end

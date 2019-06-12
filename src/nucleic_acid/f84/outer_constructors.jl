F84(κ::F,
    πA::F, πC::F, πG::F, πT::F,
    safe::Bool=true) where F <: Float64 =
  F84rel(κ, πA, πC, πG, πT, safe)


F84(κ::F, β::F,
    πA::F, πC::F, πG::F, πT::F,
    safe::Bool=true) where F <: Float64 =
  F84abs(κ, β, πA, πC, πG, πT, safe)


function F84(θ_vec::A,
             π_vec::A,
             safe::Bool=true) where A <: AbstractArray
  if safe && length(π_vec) != 4
    error("Incorrect base frequency vector length")
  end
  if length(θ_vec) == 1
    return F84rel(θ_vec[1],
                  π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
  elseif length(θ_vec) == 2
    return F84abs(θ_vec[1], θ_vec[2],
                  π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
  else
    error("Parameter vector length incompatiable with absolute or relative rate form of substitution model")
  end
end


function F84rel(θ_vec::A,
                π_vec::A,
                safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 1
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return F84rel(θ_vec[1], π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
end


function F84abs(θ_vec::A,
                π_vec::A,
                safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 2
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return F84abs(θ_vec[1], θ_vec[2], π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
end

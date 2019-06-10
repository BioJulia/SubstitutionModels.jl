HKY85(κ::F,
      πA::F, πC::F, πG::F, πT::F,
      safe::Bool=true) where F <: Float64 =
  HKY85rel(κ, πA, πC, πG, πT, safe)


HKY85(α::F, β::F,
      πA::F, πC::F, πG::F, πT::F,
      safe::Bool=true) where F <: Float64 =
  HKY85abs(α, β, πA, πC, πG, πT, safe)


function HKY85(θ_vec::A,
               π_vec::A,
               safe::Bool=true) where A <: AbstractArray
  if safe && length(π_vec) != 4
    error("Incorrect base frequency vector length")
  end
  if length(θ_vec) == 1
    return HKY85rel(θ_vec[1], π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
  elseif length(θ_vec) == 2
    return HKY85abs(θ_vec[1], θ_vec[2], π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
  else
    error("Parameter vector length incompatiable with absolute or relative rate form of substitution model")
  end
end


function HKY85rel(θ_vec::A,
                  π_vec::A,
                  safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 1
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return HKY85rel(θ_vec[1], π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
end


function HKY85abs(θ_vec::A,
                  π_vec::A,
                  safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 2
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return HKY85abs(θ_vec[1], θ_vec[2],
                  π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
end

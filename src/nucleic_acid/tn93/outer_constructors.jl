TN93(α1::F, α2::F, β::F, πA::F, πC::F, πG::F, πT::F, safe::Bool=true) where F <: Float64 = TN93abs(α1, α2, β, πA, πC, πG, πT, safe)


TN93(κ1::F, κ2::F, πA::F, πC::F, πG::F, πT::F, safe::Bool=true) where F <: Float64 = TN93rel(κ1, κ2, πA, πC, πG, πT, safe)


function TN93(θ_vec::A,
              π_vec::A,
              safe::Bool=true) where A <: AbstractArray
  if safe && length(π_vec) != 4
    error("Incorrect base frequency vector length")
  end
  if length(θ_vec) == 2
    return TN93rel(θ_vec[1], θ_vec[2],
                  π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
  elseif length(θ_vec) == 3
    return TN93abs(θ_vec[1], θ_vec[2], θ_vec[3],
                  π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
  else
    error("Parameter vector length incompatiable with absolute or relative rate form of substitution model")
  end
end


function TN93rel(θ_vec::A,
                 π_vec::A,
                 safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 2
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return TN93rel(θ_vec[1], θ_vec[2], π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
end


function TN93abs(θ_vec::A,
                 π_vec::A,
                 safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 3
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return TN93abs(θ_vec[1], θ_vec[2], θ_vec[3], π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
end

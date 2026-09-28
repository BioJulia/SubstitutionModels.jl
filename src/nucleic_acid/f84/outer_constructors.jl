F84(κ::F,
    πA::F, πC::F, πG::F, πT::F;
    safe::Bool=true) where F <: Float64 =
  F84rel(κ, πA, πC, πG, πT, safe=safe)


F84(κ::F, β::F,
    πA::F, πC::F, πG::F, πT::F;
    safe::Bool=true) where F <: Float64 =
  F84abs(κ, β, πA, πC, πG, πT, safe=safe)


function F84(θ_vec::AbstractArray,
             π_vec::AbstractArray;
             safe::Bool=true)
  if safe && length(π_vec) != 4
    error("Incorrect base frequency vector length")
  end
  if length(θ_vec) == 1
    return F84rel(θ_vec[1],
                  π_vec[1], π_vec[2], π_vec[3], π_vec[4],
                  safe=safe)
  elseif length(θ_vec) == 2
    return F84abs(θ_vec[1], θ_vec[2],
                  π_vec[1], π_vec[2], π_vec[3], π_vec[4],
                  safe=safe)
  else
    error("Parameter vector length incompatiable with absolute or relative rate form of substitution model")
  end
end


function F84rel(θ_vec::AbstractArray,
                π_vec::AbstractArray;
                safe::Bool=true)
  if safe
    if length(θ_vec) != 1
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return F84rel(θ_vec[1], π_vec[1], π_vec[2], π_vec[3], π_vec[4], safe=safe)
end


function F84abs(θ_vec::AbstractArray,
                π_vec::AbstractArray;
                safe::Bool=true)
  if safe
    if length(θ_vec) != 2
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return F84abs(θ_vec[1], θ_vec[2], π_vec[1], π_vec[2], π_vec[3], π_vec[4], safe=safe)
end

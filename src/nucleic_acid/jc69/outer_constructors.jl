JC69() = JC69rel()
JC69(λ::F, safe::Bool=true) where F <: Float64 = JC69abs(λ, safe)


function JC69(θ_vec::A,
              safe::Bool=true) where A <: AbstractArray
  if length(θ_vec) == 0
    return JC69rel()
  elseif length(θ_vec) == 1
    return JC69abs(θ_vec[1], safe)
  else
    error("Parameter vector length incompatiable with absolute or relative rate form of substitution model")
  end
end


function JC69rel(θ_vec::A,
                 safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 0
      error("Incorrect parameter vector length")
    end
  end
  return JC69rel()
end


function JC69abs(θ_vec::A,
                 safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 1
      error("Incorrect parameter vector length")
    end
  end
  return  JC69abs(θ_vec[1], safe)
end

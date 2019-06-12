K80(κ::Float64, safe::Bool=true) = K80rel(κ, safe)


K80(α::Float64, β::Float64, safe::Bool=true) = K80abs(α, β, safe)


function K80(θ_vec::A,
              safe::Bool=true) where A <: AbstractArray
  if length(θ_vec) == 1
    return K80rel(θ_vec[1], safe)
  elseif length(θ_vec) == 2
    return K80abs(θ_vec[1], θ_vec[2], safe)
  else
    error("Parameter vector length incompatiable with absolute or relative rate form of substitution model")
  end
end


function K80rel(θ_vec::A,
                safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 1
      error("Incorrect parameter vector length")
    end
  end
  return K80rel(θ_vec[1], safe)
end


function K80abs(θ_vec::A,
                safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 2
      error("Incorrect parameter vector length")
    end
  end
  return K80abs(θ_vec[1], θ_vec[2], safe)
end

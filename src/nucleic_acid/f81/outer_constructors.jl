F81(πA::F, πC::F, πG::F, πT::F,
    safe::Bool=true) where F <: Float64 =
F81rel(πA, πC, πG, πT)


F81(β::F,
    πA::F, πC::F, πG::F, πT::F,
    safe::Bool=true) where F <: Float64 =
F81abs(β, πA, πC, πG, πT, safe)


function F81(θ_vec::A,
             π_vec::A,
             safe::Bool=true) where A <: AbstractArray
  if safe && length(π_vec) != 4
    error("Incorrect base frequency vector length")
  end
  if length(θ_vec) == 0
    return F81rel(π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
  elseif length(θ_vec) == 1
    return F81abs(θ_vec[1],
                  π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T],
                  safe)
  else
    error("Parameter vector length incompatiable with absolute or relative rate form of substitution model")
  end
end


function F81rel(θ_vec::A,
                π_vec::A,
                safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 0
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return F81rel(π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
end



function F81abs(θ_vec::A,
                π_vec::A,
                safe::Bool=true) where A <: AbstractArray
  if safe
    if length(θ_vec) != 1
      error("Incorrect parameter vector length")
    elseif length(π_vec) != 4
      error("Incorrect base frequency vector length")
    end
  end
  return F81abs(θ_vec[1], π_vec[DNA_A], π_vec[DNA_C], π_vec[DNA_G], π_vec[DNA_T], safe)
end

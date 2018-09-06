struct Distance{T <: SM}
  mean::Float64
  variance::Float64
end


function mean(x::Distance)
  return x.mean
end


function variance(x::Distance)
  return x.variance
end


function show(io::IO, x::Distance{T}) where T <: SM
  print(io, "$T distance\nmean = $(mean(x))\nvariance = $(variance(x))")
end

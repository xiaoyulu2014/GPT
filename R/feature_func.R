library(plyr)

feature = function(x,n,sigmaRBF,generator){
  D = length(x)
  set.seed(generator)
  Z = matrix(rnorm(n*D)/sigmaRBF,n,D)
  b = matrix(runif(n*D),n,D)
  phi = sqrt(2/n)*cos(Z * x + b*2*pi)
  return(phi)
}

fea_notensor = function(x,m,sigmaRBF,generator){
  set.seed(generator)
  D = length(x)
  w = matrix(rnorm(m*D)/sigmaRBF,m,D)
  b = runif(m)
  return(sqrt(2/m)*cos(w %*% x + b*2*pi))
}
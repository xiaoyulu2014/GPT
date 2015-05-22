
theta_inf = function(X,y,sigma,n,sigmaRBF,q,generator){
  ##precompute features
  D = ncol(X)
  Phi = apply(X,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,D)),2)) })

  ##initialise Theta
  ##randomly select the set of non-zero entries in W
  I = sort(sample(1:(n^D))[1:(q*n^D)])  
  theta = rep(0,n^D)
  theta[I] = rnorm(length(I))  
  
  ##posterior
  Sigma = solve(1/(sigma^2) * Phi[I,] %*% t(Phi[I,]) + diag(length(I)))
  Mu = 1/(sigma^2) * Sigma %*% (Phi[I,] %*% y)
  theta[I] = Mu
  return(theta)
}

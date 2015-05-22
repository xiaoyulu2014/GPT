
tensor_inf = function(X,y,sigma,n,r,sigmaRBF,q,generator,num_iterations,burnin){
  ##precompute features
  D = ncol(X)
  res = matrix(,num_iterations+1,n^D);alpha = matrix(,num_iterations,D)
  Phi_vec = apply(X,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,D)),2)) })
  Phi_vec_array = array(,c(n^(D-1),nrow(X),D))   ###Phi_{-j}
  for (j in 1:D) {
    Phi_vec_array[,,j] = apply(X,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator)[,-j],c(n,D-1)),2)) })
  }

  b = array(,c(n,nrow(X),D))
  for (j in 1:D){
    for (i in 1:nrow(X)){
      b[,i,j] = feature(X[i,],n,sigmaRBF,generator)[,j]
    }
  }

  #initialise U's
  U_array = array(,c(n,r,ncol(X)))
  for (j in 1:ncol(X)){U_array[,,j] = rustiefel(n,r)}
  
  ##initialise W
  ##randomly select the set of non-zero entries in W
  I = sort(sample(1:(r^D))[1:(q*r^D)])
#   wprior = diag(rep(n^D/length(I),length(I)))
  w = rep(0,r^D)
  w[I] = n^D/length(I)*rnorm(length(I))
#   w[I] = rmvnorm(1,rep(0,length(I)),wprior)
  wprior_inv = diag(rep(length(I)/(n^D),length(I)))  

  
  U = Reduce(kronecker,alply(U_array,3)) 
  U_I = U[,I]
  res[1,] = U_I %*% w[I]


  ###inference
  for (m in 1:num_iterations){
    tmp0 = apply(X,1,function(x){
      Reduce(kronecker,lapply(1:D,function(z){t(U_array[,,z]) %*% feature(x,n,sigmaRBF,generator)[,z]}))[I] } )   
    Sigma = chol2inv(chol(wprior_inv + (tmp0 %*% (t(tmp0)/(sigma^2))) ))
    w[I] = Sigma %*% (t(U_I) %*% (Phi_vec %*% y)) /(sigma^2)
    W = array(w,rep(r,D))   
    
    ##posterior on the Us 
    for (j in 1:D){
      
      Wj = array(aperm(W,c(j,(1:D)[-j])),c(r,r^(D-1)))
      #form U_(minus j)
      if (D==2) {Uj = U_array[,,(1:D)[-j]]}
      else {Uj = Reduce(kronecker,alply(array(U_array[,,-j],c(n,r,D-1)),3))}
      

      a = matrix(,r,nrow(X));
      for (i in 1:nrow(X)){ a[,i] = Wj %*% (t(Uj) %*% Phi_vec_array[,i,j]) }
     
      
      log_pi = function(u) {  
        Reduce("+", lapply(c(1:nrow(X)),function(x){ tmp1 = t(b[,x,j]) %*% u %*% a[,x]
                                               return(-1/(2*sigma^2) * tmp1^2 + y[x] * tmp1 /(sigma^2)) }) ) }
      
      log_dpi = function(u) {
        Reduce("+", lapply(c(1:nrow(X)),function(x){ tmp1 = as.numeric(t(b[,x,j]) %*% u %*% a[,x])
                                               return(1/(sigma^2) * (y[x] - tmp1) * (b[,x,j] %*% t(a[,x])) ) }) ) }
      
      
      tmp2 = HMC(U_array[,,j],log_pi,log_dpi,T=20,epsilon=0.002,num_iteration=10)
      U_array[,,j] = tmp2$x; alpha[m,j] = tmp2$alpha
      } 

      U_I = Reduce(kronecker,alply(U_array,3))[,I]
      tmp3 = U_I  %*%  w[I]
     if ( m > burnin) {res[m+1,] = (tmp3 + res[m,])/2 }
     else { res[m+1,] = tmp3 }  
    }
  return(list(res = res, alpha = alpha))
}




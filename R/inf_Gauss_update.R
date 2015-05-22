# X = X_train; y = y_train ; 
# n=100; r=30 ; num_iterations=50 ;burnin=40; generator=123 ;q=0.1


Gauss_inf = function(X,y,sigma,n,r,sigmaRBF,q,generator,num_iterations,burnin){
  ##precompute features
  D = ncol(X);
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
  sigma_u = sqrt(1/r)
  U_array = array(,c(n,r,ncol(X)))
  for (j in 1:ncol(X)){U_array[,,j] = matrix(sigma_u*rnorm(n*r),n,r)}
  ##initialise W
  w = rep(0,r^D)
  I = sort(sample(1:(r^D))[1:(q*r^D)])
  w[I] = sqrt(1/q) * rnorm(length(I))
  res[1,] = Reduce(kronecker,alply(U_array[,,rev(1:D)],3))[,I] %*% w[I]
  
  U = Reduce(kronecker,alply(U_array[,,rev(1:D)],3)) 
  U_I = U[,I]
  
  ###inference
  for (m in 1:num_iterations){
    ##posterior for the core tensor W    
    tmp0 = apply(X,1,function(x){
      Reduce(kronecker,lapply(rev(1:D),function(z){t(U_array[,,z]) %*% feature(x,n,sigmaRBF,generator)[,z]}))[I] } )   
    Sigma = chol2inv(chol( q * diag(q*r^D) + (tmp0 %*% (t(tmp0)/(sigma^2))) ))
    w[I] = Sigma %*% (t(U_I) %*% (Phi_vec %*% y)) /(sigma^2)
#     W = aperm(array(w,rep(r,D)),rev(1:D))
    W = array(w,rep(r,D))
    
    
    
    
    
    ##posterior on the Us 
    for (j in 1:D){
      
      Wj = array(aperm(W,c(j,(1:D)[-j])),c(r,r^(D-1)))
      #       Wj = aperm(W,c((1:D)[-j],j))
      
      #form U_(minus j)
      if (D==2) {Uj = U_array[,,-j]}
      else {Uj = Reduce(kronecker,alply(array(U_array[,,-j],c(n,r,D-1)),3))}
      
      a = matrix(,r,nrow(X));C = matrix(,r*n,nrow(X));
      for (i in 1:nrow(X)){ a[,i] = Wj %*% (t(Uj) %*% Phi_vec_array[,i,j]) 
                            C[,i] = kronecker(a[,i],b[,i,j])
      }
      log_pi = function(u) {  
        Reduce("+", lapply(c(1:nrow(X)),function(x){ tmp1 = t(b[,x,j]) %*% u %*% a[,x]
                                                     return(-1/(2*sigma^2) * tmp1^2 + y[x] * tmp1 /(sigma^2))  } )) - 1/(2*sigma_u^2) * sum(u^2) }  
      
      print(paste("before = ", log_pi(U_array[,,j])))
      
      U_array[,,j] = matrix(solve(C %*% t(C) + (sigma/sigma_u)^2 * diag(n*r)) %*% (C %*% y),n,r)
      print(paste("after = ", log_pi(U_array[,,j])))
      
      #       logl[m] = log_pi(U_array[,,j])
    } 
    
    U_I = Reduce(kronecker,alply(U_array[,,rev(1:D)],3))[,I]
    TMP = U_I %*% w[I]
    if ( m > burnin) {res[m+1,] = (TMP + res[m,])/2 }
    else { res[m+1,] = TMP }  
    
    #     f_hat = t(Phi_vec) %*% res[m+1,]
    #     f_hat = unlist(lapply(c(1:nrow(X)),function(x){ t(b[,x,j]) %*% U_array[,,j] %*% a[,x] }))    
    #     print(sqrt(sum((f_hat-f[1:nrow(X_train)])^2)/nrow(X)))
    print(m)
  }
  return(res)
}
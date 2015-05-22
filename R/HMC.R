
HMC = function(x,log_pi,log_dpi,T,epsilon,num_iteration){
  
  n=nrow(x);r=ncol(x)
  alpha = rep(0,num_iteration)
  
  for (i in 1:num_iteration){
  
    V = matrix(rnorm(n*r),n,r)
    V = V - 0.5*x %*% (t(x) %*% V + t(V) %*% x)
    h = log_pi(x) - 0.5*sum(V^2)
    x_new = x

      for (t in 1:T){
  
          V = V + epsilon/2 * log_dpi(x_new)
          V = V - 0.5*x_new %*% (t(x_new) %*% V + t(V) %*% x_new)
          
          tmp1 = t(x_new) %*% V; tmp2 = as.matrix(expm(-epsilon*tmp1))
          tmp3 = cbind(x_new,V) %*% as.matrix(expm(epsilon* (rbind(cbind(tmp1, -t(V) %*% V),cbind(diag(r),tmp1))))) %*% as.matrix(bdiag(tmp2,tmp2))
          x_new = tmp3[,1:r]; V = tmp3[,(r+1):(2*r)]
          
          V = V + epsilon/2 * log_dpi(x_new)
          V = V - 0.5*x_new %*% (t(x_new) %*% V + t(V) %*% x_new)    
      }
    
    h_new = log_pi(x_new) - 0.5*sum(V^2)
    if (is.na(h_new)==F && runif(1) < exp(h_new-h))  {x = x_new; alpha[i] = 1}    
  }
  
  return(list(x=x,alpha=sum(alpha)/num_iteration))
  
}

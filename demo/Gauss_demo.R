
N = 30;sigma=0.1;sigmaRBF = 5;
seq_train = seq(0.5,by=1,length.out=N);seq_test = sort(N*runif(N/4))
X_train = cbind(rep(seq_train,N),rep(seq_train,each=N))
X_test = cbind(rep(seq_test,N/4),rep(seq_test,each=N/4))
f = GPdraw(rbind(X_train,X_test),sigmaRBF,1)
y = f + sigma*rnorm(length(f))
y_train = y[1:nrow(X_train)]
y_test = tail(y,nrow(X_test))
f_test = tail(f,nrow(X_test))
f_train = f[1:nrow(X_train)]
##################################################################################################################

# n=50;generator=123;q=0.1;num_iterations=30;burnin=25;r=10
# out1 = Gauss_inf(X_train,y_train,sigma,n,r,sigmaRBF,q=1,generator,num_iterations,burnin)

#RMSE
Phi_train = apply(X_train,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_train))),2)) })
Phi_test = apply(X_test,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_test))),2)) })
RMSE_train=c();RMSE_test=c()
for (i in 1:num_iterations){
  RMSE_train[i] = sqrt(sum((t(Phi_train) %*% res[i,]-y_train)^2)/nrow(X_train))
  RMSE_test[i] = sqrt(sum((t(Phi_test) %*% res[i,]-y_test)^2)/nrow(X_test))
}

#RMSE_array
generator=123;num_iterations=20;burnin=10;
RMSE_train_array = array[,c(n_vec,q_vec)]; RMSE_test_array = array[,c(n_vec,q_vec)]
n_vec = c(50,100,150);q_vec = c(100,500,1000);r=40
for (i in 1:length(n_vec)) {
  n = n_vec[i]
  Phi_train = apply(X_train,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_train))),2)) })
  Phi_test = apply(X_test,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_test))),2)) })
  for (j in 1:length(q_vec)) {
    q = q_vec[j]
    res = Gauss_inf(X_train,y_train,sigma,n,r,sigmaRBF,q,generator,num_iterations,burnin)
    RMSE_train_array[i,j] = sqrt(sum((t(Phi_train) %*% res[num_iterations,]-y_train)^2)/nrow(X_train))
    RMSE_test_array[i,j] = sqrt(sum((t(Phi_test) %*% res[num_iterations,]-y_test)^2)/nrow(X_test))   
  }
}


##########################################################################################################


######################################GP EXACT#####################################################
Kernel = function(x,y,sigmaRBF) {
  tmp = matrix(1,nrow(x),nrow(y))
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(y)){
      tmp[i,j] = exp(-1/(2*sigmaRBF^2)*(t(x[i,]-y[j,])%*%(x[i,]-y[j,])))
    }
  } 
  return(tmp)
}

tmp = solve(Kernel(X_train,X_train,sigmaRBF) + sigma^2 * diag(nrow(X_train))) %*% y_train
mean_test = Kernel(X_test,X_train,sigmaRBF) %*% tmp
mean_train = Kernel(X_train,X_train,sigmaRBF) %*% tmp

plot(y_test,type="l")
lines(mean_test,col=3)
RMSE_train_ex = sqrt(sum((mean_test-y_test)^2)/nrow(X_test))
RMSE_test_ex = sqrt(sum((mean_train-y_train)^2)/nrow(X_train))



# persp3D(seq_train,seq_train,matrix(f,N,N),theta = 50, phi = 30, 
#         expand = 0.5,xlab="x",ylab="y",zlab="z",main="GP draw")
# 
# persp3D(seq_test,seq_test,matrix(f_test,N/2,N/2),theta = 50, phi = 30, 
#         expand = 0.5,xlab="x",ylab="y",zlab="z",main="GP draw")
# 
# 
# 
# 








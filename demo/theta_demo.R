
N = 30;sigma=0.5;sigmaRBF = 5;
seq_train = seq(0.5,by=1,length.out=N);seq_test = sort(N*runif(N/4))
X_train = cbind(rep(seq_train,N),rep(seq_train,each=N))
X_test = cbind(rep(seq_test,N/4),rep(seq_test,each=N/4))
f = GPdraw(rbind(X_train,X_test),sigmaRBF,1)
y = f + sigma*rnorm(length(f))
y_train = y[1:nrow(X_train)]
f_test = tail(f,nrow(X_test))
##################################################################################################################

generator=123
n_vec = c(50,100,150)
q_vec = c(1/3,1/2,1)

  
library(parallel)
numWorkers <- 8

RMSE_train = array(,c(length(n_vec),length(q_vec)));RMSE_test = array(,c(length(n_vec),length(q_vec)))
for (i in 3:length(n_vec)) {
  n = n_vec[i]
  Phi_train = apply(X_train,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_train))),2)) })
  Phi_test = apply(X_test,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_test))),2)) })
  
  workerFunc = function(q) {theta_inf(X_train,y_train,sigma,n,sigmaRBF,q,generator)}
  cl <- makeCluster(numWorkers, type = "PSOCK")
  clusterExport(cl,"theta_inf")
  clusterExport(cl=cl, varlist=c("X_train", "y_train", "sigma", "n", "sigmaRBF","generator"), envir=environment())
  res <- Reduce(cbind,parLapply(cl, q_vec, workerFunc))
  stopCluster(cl)

  RMSE_train[i,] =   sqrt(colSums((t(Phi_train) %*% res - f[1:nrow(X_train)])^2)/nrow(X_train))
  RMSE_test[i,] =   sqrt(colSums((t(Phi_test) %*% res - f_test)^2)/nrow(X_test))
  print(paste("i= ",i))
}


plot(q_vec,RMSE_test[1,],type="n",ylim=c(0.1,0.21),main="RMSE on test set")
for (i in 1:3) {lines(q_vec,RMSE_test[i,],col=i)}; abline(h = RMSE1,col=4)
legend("topright",paste("n= ",n_vec),col=1:3,lwd=rep(2,3));text(0.4,0.115,paste("RMSE exact = ", round(RMSE1,2)))

plot(q_vec,RMSE_train[1,],type="n",ylim=c(0,0.21),main="RMSE on training set")
for (i in 1:3) {lines(q_vec,RMSE_train[i,],col=i)}; abline(h = RMSE_train,col=4)
legend("topright",paste("n= ",n_vec),col=1:3,lwd=rep(2,3));text(0.4,0.115,paste("RMSE exact = ", round(RMSE1,2)))


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

# X_test = X_train
pos_ker =  Kernel(X_test,X_train,sigmaRBF)
pos_tmp =  pos_ker %*% solve(Kernel(X_train,X_train,sigmaRBF) + sigma^2 * diag(nrow(X_train)))
pos_mean = pos_tmp %*% y_train
pos_cov =  Kernel(X_test,X_test,sigmaRBF) -pos_tmp %*% t(pos_ker)

plot(f_test,type="l")
lines(pos_mean,col=3)
RMSE_test = sqrt(sum((pos_mean-f_test)^2)/nrow(X_test))
RMSE_train = sqrt(sum((pos_mean-f_test)^2)/nrow(X_test))


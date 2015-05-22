
N = 30;sigma=0.5;sigmaRBF = 5;
seq_train = seq(0.5,by=1,length.out=N);seq_test = sort(N*runif(N/4))
X_train = cbind(rep(seq_train,N),rep(seq_train,each=N))
X_test = cbind(rep(seq_test,N/4),rep(seq_test,each=N/4))
f = GPdraw(rbind(X_train,X_test),sigmaRBF,1)
y = f + sigma*rnorm(length(f))
y_train = y[1:nrow(X_train)]
f_test = tail(f,nrow(X_test))
##################################################################################################################
persp3D(seq_train,seq_train,matrix(f,N,N),theta = 50, phi = 30, 
        expand = 0.5,xlab="x",ylab="y",zlab="z",main="GP draw")

persp3D(seq_test,seq_test,matrix(f_test,N/2,N/2),theta = 50, phi = 30, 
        expand = 0.5,xlab="x",ylab="y",zlab="z",main="GP draw")

n=100;generator=123;q=0.1;num_iterations=30;burnin=25;r=30
out1 = Gauss_inf(X_train,y_train,sigma,n,r,sigmaRBF,q=0.1,generator,num_iterations,burnin)



#training RMSE
Phi_vec = apply(X_train,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_train))),2)) })
plot(f[1:nrow(X_train)],type="l",main=paste("n=",n, ", r= ", r),ylim=c(-4,3));
RMSE_train2=c()
for (i in 1:22){
  f_fit = t(Phi_vec) %*% res[i,]
  #   lines(f_fit,col=i+1)
  RMSE_train2[i] = sqrt(sum((f_fit-f[1:nrow(X_train)])^2)/nrow(X_train))
}
#test RMSE
Phi_vec = apply(X_test,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_test))),2)) })
plot(f_test,type="l",main=paste("n=",n, ", r= ", r),ylim=c(-4,3));
RMSE_test2=c()
for (i in 1:22){
  f_fit = t(Phi_vec) %*% res[i,]
  lines(f_fit,col=i+1)
  RMSE_test2[i] = sqrt(sum((f_fit-f_test)^2)/nrow(X_test))
}

par(mfrow=c(1,2))
plot(2:22,RMSE_train1[2:22],type="l",xlab="number of iterations",main="RMSE_train for r=30, q=2/3",ylim=c(0.11,0.16));
lines(2:22,RMSE_train2[2:22],col=3)
abline(h=RMSE2,col=2);text(8,0.12,paste("RMSE exact = ",round(RMSE2,3)),col=2)
legend("topright",c("n=100","n=50"),col=c(1,3),lwd=c(2,2))
plot(2:22,RMSE_test1[2:22],type="l",xlab="number of iterations",main="RMSE_test for r=30, q=2/3",ylim=c(0.09,0.15));
lines(2:22,RMSE_test2[2:22],col=3)
abline(h=RMSE1,col=2);text(8,0.096,paste("RMSE exact = ",round(RMSE1,3)),col=2)
legend("topright",c("n=100","n=50"),col=c(1,3),lwd=c(2,2))



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


pos_ker =  Kernel(X_test,X_train,sigmaRBF)
pos_tmp =  pos_ker %*% solve(Kernel(X_train,X_train,sigmaRBF) + sigma^2 * diag(nrow(X_train)))
pos_mean = pos_tmp %*% y_train
pos_cov =  Kernel(X_test,X_test,sigmaRBF) -pos_tmp %*% t(pos_ker)

plot(f_test,type="l")
lines(pos_mean,col=3)
RMSE1 = sqrt(sum((pos_mean-f_test)^2)/nrow(X_test))
RMSE2 = sqrt(sum((pos_mean-f[1:nrow(X_train)])^2)/nrow(X_train))















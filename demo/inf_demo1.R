
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
# out3 = tensor_inf(X_train,y_train,sigma,n,r=50,sigmaRBF,q=0.01,generator,num_iterations=10,burnin)
out1 = Gauss_inf(X_train,y_train,sigma,n,r,sigmaRBF,q=0.1,generator,num_iterations,burnin)


n_vec = c(50,100,150)
r_vec = c(1/3,1/2,1)
q_vec = c(1/3,1/2,1)
generator=123;num_iterations=8;burnin=6;

library(parallel)
library(doParallel)
registerDoParallel(cores=8)
numWorkers <- 8

RMSE_train = array(,c(length(n_vec),length(r_vec),length(q_vec)));RMSE_test = array(,c(length(n_vec),length(r_vec),length(q_vec)))
for (i in 1:3){
  n = n_vec[i]
  Phi_train = apply(X_train,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_train))),2)) })
  Phi_test = apply(X_test,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_test))),2)) })
 
  out <- foreach(j=1:length(r_vec), .combine=cbind) %dopar% {
    r = round(n*r_vec[j])
    
    workerFunc = function(q) {Gauss_inf(X_train,y_train,sigma,n,r,sigmaRBF,q,generator,num_iterations,burnin)[num_iterations,]}
    cl <- makeCluster(numWorkers, type = "PSOCK")
    clusterExport(cl,list("Gauss_inf","alply","feature"))
    clusterExport(cl=cl, varlist=c("X_train", "y_train", "sigma", "n","r", "sigmaRBF","generator","num_iterations","burnin"), envir=environment())
    res <- Reduce(cbind,parLapply(cl, q_vec, workerFunc))
    
    append(sqrt(colSums((t(Phi_train) %*% res - f[1:nrow(X_train)])^2)/nrow(X_train)),
    sqrt(colSums((t(Phi_test) %*% res - f_test)^2)/nrow(X_test)) )
  }
}
stopCluster(cl)

par(mfrow=c(1,2))
plot(q_vec,RMSE_test[1,1,],type="n",ylim=c(0.14,0.17),main="RMSE on test set when n=50",xlab="q")
for (i in 1:3) {lines(q_vec,RMSE_test[1,i,],col=i)}; abline(h=round(RMSE1,3),col=4)
legend("topright",c("r=1/3","r=2/3","r=1"), lwd=c(2,2,2),col=1:3); text(0.5,0.142,paste("RMSE exact = ",round(RMSE1,3)))

plot(q_vec,RMSE_train[1,1,],type="n",ylim=c(0.12,0.17),main="RMSE on training set when n=50",xlab="q")
for (i in 1:3) {lines(q_vec,RMSE_train[1,i,],col=i)}; abline(h=round(RMSE2,3),col=4)
legend("topright",c("r=1/3","r=2/3","r=1"), lwd=c(2,2,2),col=1:3);text(0.5,0.125,paste("RMSE exact = ",round(RMSE2,3)))


# library(parallel)
# numWorkers <- 8
# 
# RMSE_train = array(,c(length(n_vec),length(r_vec),length(q_vec)));RMSE_test = array(,c(length(n_vec),length(r_vec),length(q_vec)))
# for (i in 1:3){
#   n = n_vec[i]
#   Phi_train = apply(X_train,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_train))),2)) })
#   Phi_test = apply(X_test,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_test))),2)) })
#   for (j in 1:length(r_vec)){
#     r = round(n*r_vec[j])
#     
#     workerFunc = function(q) {Gauss_inf(X_train,y_train,sigma,n,r,sigmaRBF,q,generator,num_iterations,burnin)[num_iterations,]}
#     cl <- makeCluster(numWorkers, type = "PSOCK")
#     clusterExport(cl,list("Gauss_inf","alply","feature"))
#     clusterExport(cl=cl, varlist=c("X_train", "y_train", "sigma", "n","r", "sigmaRBF","generator","num_iterations","burnin"), envir=environment())
#     res <- Reduce(cbind,parLapply(cl, q_vec, workerFunc))
#     stopCluster(cl)
#     
#     RMSE_train[i,j,] = sqrt(colSums((t(Phi_train) %*% res - f[1:nrow(X_train)])^2)/nrow(X_train))
#     RMSE_test[i,j,] =   sqrt(colSums((t(Phi_test) %*% res - f_test)^2)/nrow(X_test))
#     print(paste("j= ", j))
#   }
#   print(paste("i= ", i))
# }



RMSE_train = array(,c(length(n_vec),length(r_vec),length(q_vec)));RMSE_test = array(,c(length(n_vec),length(r_vec),length(q_vec)))
for (i in 1:3){
  n = n_vec[i]
  Phi_train = apply(X_train,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_train))),2)) })
  Phi_test = apply(X_test,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X_test))),2)) })
  for (j in 1:length(r_vec)){
    r = round(n*r_vec[j])
    res = list()
    for (k in 1:3){
      q = q_vec[k]
      res[[k]] = Gauss_inf(X_train,y_train,sigma,n,r,sigmaRBF,q,generator,num_iterations,burnin)
      RMSE_train[i,j,k] = sqrt(colSums((t(Phi_train) %*% res[[k]][num_iterations,] - f[1:nrow(X_train)])^2)/nrow(X_train))
      RMSE_test[i,j,k] =   sqrt(colSums((t(Phi_test) %*% res[[k]][num_iterations,] - f_test)^2)/nrow(X_test))
      print(paste("k= ", k))
    }
    print(paste("j= ", j))
  }
  print(paste("i= ", i))
}
}



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





q_vec = c(0.01,0.05,0.1,0.2,0.5,1);RMSE_150 = c()
for (i in 1:length(q_vec)){
  q = q_vec[i]
  out5 = theta_inf(X_train,y_train,sigma,n,sigmaRBF,q,generator)
  f_fit = t(Phi_vec) %*% out5
  RMSE_150[i] = sqrt(sum((f_fit-f_test)^2)/nrow(X_test))
  print(i)
}
plot(q_vec,RMSE_50,type="l")



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




# kernelRBF = function(x,y)  exp(-1/(2*sigmaRBF^2)*(t(x-y)%*%(x-y)))
# nr = nrow(X)
# cov = matrix(0,nr,nr)
# for (i in 1:nr){
#   for (j in i:nr){
#     cov[i,j] = kernelRBF(X[i,],X[j,])
#   }
# }
# K = cov + t(cov) - diag(nr)
# 













n_vec = c(50,80,100,150)
r_vec = c(1/4,1/2,3/4,1)
q_vec = c(0.1,0.4,0.7,1)

f_array = array(,c(N^ncol(X),length(n_vec),length(r_vec),length(q_vec)));RMSE_array = array(,c(length(n_vec),length(r_vec),length(q_vec)))
for (i in 3:3){
  n = n_vec[i]
  Phi_vec = apply(X_test,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X))),2)) })
  for (j in 1:length(r_vec)){
    r = round(n*r_vec[j])
      for (k in 1:3){
        q = q_vec[k]
        theta = tensor_inf(X,y,sigma,n,r,sigmaRBF,q,generator,num_iterations,burnin)[num_iterations,]
        f_array[,i,j,k] = t(Phi_vec) %*% theta
        RMSE_array[i,j,k] = sqrt(sum((f_array[,i,j,k]-fX)^2)/nrow(X))
      }
    }
  }
}






par(mfrow=c(1,2))
par(mar=rep(5,4))

Phi_vec = apply(X,1,function(x){Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X))),2)) })
persp3D(sequence,sequence,matrix(fX,N,N),theta = 50, phi = 30, 
        expand = 0.5,xlab="x",ylab="y",zlab="z",main="GP draw")

plot(fX,type="l",main=paste("n=",n, ", r= ", r),ylim=c(-4,3));
RMSE=c()
for (i in 1:num_iterations){
   f_fit = t(Phi_vec) %*% f_array[i,]
   persp3D(sequence,sequence,matrix(f_array[,i,j,k],N,N),theta = 50, phi = 30,
           expand = 0.5,xlab="x",ylab="y",zlab="z",main=paste("fitted GP when n=",n,", r=",r_vec[i]," ,q=",q))
   lines(f_fit,col=i+1)
   RMSE[i] = sqrt(sum((f_fit-fX)^2)/nrow(X))
}

persp3D(sequence,sequence,matrix(pos_mean,N,N),theta = 50, phi = 30,
        expand = 0.5,xlab="x",ylab="y",zlab="z",main="fitted GP using traditional method")


plot(r_vec[1:4],RMSE[1:4],type="l",main="RMSE vs r when n=100, q=0.1",xlab="r")

##fix q
par(mfrow=c(2,2))
for (k in 1:length(q_vec)){
  plot(r_vec,RMSE_array[1,,k],type="n",ylim=c(0.1,max(RMSE_array[1,,k])),xlab="r as a percentage of n",ylab="RMSE",
       main=paste("RMSE on a 50*50 grid for q= ",q_vec[k]))
  for (i in 1:4) {lines(r_vec,RMSE_array[i,,k],col=i)}; 
  abline(h=RMSE1,col=6); text(0.35,0.1,paste("RMSE_exact= ", round(RMSE1,2)))
  legend("topright",c(paste("n= ",n_vec)),col=1:4,lwd=rep(1,4))
}

#fix r
par(mfrow=c(2,2))
for (j in 1:length(r_vec)){
  plot(q_vec,RMSE_array[1,j,],type="n",ylim=c(0.1,max(RMSE_array[1,j,])),xlab="q",ylab="RMSE",
       main=paste("RMSE on a 50*50 grid for r= ",r_vec[j]))
  for (i in 1:4) {lines(q_vec,RMSE_array[i,j,],col=i)}; 
  abline(h=RMSE1,col=6); text(0.35,0.1,paste("RMSE_exact= ", round(RMSE1,2)))
  legend("topright",c(paste("n= ",n_vec)),col=1:4,lwd=rep(1,4))
}

#fix n
par(mfrow=c(2,2))
for (i in 1:length(n_vec)){
  plot(q_vec,RMSE_array[i,1,],type="n",ylim=c(0.1,max(RMSE_array[i,1,])),xlab="q",ylab="RMSE",
       main=paste("RMSE on a 50*50 grid for n= ",n_vec[i]))
  for (j in 1:4) {lines(r_vec,RMSE_array[i,j,],col=i)}; 
  abline(h=RMSE1,col=6); text(0.35,0.1,paste("RMSE_exact= ", round(RMSE1,2)))
  legend("topright",c(paste("r= ",r_vec)),col=1:4,lwd=rep(1,4))
}
##surface plot
saveGIF({
  ani.options(interval = 0.5)

  par(bg = "white") 
  split.screen(c(1,2)) # split display into two screens
  screen(1) # prepare screen 1 for output
  persp3D(sequence,sequence,matrix(fX,N,N),theta = 50, phi = 30, 
          expand = 0.5,xlab="x",ylab="y",zlab="z",main="GP draw")
  
  screen(2) # prepare screen 1 for output
  for (i in 1:5) {
    f_fit = t(Phi_vec) %*% Theta_varyr[,,1][i,]
    persp3D(sequence,sequence,matrix(f_fit,N,N),theta = 50, phi = 30,
          expand = 0.5,xlab="x",ylab="y",zlab="z",main=paste("fitted GP when n=",n,", r=",30," ,q=",q))
    ani.pause()
  }
}, img.name = "sur_plot",autoplay = FALSE, interval = 0.5, movie.name = "surface.gif")



# ,col=rgb(1, 1, 1, 0.5)
# 


# func = function(x) {Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X))),2)) %*% Theta}
# f_fit = apply(X,1,func)


# 
# ##profiling
# Rprof("prof.out")
# x = profr(tensor_inf(X,y,sigma,n=10,r=10,sigmaRBF,q=0.5,generator=1,num_iterations=50,burnin=10))
# Rprof(NULL)
# ggplot.profr(parse_rprof("prof.out"))

# par(mfrow=c(2,3))
# nvec= c(10,30,50,70,100)
# for (i in 1:5) {
#   n = nvec[i];r=n
#   Theta = tensor_inf(X,y,sigma,n,r,sigmaRBF,q=0.5,generator,num_iterations=50,burnin=10)
#   func = function(x) {Reduce(kronecker,alply(array(feature(x,n,sigmaRBF,generator),c(n,ncol(X))),2)) %*% Theta}
#   f_fit = apply(X,1,func)
#   plot(fX,type="l",main=paste("n=",n, ", r= ", r));lines(f_fit,col=2)
#   theta_list[[i]] = Theta
# }
##compare with traditional way
# kernelRBF = function(x,y)  exp(-1/(2*sigmaRBF^2)*(t(x-y)%*%(x-y)))
# nr = nrow(X)
# cov = matrix(0,nr,nr)
# for (i in 1:nr){
#   for (j in i:nr){
#     cov[i,j] = kernelRBF(X[i,],X[j,])
#   }
# }
# K = cov + t(cov) - diag(nr)
# 

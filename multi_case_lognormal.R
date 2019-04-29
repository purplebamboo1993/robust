setwd("C:/Users/weichenw/Dropbox/robust/code_2019_Weichen")
source("utl.R")
set.seed(111)

N = 1
f11 = matrix(0, N, 8)
f12 = matrix(0, N, 8)
f21 = matrix(0, N, 8)
f22 = matrix(0, N, 8)
f31 = matrix(0, N, 8)
f32 = matrix(0, N, 8)

f1_para = matrix(0, N, 2)
f2_para = matrix(0, N, 2)
f3_para = matrix(0, N, 2)

genTrueTheta = function(d){
  r0=5    
  n0=100
  Z=matrix(rnorm(n0*d,0,1), n0, d)    
  S=t(Z)%*%Z
  V=svd(S, nu=r0, nv=r0)$u[,1:r0]
  Theta_t=V%*%t(V)
  return(Theta_t)
}


CV = function(Y,X,lambda){
  fold = 5
  nn = floor(dim(Y)[1] / fold)
  CV_error = rep(0, fold)
  for (h in 1:fold) {
    if (h < fold)  testSam = 1:nn + (h-1)*nn
    else  testSam = ((fold-1)*nn+1):(dim(Y)[1])
    Y_train = Y[-testSam,]
    X_train = X[-testSam,]
    lambda_train = lambda * sqrt(dim(Y)[1]) / sqrt(dim(Y_train)[1])
    
    Theta = prsm_multi(Y_train, X_train, alpha=.9, beta=1, lambda=lambda_train)
    Y_target = matrix(0, length(testSam), d)
    for (dd in 1:d) {
      Y_th = quantile(abs(Y[testSam,dd]),0.95)
      Y_target[,dd] = truncate_abs(Y[testSam,dd], Y_th)
    }
    CV_error[h] = sum((Y_target - X[testSam,] %*% Theta)^2) / length(testSam)
  }
  return(c(mean(CV_error), sd(CV_error) / sqrt(fold)))
}


CV_thresh = function(Y,X,th,lambda1){
  fold = 5
  nn = floor(dim(Y)[1] / fold)
  CV_error = rep(0, fold)
  for (h in 1:fold) {
    if (h < fold)  testSam = 1:nn + (h-1)*nn
    else  testSam = ((fold-1)*nn+1):(dim(Y)[1])
    Y_train = Y[-testSam,]
    X_train = X[-testSam,]
    lambda_train = lambda1 * sqrt(dim(Y)[1]) / sqrt(dim(Y_train)[1])
    th_train = th/sqrt(dim(Y)[1])*sqrt(dim(Y_train)[1])
    
    Theta = prsm_multi(truncate_l2(Y_train, tau=th_train), X_train, alpha=.9, beta=1, lambda=lambda_train)
    Y_target = matrix(0, length(testSam), d)
    for (dd in 1:d) {
      Y_th = quantile(abs(Y[testSam,dd]),0.95)
      Y_target[,dd] = truncate_abs(Y[testSam,dd], Y_th)
    }
    CV_error[h] = sum((Y_target - X[testSam,] %*% Theta)^2) / length(testSam)
  }
  return(c(mean(CV_error), sd(CV_error) / sqrt(fold)))
}


#for (d in c(50)) {
for (d in c(50, 100, 150)) {
  cat("### d =", d, "### \n")
  
  f1 = matrix(0, N, 8)
  f2 = matrix(0, N, 8)
  para = matrix(0, N, 2)
  
  for (j in 1:N){
    n=500
    r=5
    Theta_t = genTrueTheta(d)
    
    #sigma = .5
    #epsilon=matrix(rnorm(n*d, 0, sigma), n, d)
    #epsilon=randomt(n, 1)/128
    sigma=2.8
    epsilon=matrix((exp(rnorm(n*d, 0, sigma))-exp(sigma^2/2)*rep(1,n*d)), n, d)/500
    X=matrix(rnorm(n*d, 0, 1), n, d)
    Y=X%*%Theta_t+epsilon

    lambda_cv_list = c(0:10) * 8 + 8
    f11_cv = rep(0, length(lambda_cv_list))
    f11_cv_se = rep(0, length(lambda_cv_list))
    
    for (k in 1:length(lambda_cv_list)){
      lambdac = lambda_cv_list[k]
      
      lambda = lambdac * sqrt(log(d) / (n*d))
      obj = CV(Y,X,lambda)
      f11_cv[k] = obj[1]
      f11_cv_se[k] = obj[2]
      
      cat("======= lambdac CV",  k, "=======", "\n")
    }
    
    k_min = which(f11_cv == min(f11_cv))
    k_min = k_min[length(k_min)]
    kvec = which(f11_cv < f11_cv[k_min] + f11_cv_se[k_min])
    k_min_1sd = kvec[length(kvec)]
    lambdac1 = lambda_cv_list[k_min_1sd]
    lambdac2 = lambdac1
    para[j,1] = lambdac1
    
    th_cv_list = c(0:10) * 2 + 2
    f12_cv = rep(0, length(th_cv_list))
    f12_cv_se = rep(0, length(th_cv_list))
    
    for (k in 1:length(th_cv_list)){
      thc = th_cv_list[k]
      
      lambda1=lambdac2 * sqrt(log(d) / (n*d))
      th=thc*sqrt(n/(d*log(d)))
      obj = CV_thresh(Y,X,th,lambda1)
      f12_cv[k] = obj[1]
      f12_cv_se[k] = obj[2]
      
      cat("======= thc CV",  k, "=======", "\n")
    } 
    
    k_min = which(f12_cv == min(f12_cv))
    k_min = k_min[1]
    kvec = which(f12_cv < f12_cv[k_min] + f12_cv_se[k_min])
    k_min_1sd = kvec[1]
    thc = th_cv_list[k_min_1sd]
    para[j,2] = thc
    
    
    for (i in 1:8){
      n=500*i
      #sigma = .5
      #epsilon=matrix(rnorm(n*d, 0, sigma), n, d)
      #epsilon=randomt(n, 1)/128
      sigma=2.8
      epsilon=matrix((exp(rnorm(n*d, 0, sigma))-exp(sigma^2/2)*rep(1,n*d)), n, d)/500
      X=matrix(rnorm(n*d, 0, 1), n, d)
      Y=X%*%Theta_t+epsilon
      
      lambda1=lambdac1 * sqrt(log(d) / (n*d))
      Theta1=prsm_multi(Y, X, alpha=.9, beta=1, lambda=lambda1)
      f1[j,i]=sqrt(sum((Theta1-Theta_t)^2))
      
      lambda2=lambdac2 * sqrt(log(d) / (n*d))
      th=thc*sqrt(n/(d*log(d)))
      Theta2=prsm_multi(truncate_l2(Y, tau=th), X, alpha=.9, beta=1, lambda=lambda2)        
      f2[j,i]=sqrt(sum((Theta2-Theta_t)^2))
      
      cat("=======", j, ",", i, "=======", "\n")
    }
  }
  
  print(f1)
  print(f2)
  print(para)
  
  if (d == 50) {
    f11 = f1
    f12 = f2
    f1_para = para
  }
  if (d == 100) {
    f21 = f1
    f22 = f2
    f2_para = para
  }
  if (d == 150) {
    f31 = f1
    f32 = f2
    f3_para = para
  }
}
  
save(f11, f12, f21, f22, f31, f32, f1_para, f2_para, f3_para, file="multi_CV.RData")


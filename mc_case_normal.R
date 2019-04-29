setwd("/mhome/stats/t/zz385/robust")
source("utl.R")
set.seed(111)


N = 1
f11 = matrix(0, N, 6)
f12 = matrix(0, N, 6)
f21 = matrix(0, N, 6)
f22 = matrix(0, N, 6)
f31 = matrix(0, N, 6)
f32 = matrix(0, N, 6)

f1_para = matrix(0, N, 2)
f2_para = matrix(0, N, 2)
f3_para = matrix(0, N, 2)


genTrueTheta = function(d){
    repeat{
      r0=5    
      n0=100
      Z=matrix(rnorm(n0*d,0,1), n0, d)    
      S=t(Z)%*%Z
      V=svd(S, nu=r0, nv=r0)$u[,1:r0]
      Theta_t=as.vector(V%*%t(V))/sqrt(r0)
      max(abs(Theta_t))*d
      if (max(abs(Theta_t))*d<=6) {
        break    
      }
    }
    return(Theta_t)
}

CV = function(Y,X,lambda){
    fold = 5
    nn = floor(length(Y) / fold)
    CV_error = rep(0, fold)
    #CV_error_75 = rep(0, fold)
    for (h in 1:fold) {
        if (h < fold)  testSam = 1:nn + (h-1)*nn
        else  testSam = ((fold-1)*nn+1):(length(Y))
        Y_train = Y[-testSam]
        X_train = X[-testSam,]
        lambda_train = lambda * sqrt(length(Y)) / sqrt(length(Y_train))
      
        Theta=admm(Y=Y_train, X=X_train, rho=.1, mu=lambda_train, alpha=alpha)
        Y_th = quantile(abs(Y[testSam]),0.95)
        Y_target = truncate_abs(Y[testSam], Y_th) 
        for (ii in 1:length(testSam)){
            CV_error[h] = CV_error[h] + (Y_target[ii] - Theta[X[testSam[ii],1],X[testSam[ii],2]])^2
        }
        CV_error[h] = CV_error[h] / length(testSam)
        #CV_error_each = rep(0, length(testSam))
        #for (ii in 1:length(testSam)) {
        #    CV_error_each[ii] = abs(Y[testSam[ii]] - Theta[X[testSam[ii],1],X[testSam[ii],2]])
        #}
        #CV_error[h] = median(CV_error_each)
        #CV_error_75[h] = quantile(CV_error_each, 0.75)
    }
    return(c(mean(CV_error), sd(CV_error) / sqrt(fold)))
    #return(c(mean(CV_error), (mean(CV_error_75) - mean(CV_error))/sqrt(fold)))
}

CV_thresh = function(Y,X,th,lambda1){
    fold = 5
    nn = floor(length(Y) / fold)
    CV_error = rep(0, fold)
    #CV_error_75 = rep(0, fold)
    for (h in 1:fold) {
        if (h < fold)  testSam = 1:nn + (h-1)*nn
        else  testSam = ((fold-1)*nn+1):(length(Y))
        Y_train = Y[-testSam]
        X_train = X[-testSam,]
        lambda_train = lambda1 * sqrt(length(Y)) / sqrt(length(Y_train))
        th_train = th/sqrt(length(Y))*sqrt(length(Y_train))
        
        Theta=admm(Y=truncate_abs(Y_train, th_train), X=X_train, rho=.1, mu=lambda_train, alpha=alpha)
        Y_th = quantile(abs(Y[testSam]),0.95)
        Y_target = truncate_abs(Y[testSam], Y_th) 
        for (ii in 1:length(testSam)){
            CV_error[h] = CV_error[h] + (Y_target[ii] - Theta[X[testSam[ii],1],X[testSam[ii],2]])^2
        }
        CV_error[h] = CV_error[h] / length(testSam)
        #CV_error_each = rep(0, length(testSam))
        #for (ii in 1:length(testSam)) {
        #    CV_error_each[ii] = abs(Y[testSam[ii]] - Theta[X[testSam[ii],1],X[testSam[ii],2]])
        #}
        #CV_error[h] = median(CV_error_each)
        #CV_error_75[h] = quantile(CV_error_each, 0.75)
    }
    return(c(mean(CV_error), sd(CV_error) / sqrt(fold)))
    #return(c(mean(CV_error), (mean(CV_error_75) - mean(CV_error))/sqrt(fold)))
}



for (d in c(50, 100, 150)) {
  cat("### d =", d, "### \n")
  
  f1 = matrix(0, N, 6)
  f2 = matrix(0, N, 6)
  para = matrix(0, N, 2)
  
  for (j in 1:N){
      n=2000
      r=5
      Theta_t = genTrueTheta(d)
      alpha=6/d
      
      #epsilon=randomt(n, 1)/(64*d)
      #sigma = 3
      #epsilon=(exp(rnorm(n, 0, sigma))-exp(sigma^2/2)*rep(1,n))/(1000*d)
      sigma=.5
      epsilon=rnorm(n, mean=0, sd=sigma)/d
      obs=sample(x=1:(d^2), size=n, replace=T)
      X=matrix(0, n, 2)    
      X[, 1]=(obs-.1) %/% d+1
      X[, 2]=obs %% d
      X[X[, 2]==0, 2]=d
      Y=Theta_t[obs]+epsilon
  
      lambda_cv_list = c(0:10) * 0.1
      f11_cv = rep(0, length(lambda_cv_list))
      f11_cv_se = rep(0, length(lambda_cv_list))
      
      for (k in 1:length(lambda_cv_list)){
          lambdac = lambda_cv_list[k]
          
          lambda=lambdac*sqrt(d*log(d)/n)/(d^2)
          obj = CV(Y,X,lambda)
          f11_cv[k] = obj[1]
          f11_cv_se[k] = obj[2]
          
          cat("======= lambdac CV",  k, "=======", "\n")
      }
  
      k_min = which(f11_cv == min(f11_cv))
      kvec = which(f11_cv < f11_cv[k_min] + f11_cv_se[k_min])
      k_min_1sd = kvec[length(kvec)]
      lambdac1 = lambda_cv_list[k_min_1sd]
      #lambdac1 = lambda_cv_list[k_min]
      lambdac2 = lambdac1
      para[j,1] = lambdac1
  
  
      th_cv_list = c(0:10) * 0.1 + 0.7
      f12_cv = rep(0, length(th_cv_list))
      f12_cv_se = rep(0, length(th_cv_list))
  
      for (k in 1:length(th_cv_list)){
          thc = th_cv_list[k]
            
          lambda1=lambdac2*sqrt(d*log(d)/n)/(d^2)
          th=thc*sqrt(n/(d*log(d)))/d
          obj = CV_thresh(Y,X,th,lambda1)
          f12_cv[k] = obj[1]
          f12_cv_se[k] = obj[2]
          
          cat("======= thc CV",  k, "=======", "\n")
      }    
      
      k_min = which(f12_cv == min(f12_cv))
      kvec = which(f12_cv < f12_cv[k_min] + f12_cv_se[k_min])
      k_min_1sd = kvec[1]
      thc = th_cv_list[k_min_1sd]
      #thc = th_cv_list[k_min]
      para[j,2] = thc
  
      
      for (i in 1:6){
          n=2000*i
          #epsilon=randomt(n, 1)/(64*d)
          #sigma = 3
          #epsilon=(exp(rnorm(n, 0, sigma))-exp(sigma^2/2)*rep(1,n))/(1000*d)
          sigma=.5
          epsilon=rnorm(n, mean=0, sd=sigma)/d
          obs=sample(x=1:(d^2), size=n, replace=T)
          X=matrix(0, n, 2)    
          X[, 1]=(obs-.1) %/% d+1
          X[, 2]=obs %% d
          X[X[, 2]==0, 2]=d
          Y=Theta_t[obs]+epsilon
              
          lambda=lambdac1*sqrt(d*log(d)/n)/(d^2)
          Theta=admm(Y=Y, X=X, rho=.1, mu=lambda, alpha=alpha)
          f1[j,i]=norm(Theta-Theta_t, "f")
              
          lambda1=lambdac2*sqrt(d*log(d)/n)/(d^2)
          th=thc*sqrt(n/(d*log(d)))/d
          Theta1=admm(Y=truncate_abs(Y, th), X=X, rho=.1, mu=lambda1, alpha=alpha)
          f2[j,i]=norm(Theta1-Theta_t, "f")
              
          cat("=======", j, ",",  i, "=======", "\n")
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

fn = paste("/mhome/stats/t/zz385/robust/data_042019/mc_CV_cauchy_", seed, ".RData", sep="")
save(f11, f12, f21, f22, f31, f32, f1_para, f2_para, f3_para, file=fn)
save(f11, f12, f21, f22, f31, f32, f1_para, f2_para, f3_para, file="mc_CV.RData")

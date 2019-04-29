setwd("C:/Users/weichenw/Dropbox/robust/code_2019_Weichen")
source("utl.R")
set.seed(111)

N = 3
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
    r0=5    
    n0=100
    Z=matrix(rnorm(n0*d,0,1), n0, d)    
    S=t(Z)%*%Z
    V=svd(S, nu=r0, nv=r0)$u[,1:r0]
    Theta_t=as.vector(V%*%t(V))
    return(Theta_t)
}


CV = function(Y,X,lambda){
    fold = 5
    nn = floor(length(Y) / fold)
    CV_error = rep(0, fold)
    for (h in 1:fold) {
        if (h < fold)  testSam = 1:nn + (h-1)*nn
        else  testSam = ((fold-1)*nn+1):(length(Y))
        Y_train = Y[-testSam]
        X_train = X[-testSam,]
        lambda_train = lambda * sqrt(length(Y)) / sqrt(length(Y_train))
        
        Theta = prsm(Y_train, X_train, alpha=.9, beta=1, lambda=lambda_train)
        Y_th = quantile(abs(Y[testSam]),0.95)
        Y_target = truncate_abs(Y[testSam], Y_th) 
        CV_error[h] = sum((Y_target - X[testSam,] %*% Theta)^2) / length(testSam)
    }
    return(c(mean(CV_error), sd(CV_error) / sqrt(fold)))
}


CV_thresh = function(Y,X,th,lambda1){
    fold = 5
    nn = floor(length(Y) / fold)
    CV_error = rep(0, fold)
    for (h in 1:fold) {
        if (h < fold)  testSam = 1:nn + (h-1)*nn
        else  testSam = ((fold-1)*nn+1):(length(Y))
        Y_train = Y[-testSam]
        X_train = X[-testSam,]
        lambda_train = lambda1 * sqrt(length(Y)) / sqrt(length(Y_train))
        th_train = th/sqrt(length(Y))*sqrt(length(Y_train))
        
        Theta = prsm(truncate_abs(Y_train, tau=th_train), X_train, alpha=.9, beta=1, lambda=lambda_train)
        Y_th = quantile(abs(Y[testSam]),0.95)
        Y_target = truncate_abs(Y[testSam], Y_th) 
        CV_error[h] = sum((Y_target - X[testSam,] %*% Theta)^2) / length(testSam)
    }
    return(c(mean(CV_error), sd(CV_error) / sqrt(fold)))
}


for (d in c(20, 40, 60)) {
    cat("### d =", d, "### \n")
    
    f1 = matrix(0, N, 6)
    f2 = matrix(0, N, 6)
    para = matrix(0, N, 2)
  
    for (j in 1:N){
        n=300
        r=5
        Theta_t = genTrueTheta(d)
        
        sigma=.5
        epsilon=rnorm(n, mean=0, sd=sigma)
        #epsilon=randomt(n, 1)/64
        #sigma=2.5
        #epsilon=(exp(rnorm(n, 0, sigma))-exp(sigma^2/2)*rep(1,n))/1000
        X=matrix(rnorm(n*d*d, 0, 1), n, d^2)
        Y=X%*%Theta_t+epsilon
        
        lambda_cv_list = c(0:10) * 0.2 + 0.2
        f11_cv = rep(0, length(lambda_cv_list))
        f11_cv_se = rep(0, length(lambda_cv_list))
        
        for (k in 1:length(lambda_cv_list)){
            lambdac = lambda_cv_list[k]
            
            lambda = lambdac * sqrt(d / n)
            obj = CV(Y,X,lambda)
            f11_cv[k] = obj[1]
            f11_cv_se[k] = obj[2]
            
            cat("======= lambdac CV",  k, "=======", "\n")
        }
        
        k_min = which(f11_cv == min(f11_cv))
        if (length(k_min) > 1) k_min_1sd = k_min[1]
        else {
            kvec = which(f11_cv < f11_cv[k_min] + f11_cv_se[k_min])
            k_min_1sd = kvec[length(kvec)]
        }
        lambdac1 = lambda_cv_list[k_min_1sd]
        lambdac2 = lambdac1
        para[j,1] = lambdac1
            
        
        th_cv_list = c(0:10) * 0.1 + 0.7
        f12_cv = rep(0, length(th_cv_list))
        f12_cv_se = rep(0, length(th_cv_list))
        
        for (k in 1:length(th_cv_list)){
            thc = th_cv_list[k]
            
            lambda1=lambdac2*sqrt(d/n)
            th=thc*sqrt(n/d)
            obj = CV_thresh(Y,X,th,lambda1)
            f12_cv[k] = obj[1]
            f12_cv_se[k] = obj[2]
            
            cat("======= thc CV",  k, "=======", "\n")
        } 
            
        k_min = which(f12_cv == min(f12_cv))
        if (length(k_min) > 1) k_min_1sd = k_min[length(k_min)]
        else {
            kvec = which(f12_cv < f12_cv[k_min] + f12_cv_se[k_min])
            k_min_1sd = kvec[1]
        }
        thc = th_cv_list[k_min_1sd]
        para[j,2] = thc
        
        for (i in 1:6){
            n=300*i
            sigma=.5
            epsilon=rnorm(n, mean=0, sd=sigma)
            #epsilon=randomt(n, 1)/64
            #sigma=2.5
            #epsilon=(exp(rnorm(n, 0, sigma))-exp(sigma^2/2)*rep(1,n))/1000
            X=matrix(rnorm(n*d*d, 0, 1), n, d^2)
            Y=X%*%Theta_t+epsilon
              
            lambda1=lambdac1*sqrt(d/n)
            Theta1=prsm(Y, X, alpha=.9, beta=1, lambda=lambda1)
            f1[j,i]=sqrt(sum((Theta1-Theta_t)^2))
            
            lambda2=lambdac2*sqrt(d/n)
            th=thc*sqrt(n/d)
            Theta2=prsm(truncate_abs(Y, tau=th), X, alpha=.9, beta=1, lambda=lambda2)        
            f2[j,i]=sqrt(sum((Theta2-Theta_t)^2))
            
            cat("=======", j, ",", i, "=======", "\n")
        }
    }
    
    print(f1)
    print(f2)
    print(para)
    
    if (d == 20) {
        f11 = f1
        f12 = f2
        f1_para = para
    }
    if (d == 40) {
        f21 = f1
        f22 = f2
        f2_para = para
    }
    if (d == 60) {
        f31 = f1
        f32 = f2
        f3_para = para
    }
}

save(f11, f12, f21, f22, f31, f32, f1_para, f2_para, f3_para, file="cs_CV_normal.RData")





set.seed(222)

N = 10
g11 = matrix(0, N, 6)
g12 = matrix(0, N, 6)
g21 = matrix(0, N, 6)
g22 = matrix(0, N, 6)
g31 = matrix(0, N, 6)
g32 = matrix(0, N, 6)



d=20
lambdac1 = mean(f1_para[,1])
lambdac2 = lambdac1
thc = mean(f1_para[,2])
for (i in 1:6){
    r=5
    n=300*i
    for (j in 1:N){
        Theta_t = genTrueTheta(d)
        sigma=.5
        epsilon=rnorm(n, mean=0, sd=sigma)
        #epsilon=randomt(n, 1)/64
        #sigma=2.5
        #epsilon=(exp(rnorm(n, 0, sigma))-exp(sigma^2/2)*rep(1,n))/1000
        X=matrix(rnorm(n*d*d, 0, 1), n, d^2)
        Y=X%*%Theta_t+epsilon
        
        lambda1=lambdac1*sqrt(d/n)
        Theta1=prsm(Y, X, alpha=.9, beta=1, lambda=lambda1)
        g11[j,i]=sqrt(sum((Theta1-Theta_t)^2))
        
        lambda2=lambdac2*sqrt(d/n)
        th=thc*sqrt(n/d)
        Theta2=prsm(truncate_abs(Y, tau=th), X, alpha=.9, beta=1, lambda=lambda2)        
        g12[j,i]=sqrt(sum((Theta2-Theta_t)^2))

        cat("=======", i, ",", j, "=======", "\n")
    }    
}

d=40
lambdac1 = mean(f2_para[,1])
lambdac2 = lambdac1
thc = mean(f2_para[,2])
for (i in 1:6){
    r=5
    n=300*i
    for (j in 1:N){
        Theta_t = genTrueTheta(d)
        sigma=.5
        epsilon=rnorm(n, mean=0, sd=sigma)
        #epsilon=randomt(n, 1)/64
        #sigma=2.5
        #epsilon=(exp(rnorm(n, 0, sigma))-exp(sigma^2/2)*rep(1,n))/1000
        X=matrix(rnorm(n*d*d, 0, 1), n, d^2)
        Y=X%*%Theta_t+epsilon
        
        lambda1=lambdac1*sqrt(d/n)
        Theta1=prsm(Y, X, alpha=.9, beta=1, lambda=lambda1)
        g21[j,i]=sqrt(sum((Theta1-Theta_t)^2))
        
        lambda2=lambdac2*sqrt(d/n)
        th=thc*sqrt(n/d)
        Theta2=prsm(truncate_abs(Y, tau=th), X, alpha=.9, beta=1, lambda=lambda2)        
        g22[j,i]=sqrt(sum((Theta2-Theta_t)^2))
        
        cat("=======", i, ",", j, "=======", "\n")
    }    
}

d=60
lambdac1 = mean(f3_para[,1])
lambdac2 = lambdac1
thc = mean(f3_para[,2])
for (i in 1:6){
    r=5
    n=300*i
    for (j in 1:N){
        Theta_t = genTrueTheta(d)
        sigma=.5
        epsilon=rnorm(n, mean=0, sd=sigma)
        #epsilon=randomt(n, 1)/64
        #sigma=2.5
        #epsilon=(exp(rnorm(n, 0, sigma))-exp(sigma^2/2)*rep(1,n))/1000
        X=matrix(rnorm(n*d*d, 0, 1), n, d^2)
        Y=X%*%Theta_t+epsilon
        
        lambda1=lambdac1*sqrt(d/n)
        Theta1=prsm(Y, X, alpha=.9, beta=1, lambda=lambda1)
        g31[j,i]=sqrt(sum((Theta1-Theta_t)^2))
        
        lambda2=lambdac2*sqrt(d/n)
        th=thc*sqrt(n/d)
        Theta2=prsm(truncate_abs(Y, tau=th), X, alpha=.9, beta=1, lambda=lambda2)        
        g32[j,i]=sqrt(sum((Theta2-Theta_t)^2))
        
        cat("=======", i, ",", j, "=======", "\n")
    }    
}

save(g11, g12, g21, g22, g31, g32, file="cs_normal.RData")


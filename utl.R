prsm_multi=function(Y, X, alpha, beta, lambda){
    error=10^-3
    Theta_x=matrix(0, d, d)
    Theta_y=matrix(0, d, d)
    P=matrix(0, d, d)
    iter=0
    maxiter=10^3
    repeat{
        Theta_x=solve(t(X)%*%X/n+beta*diag(d))%*%(beta*Theta_y+P+t(X)%*%Y/n)
        P=P-alpha*beta*(Theta_x-Theta_y)
        W=Theta_x-P/beta
        tmp=svd(W)
        Theta_y1=tmp$u %*% diag(soft(tmp$d, tau=lambda/beta)) %*% t(tmp$v)		
        P=P-alpha*beta*(Theta_x-Theta_y)
        if (norm(Theta_y1-Theta_y, "F")<error){
            return(Theta_y)
            break
        }
        iter=iter+1
        if (iter>maxiter){
            return(NA)
            break
        }
        # if (iter %% 10==0){
        # cat("=====", iter/10, "=====", "\n")
        # }
        Theta_y=Theta_y1
    }
}


prsm = function(Y, X, alpha, beta, lambda){
    # This PRSM algorithm is for compressed sensing and multi_task regression. 
    # alpha, beta: some fixed parameters 
    # lambda: penalty coefficient 


    error = 10 ^ -3    
    Theta_x = rep(0, d ^ 2)
    Theta_y = rep(0, d ^ 2)   
    P = rep(0, d ^ 2)
    iter = 0
    maxiter = 10 ^ 3
    n = length(Y)
    repeat{
        if (d ^ 2 < n){
            Theta_x = solve(2 * t(X) %*% X / n + beta * diag(d ^ 2)) %*% (beta * Theta_y + P + 2 * t(X) %*% Y / n)       
        }
        else {
            # X=matrix(rnorm(n*d*d, 0, 1), n, d^2)
            # Y=rnorm(n, mean=0, sd=1)
            # t1=proc.time()
            S = svd(X, nu=0, nv=n)
            Theta_x = (S$v %*% diag((2 * (S$d) ^ 2 / n + beta * rep(1,n)) ^ (-1) - beta^(-1) * rep(1,n)) %*% t(S$v) 
                + (beta ^ -1) * diag(d ^ 2) ) %*% (beta * Theta_y + P + 2 * t(X) %*% Y / n)
            # t2=proc.time()-t1
            # t1=proc.time()
            # Theta_x=solve(2*t(X)%*%X/n+beta*diag(d^2))%*%(beta*Theta_y+P+2*t(X)%*%Y/n)
            # t2=proc.time()-t1
        }
        # Theta_x=solve(2*t(X)%*%X/n+beta*diag(d^2))%*%(beta*Theta_y+P+2*t(X)%*%Y/n)
        P = P - alpha * beta * (Theta_x - Theta_y)
        W = matrix(Theta_x - P / beta, d, d)
        tmp = svd(W)
        Theta_y1 = as.vector(tmp$u %*% diag(soft(tmp$d, tau=lambda / beta)) %*% t(tmp$v))      
        P = P - alpha * beta * (Theta_x - Theta_y1)
        if (sqrt(sum((Theta_y1 - Theta_y) ^ 2)) < error){
            return(Theta_y)
            break
        }
        iter = iter + 1
        if (iter > maxiter){
            return(NA)
            break
        }   
        # cat("=====", sqrt(sum((Theta_y1-Theta_y)^2)) , "=====", "\n")
        Theta_y=Theta_y1
    }
}


admm = function(Y, X, rho, mu, alpha){
    # This ADMM algorithm is for matrix completion. 
    # mu: coeffcient in front of the nuclear-norm penalty
    # alpha: elementwise absolute value bound


    tau = 1.618
    iter = 0
    nTheta = matrix(0, d, d)
    sTheta = matrix(0, d, d)
    n = length(Y)
    for (i in 1:n){
        nTheta[X[i,1], X[i,2]] = nTheta[X[i,1], X[i,2]] + 1
        sTheta[X[i,1], X[i,2]] = sTheta[X[i,1], X[i,2]] + Y[i]
    }
    Z = matrix(0, 2 * d, 2 * d)
    W = matrix(0, 2 * d, 2 * d)
    repeat{
        s = eigen(Z - (rho ^ -1) * (W + mu * diag(2 * d)))
        V = s$vectors %*% diag(pmax(s$values, 0)) %*% t(s$vectors)
        C = V + W / rho
        C12 = C[1:d, (d + 1):(2 * d)]
        Z12 = (rho * C12 + 2 / n * sTheta) / (rho + 2 / n * nTheta)
        supp = abs(Z12) > alpha
        Z12[supp] = Z12[supp] / abs(Z12[supp]) * alpha
        Z.new = C
        Z.new[1:d, (d + 1):(2 * d)] = Z12
        Z.new[(d+1):(2*d), 1:d]=t(Z12)
        W = W + tau * rho * (V - Z.new)
        dW = (tau - 1) * rho * (V - Z.new)
        Rp = norm(Z.new - V, "f")
        Rd = max(norm(rho * (Z - Z.new) + dW, "f"), norm(dW, "f") )
        #Rd = max(norm(tau*rho*(V-Z.new),"f"),norm(Z-Z.new,"f"))
        iter = iter + 1
        if (max(Rp, Rd) < 10 ^ -5){
            break
        }
        if (iter > 500){
            break
        }
        if (iter %% 10 == 0){
            if (Rp < .5 * Rd) 
                rho = .7 * rho
            if (Rd< .5 * Rp) 
                rho = 1.3 * rho
        }
        Z = Z.new
        # cat("=====", iter, " ", Rp, " ", Rd, "=====", "\n")
    }
    return(Z[1:d, (d + 1):(2 * d)])
}


truncate_abs = function(Y, tau){
    # Truncate the absolute value of each entry of Y by tau

    Y1 = Y
    supp = abs(Y1) > tau
    Y1[supp] = Y1[supp] / abs(Y1[supp]) * tau
    return(Y1)
}


truncate_l4 = function(Y, tau){
    # Shrink each row of Y by its l4 norm. 


    nrow = dim(Y)[1]
    ncol = dim(Y)[2]
    Y1 = matrix(0, nrow, ncol)
    for (i in 1:nrow){
        y = Y[i, ]
        if (sum(y ^ 4)^(1 / 4) >= tau){
            Y1[i, ] = y / ((sum(y ^ 4)) ^ (1 / 4)) * tau
        }
        else Y1[i, ] = y
    }
    return(Y1)
}


truncate_l2 = function(Y, tau){
    # Truncate the l2-norm of each row of Y by tau. 


    nrow = dim(Y)[1]
    ncol = dim(Y)[2]
    Y1 = matrix(0, nrow, ncol)
    for (i in 1:nrow){
        y = Y[i, ]
        if (sqrt(sum(y ^ 2)) >= tau){
            Y1[i, ] = y / (sqrt(sum(y ^ 2))) * tau
        }
        else Y1[i, ] = y
    }
    return(Y1)
}



soft = function(x, tau){
    # Apply soft-thresholding to x by tau. 


    sgn = sign(x)
    x = abs(x) - tau
    x[x < 0] = 0
    return(sgn * x)   
}


randomt = function(n, factor){
# Generate n independent Cauchy variables. 


    y = sapply(1:n, function (o){
        repeat{
            x = rt(1, df=1) / factor
            if (abs(x) < 10 ^ 4) break
        }
        return(x)
    })
    return(y)
}

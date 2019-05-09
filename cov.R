Distribution = 't'
#Distribution = 'normal'


#### find constant C
set.seed(111)

getEta = function(vals, tau) {
  return(mean(pmin(tau / vals, 1))) 
}

solveForTau_inner = function(vals,eta,lo,hi) {
  mid = (lo+hi) / 2
  etaNow = getEta(vals, mid)
  if (abs(etaNow - eta) < 0.01) {
    return(mid)
  } else {
    if (etaNow > eta) return(solveForTau_inner(vals, eta, lo, mid))
    else return(solveForTau_inner(vals, eta, mid, hi))
  }
}

solveForTau = function(vals, eta) {
  lo = 0
  hi = max(vals)
  return(solveForTau_inner(vals,eta,lo,hi))
}

n = 100
d = 100
N = 100
C_select = matrix(0,3,N)
for (i in 1:N){
  if (Distribution == 't') {
    X = sapply(1:d, function(o) rt(n, df=3)/sqrt(3))
  } else {
    X = sapply(1:d, function(o) rnorm(n))
  }
  D = diag(c(2,  rep(1, d - 1)))
  Y = X %*% D

  for (k in 1:3) {
    eta = 0.7 + 0.1*(k-1)
    tauTarget = solveForTau(apply(Y,1,function(a) sum(a^4)^(1/4)), eta)
    C_select[k,i]= tauTarget * (log(d)/n)^(1/4)
  }
}
C_final = apply(C_select, 1, mean)
print(C_final)


### cov estimation
set.seed(111)

truncate=function(Y, tau){
  nrow = dim(Y)[1]
  ncol = dim(Y)[2]
  Y1 = matrix(0, nrow, ncol)
  dataRatio = 0
  for (i in 1:nrow){
    y=Y[i,]
    if (sum(y^4)^(1/4)>=tau){
      ratio = tau / ((sum(y^4))^(1/4))
      dataRatio = dataRatio + ratio
      Y1[i,]=y * ratio
    }
    else {
      dataRatio = dataRatio + 1
      Y1[i,]=y
    }
  }
  return(list(Y1, dataRatio / n))
}

errtot1 = matrix(0, nrow=3, ncol=5)
errtot2_7 = matrix(0, nrow=3, ncol=5)
errtot2_8 = matrix(0, nrow=3, ncol=5)
errtot2_9 = matrix(0, nrow=3, ncol=5)
dataEffPerc_7 = matrix(0, nrow=3, ncol=5)
dataEffPerc_8 = matrix(0, nrow=3, ncol=5)
dataEffPerc_9 = matrix(0, nrow=3, ncol=5)

N = 100
ratio.all = c(.2, .5, 1)

for (k in 1:3){
    
    ratio = ratio.all[k]
        
    for (j in 1:5){

        n = 100 * j
        d = ratio * n
        D = diag(c(2,  rep(1, d - 1)))
        S = D %*% D
        
        error1 = rep(0, N)
        error2_7 = rep(0, N)
        error2_8 = rep(0, N)
        error2_9 = rep(0, N)
        dataEffRatio_7 = rep(0,N)
        dataEffRatio_8 = rep(0,N)
        dataEffRatio_9 = rep(0,N)
        for (i in 1:N){
            if (Distribution == 't')
              X = sapply(1:d, function(o) rt(n, df=3)/sqrt(3))
            else
              X = sapply(1:d, function(o) rnorm(n))
            
            Y = X %*% D

            Sn = t(Y) %*% Y / n
            error1[i] = norm(Sn - S, "2")
            
            for (ci in 1:3) {
              c0 = C_final[ci]
              tau = c0 * (n / log(d)) ^ (1 / 4)
              li = truncate(Y, tau)
              Y_trunc = li[[1]]
              Sn_trunc = t(Y_trunc) %*% Y_trunc/n
              
              if (ci == 1) {
                dataEffRatio_7[i] = li[[2]]
                error2_7[i]= norm(Sn_trunc - S, "2")
              }
              if (ci == 2) {
                dataEffRatio_8[i] = li[[2]]
                error2_8[i]= norm(Sn_trunc - S, "2")
              }
              if (ci == 3) {
                dataEffRatio_9[i] = li[[2]]
                error2_9[i]= norm(Sn_trunc - S, "2")
              }
            }
        }
        errtot1[k, j] = mean(error1)
        errtot2_7[k, j] = mean(error2_7)
        errtot2_8[k, j] = mean(error2_8)
        errtot2_9[k, j] = mean(error2_9)
        dataEffPerc_7[k, j] = mean(dataEffRatio_7)
        dataEffPerc_8[k, j] = mean(dataEffRatio_8)
        dataEffPerc_9[k, j] = mean(dataEffRatio_9)
        print(j)
    }
}

save(errtot1, errtot2_7,errtot2_8,errtot2_9, dataEffPerc_7,dataEffPerc_8,dataEffPerc_9, file=paste("cov_", Distribution,".RData",sep=''))





lo = 0 # min(min(errtot1),min(errtot2))
hi = max(max(errtot1),max(errtot2_7),max(errtot2_8),max(errtot2_9))

pdf(paste("cov_", Distribution,".pdf",sep='')) 
#plot.new()
for (i in 1:3){
  e1=errtot1[i,]
  e2_7=errtot2_7[i,]
  e2_8=errtot2_8[i,]
  e2_9=errtot2_9[i,]
  if (i==1){
    plot((1:5)*100, e1, lwd=1, col=1+i, pch=0, type="b", ylim=c(lo,hi), ylab="Error in Operator Norm", xlab="n")	
  } else {
    lines((1:5)*100, e1, lwd=1, col=1+i, pch=0, type="b")
  }
  lines((1:5)*100, e2_8, lwd=1, col=1+i, pch=3, type="b")
  lines((1:5)*100, e2_7, lwd=1, col=1+i, pch=3, type="b",lty='dotted')
  lines((1:5)*100, e2_9, lwd=1, col=1+i, pch=3, type="b",lty='dotdash')
}

if (Distribution == 't') {
  legend("topleft", legend=c("d=0.2n, Standard", "d=0.2n, Robust",
                             "d=0.5n, Standard", "d=0.5n, Robust",
                             "d=n, Standard", "d=n, Robust"), col=c(2,2,3,3,4,4), pch=c(0,3,0,3,0,3))
} else {
  legend("bottomright", legend=c("d=0.2n, Standard", "d=0.2n, Robust",
                             "d=0.5n, Standard", "d=0.5n, Robust",
                             "d=n, Standard", "d=n, Robust"), col=c(2,2,3,3,4,4), pch=c(0,3,0,3,0,3))
}
dev.off()



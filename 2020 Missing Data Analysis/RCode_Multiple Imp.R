

################# theoretical simulation data ########################
library(MASS)
set.seed(1)
r = matrix(runif(n = 9^2,min = -1,max=1), ncol=9)
Cor = cov2cor(r%*%t(r))
aaa= mvrnorm(n = 1000, mu = rep(0,9), Sigma = Cor) # 900*9 matrix
X = matrix(c(rep(1,1000),aaa),ncol=10)

set.seed(123)
beta = runif(n = 10, min = -2, max = 2)
sig = runif(n = 1, min = 0.5, max = 2)
eps= rnorm(n=1000 , 0, sig)
y =matrix(rep(0,1000*9),ncol=9)
for(i in 2:10){
  y[,i-1] = X[,1:i]%*%beta[1:i]+eps
}
set.seed(3)
u = matrix(runif(9000,0,1),ncol=9)
X_mis = X
for(i in 2:10){
  X_mis[which(u[,i-1]<0.5 & y[,i-1]<quantile(y[,i-1],0.5) | u[,i-1]>0.9 & y[,i-1]>quantile(y[,i-1],0.5)),i] = NA
}

summary(X_mis)

################ check MAR ###############################
par(mfrow=c(1,1))
p1 =hist(y[,8], breaks = 20, probability = TRUE)
p2 = hist(y[-which(is.na(X_mis[,9])),8], probability = TRUE, breaks = 20)
plot(p1, col=rgb(0,0,1,1/4), main="Change of Distribution in the 8th Variable")
plot(p2, col=rgb(1,0,0,1/4),add=T)
wilcox.test(y[,1], y[-which(is.na(X_mis[,2])),1])

cor(X)
cor(na.omit(X_mis[,1:10]))

################ CD ########################
CD = function(X){
  set.seed(1)
  r = matrix(runif(n = 9^2,min = -1,max=1), ncol=9)
  Cor = cov2cor(r%*%t(r))
  seCD = ciCD = betaCD = tempCD = rmseCD = tempse = cicovCD = cilenCD = biasCD = tempbias = NULL
  for(i in 1:9){
    print(i)
    for (n in 1:500){
      set.seed(n)
      
      aaa= mvrnorm(n = 500, mu = rep(0,9), Sigma = Cor) # 900*6 matrix
      X = matrix(c(rep(1,500),aaa),ncol=10)
      
      beta = runif(n = 10, min = -2, max = 2)
      sig = runif(n = 1, min = 0.5, max = 2)
      eps= rnorm(n=500 , 0, sig)
      y =matrix(rep(0,500*9),ncol=9)
      for(kk in 2:10){
        y[,kk-1] = X[,1:kk]%*%beta[1:kk]+eps
      }
      
      u = matrix(runif(4500,0,1),ncol=9)
      betaCD[n] = coef(lm(y[,i]~X[,2:c(i+1)]))[2]
      seCD[n] = as.matrix(summary(lm(y[,i]~X[,2:c(i+1)]))$coefficients)[2,2]
      ciCD[n] = ifelse(betaCD[n]-1.96*seCD[n]<= beta[2] & betaCD[n]+1.96*seCD[n] >= beta[2],1,0)
      tempbias[n] = betaCD[n] - beta[2]
    }
    biasCD[i] = mean(tempbias)
    rmseCD[i] = sqrt(mean((tempbias)^2))
    cicovCD[i] = mean(ciCD)
    cilenCD[i] = 2*1.96*mean(seCD)
  }
  return(data.frame(biasCD = biasCD, rmseCD = rmseCD, cicovCD = cicovCD, cilenCD = cilenCD))
}

CD = CD(X_mis)
write.csv(t(CD), 'CD.csv')

############### LD #########################

LD = function(X_mis){
  set.seed(1)
  r = matrix(runif(n = 9^2,min = -1,max=1), ncol=9)
  Cor = cov2cor(r%*%t(r))
  seLD = ciLD = betaLD = tempLD = rmseLD = tempse = cicovLD = cilenLD = temprmse = biasLD = tempbias = NULL
  for(i in 1:9){
    print(i)
    for (n in 1:500){
      set.seed(n)
      
      aaa= mvrnorm(n = 500, mu = rep(0,9), Sigma = Cor) # 900*6 matrix
      X = matrix(c(rep(1,500),aaa),ncol=10)
      
      beta = runif(n = 10, min = -2, max = 2)
      sig = runif(n = 1, min = 0.5, max = 2)
      eps= rnorm(n=500 , 0, sig)
      y =matrix(rep(0,500*9),ncol=9)
      for(kk in 2:10){
        y[,kk-1] = X[,1:kk]%*%beta[1:kk]+eps
      }
      
      u = matrix(runif(4500,0,1),ncol=9)
      X_mis = X
      
      for(k in 2:10){
        X_mis[which(u[,k-1]<0.5 & y[,k-1]<quantile(y[,k-1],0.5) | u[,k-1]>0.9 & y[,k-1]>quantile(y[,k-1],0.5)),k] = NA
      }
      betaLD[n] = coef(lm(y[,i]~X_mis[,2:c(i+1)],na.action=na.omit))[2]
      seLD[n] = as.matrix(summary(lm(y[,i]~X_mis[,2:c(i+1)]))$coefficients)[2,2]
      ciLD[n] = ifelse(betaLD[n]-1.96*seLD[n]<= beta[2] & betaLD[n]+1.96*seLD[n] >= beta[2],1,0)
      tempbias[n] = betaLD[n]- beta[2]
    }
    biasLD[i] = mean(tempbias) 
    rmseLD[i] = sqrt(mean((tempbias)^2))
    cicovLD[i] = mean(ciLD)
    cilenLD[i] = 2*1.96*mean(seLD)
  }
  return(data.frame(biasLD = biasLD, rmseLD = rmseLD, cicovLD = cicovLD, cilenLD = cilenLD))
}

LD = LD(X_mis)
write.csv(t(LD), 'LD.csv')

########SSI##############
library(mice)

SSI = function(X_mis){
  set.seed(1)
  r = matrix(runif(n = 9^2,min = -1,max=1), ncol=9)
  Cor = cov2cor(r%*%t(r))
  seSSI = ciSSI = time = comptime = betaSSI = tempSSI = rmseSSI = tempse = cicovSSI = cilenSSI = temprmse = biasSSI = tempbias = NULL
  for(i in 1:9){
    print(i)
    for (n in 1:500){
      set.seed(n)
      
      aaa= mvrnorm(n = 500, mu = rep(0,9), Sigma = Cor) # 900*6 matrix
      X = matrix(c(rep(1,500),aaa),ncol=10)
      
      beta = runif(n = 10, min = -2, max = 2)
      sig = runif(n = 1, min = 0.5, max = 2)
      eps= rnorm(n=500 , 0, sig)
      y =matrix(rep(0,500*9),ncol=9)
      for(kk in 2:10){
        y[,kk-1] = X[,1:kk]%*%beta[1:kk]+eps
      }
      
      u = matrix(runif(4500,0,1),ncol=9)
      X_mis = X
      
      for(k in 2:10){
        X_mis[which(u[,k-1]<0.5 & y[,k-1]<quantile(y[,k-1],0.5) | u[,k-1]>0.9 & y[,k-1]>quantile(y[,k-1],0.5)),k] = NA
      }
      start = Sys.time()      
      imp = mice(data.frame(y[,i],X_mis[,2:c(i+1)]), m=1, printFlag=FALSE, method="norm.nob", seed=n)
      Xdep = as.matrix(complete(imp)[-1])
      betaSSI[n] = coef(lm(y[,i]~Xdep))[2]
      seSSI[n] = as.matrix(summary(lm(y[,i]~Xdep))$coefficients)[2,2]
      ciSSI[n] = ifelse(betaSSI[n]-1.96*seSSI[n]<= beta[2] & betaSSI[n]+1.96*seSSI[n] >= beta[2],1,0)
      end = Sys.time()
      time[n] = as.numeric(end-start)
      tempbias[n] = betaSSI[n]- beta[2]
    }
    comptime[i] = mean(time)
    biasSSI[i] = mean(tempbias) 
    rmseSSI[i] = sqrt(mean((tempbias)^2))
    cicovSSI[i] = mean(ciSSI)
    cilenSSI[i] = 2*1.96*mean(seSSI)
  }
  return(data.frame(biasSSI = biasSSI, rmseSSI = rmseSSI, cicovSSI = cicovSSI, cilenSSI = cilenSSI))
}
SSI = SSI(X_mis)
write.csv(t(SSI), 'SSI.csv')


############## EMB ########################
library(Amelia)




EMB = function(X_mis){
  set.seed(1)
  r = matrix(runif(n = 9^2,min = -1,max=1), ncol=9)
  Cor = cov2cor(r%*%t(r))
  seEM = ciEM = time = comptime = betaEM = tempEM = rmseEM = tempse = cicovEM = cilenEM = temprmse = biasEM = tempbias = NULL
  for(i in 1:9){
    print(i)
    for (n in 1:500){
      set.seed(n)
      
      aaa= mvrnorm(n = 500, mu = rep(0,9), Sigma = Cor) # 900*6 matrix
      X = matrix(c(rep(1,500),aaa),ncol=10)
      
      beta = runif(n = 10, min = -2, max = 2)
      sig = runif(n = 1, min = 0.5, max = 2)
      eps= rnorm(n=500 , 0, sig)
      y =matrix(rep(0,500*9),ncol=9)
      for(kk in 2:10){
        y[,kk-1] = X[,1:kk]%*%beta[1:kk]+eps
      }
      
      u = matrix(runif(4500,0,1),ncol=9)
      X_mis = X
      
      for(k in 2:10){
        X_mis[which(u[,k-1]<0.5 & y[,k-1]<quantile(y[,k-1],0.5) | u[,k-1]>0.9 & y[,k-1]>quantile(y[,k-1],0.5)),k] = NA
      }
      start = Sys.time()      
      emb = amelia(data.frame(y[,i],X_mis[,2:c(i+1)]), m=5,p2s=0)
      for(j in 1:5){
        Xdep = as.matrix(emb$imputations[[j]][-1])
        tempEM[j] = coef(lm(y[,i]~Xdep))[2]
        tempse[j] = as.matrix(summary(lm(y[,i]~Xdep))$coefficients)[2,2]
      }
      betaEM[n] = mean(tempEM)
      seEM[n] = sqrt(mean(tempse^2)+var(tempEM)*(1+1/5))
      ciEM[n] = ifelse(betaEM[n]-1.96*seEM[n] <= beta[2] & betaEM[n]+1.96*seEM[n] >= beta[2],1,0)
      end = Sys.time()
      time[n] = as.numeric(end-start)
      tempbias[n] = betaEM[n] - beta[2]
      temprmse[n] = sqrt(mean((tempEM-beta[2])^2))
    }
    comptime[i] = mean(time)
    biasEM[i] = mean(tempbias)
    rmseEM[i] = mean(temprmse)
    cicovEM[i] = mean(ciEM)
    cilenEM[i] = 2*1.96*mean(seEM)
  }
  return(data.frame(biasEM = biasEM, rmseDA = rmseEM, cicovEM = cicovEM, cilenEM = cilenEM, comptime = comptime))
}
EMB = EMB(X_mis)
write.csv(t(EMB), 'EMB.csv')


summary(emb)

library(mice)


# Imputed Values
emb = amelia(data.frame(y[,9],X_mis[,2:10]), m=10)
hist(emb$imputations[[1]][,9], col="grey", border="white", main = "Imputed Values")


# MCMC(DA)
p = ncol(y)
m = 5
t = 20
mc = 500
type = 1
str(y)
summary(lm(y[,1]~X[,2]))

i=j=1
DA = function(m = 5, t = num_iter, mc = 500, type = 1){
  set.seed(1)
  r = matrix(runif(n = 9^2,min = -1,max=1), ncol=9)
  Cor = cov2cor(r%*%t(r))
  p = ncol(r)
  seDA = ciDA = time = comptime = betaDA = tempDA = rmseDA = tempse = cicovDA = cilenDA = biasDA = NULL
  iter_biasDA = iter_rmse = NULL
  for (i in 1:p) {
    for (j in 1:mc){
      set.seed(j)
      aaa= mvrnorm(n = 1000, mu = rep(0,9), Sigma = Cor) # 900*6 matrix
      X = matrix(c(rep(1,1000),aaa),ncol=10)
      # set.seed(j)
      beta = runif(n = 10, min = -2, max = 2)
      # set.seed(j)
      sig = runif(n = 1, min = 0.5, max = 2)
      # set.seed(j)
      eps= rnorm(n=1000 , 0, sig)
      y = matrix(rep(0,1000*9),ncol=9)
      for(kk in 2:10){
        y[,kk-1] = X[,1:kk]%*%beta[1:kk]+eps
      }
      
      u = matrix(runif(9000,0,1),ncol=9)
      X_mis = X
      
      for(k in 2:10){
        X_mis[which(u[,k-1]<0.5 & y[,k-1]<quantile(y[,k-1],0.5) | u[,k-1]>0.9 & y[,k-1]>quantile(y[,k-1],0.5)),k] = NA
      }
      
      if(type == 1){
        init_beta=matrix(rep(0,i+1), nrow=1, ncol=i+1)
        init_sigma=matrix(rep(0,(i+1)*(i+1)), nrow=i+1, ncol=i+1)
        diag(init_sigma)=1
        
        start = Sys.time()
        da = mcmcNorm(obj=cbind(y[,i], X_mis[,2:c(i+1)]), starting.values=list(beta=init_beta, sigma=init_sigma), 
                      iter=1*m, impute.every=1, seeds = c(j,j))
        tempDA = sapply(da$imp.list, function(x){coef(lm(x[,1]~x[,2:(i+1)]))[2]})
        tempse = sapply(da$imp.list, function(x){summary(lm(x[,1]~x[,2:(i+1)]))$coefficients[2,2]})
        betaDA[j] = mean(tempDA)
        iter_biasDA[j] = betaDA[j] - beta[2]
        iter_rmse[j] = sqrt(mean((betaDA[j]-beta[2])^2))
        seDA[j] = sqrt(mean(tempse^2)+var(tempDA)*(1+1/m))
        ciDA[j] = ifelse(betaDA[j]-1.96*seDA[j] <= beta[2] & betaDA[j]+1.96*seDA[j] >= beta[2],1,0)
        end = Sys.time() 
        
        time[j] = as.numeric(end-start)
      } else {
        init_beta=matrix(rep(0,i+1), nrow=1, ncol=i+1)
        init_sigma=matrix(rep(0,(i+1)*(i+1)), nrow=i+1, ncol=i+1)
        diag(init_sigma)=1
        
        start = Sys.time()
        da = mcmcNorm(obj=cbind(y[,i], X_mis[,2:c(i+1)]), starting.values=list(beta=init_beta, sigma=init_sigma), 
                      iter=num_iter[i]*2*m, impute.every=t, seeds = c(j,j))
        tempDA = sapply(da$imp.list, function(x){coef(lm(x[,1]~x[,2:(i+1)]))[2]})
        tempse = sapply(da$imp.list, function(x){summary(lm(x[,1]~x[,2:(i+1)]))$coefficients[2,2]})
        betaDA[j] = mean(tempDA)
        iter_biasDA[j] = betaDA[j] - beta[2]
        iter_rmse[j] = sqrt(mean((betaDA[j]-beta[2])^2))
        seDA[j] = sqrt(mean(tempse^2)+var(tempDA)*(1+1/m))
        ciDA[j] = ifelse(betaDA[j]-1.96*seDA[j] <= beta[2] & betaDA[j]+1.96*seDA[j] >= beta[2],1,0)
        end = Sys.time() 
        
        time[j] = as.numeric(end-start)
      }
    }
    comptime[i] = mean(time)
    biasDA[i] = mean(iter_biasDA) 
    rmseDA[i] = mean(iter_rmse)
    cicovDA[i] = mean(ciDA)
    cilenDA[i] = 2*1.96*mean(seDA)
    print(i)
  }
  data.frame(biasDA = biasDA, rmseDA = rmseDA, cicovDA = cicovDA, cilenDA = cilenDA, comptime = comptime)
}

da1 = DA(type = 1, mc=10)
da2 = DA(type = 2, mc=50)


# FCS

FCS = function(m = 10, t = t, mc = 1000, type = 1){
  set.seed(1)
  r = matrix(runif(n = 9^2,min = -1,max=1), ncol=9)
  Cor = cov2cor(r%*%t(r))
  
  seFCS = ciFCS = time = comptime = betaFCS = tempFCS = rmseFCS = tempse = cicovFCS = cilenFCS = temprmse = biasFCS = tempbias = NULL
  
  for(i in 1:9){
    print(i)
    for (n in 1:mc){
      set.seed(n)      
      aaa= mvrnorm(n = 1000, mu = rep(0,9), Sigma = Cor) # 900*6 matrix
      X = matrix(c(rep(1,1000),aaa),ncol=10)
      
      beta = runif(n = 10, min = -2, max = 2)
      sig = runif(n = 1, min = 0.5, max = 2)
      eps= rnorm(n=1000 , 0, sig)
      y =matrix(rep(0,1000*9),ncol=9)
      for(kk in 2:10){
        y[,kk-1] = X[,1:kk]%*%beta[1:kk]+eps
      }
      
      u = matrix(runif(9000,0,1),ncol=9)
      X_mis = X
      
      for(k in 2:10){
        X_mis[which(u[,k-1]<0.5 & y[,k-1]<quantile(y[,k-1],0.5) | u[,k-1]>0.9 & y[,k-1]>quantile(y[,k-1],0.5)),k] = NA
      }
      
      if(type == 1){
        ##FCS1 :with no iteration
        start = Sys.time()
        
        fcs <- cbind(y[,i],X_mis[,2:(i+1)])%>%
          mice( m =m, maxit = 1, seed = 123,print = FALSE) %>%
          mice::complete("all") 
        time[n] = as.numeric(end-start)
        end = Sys.time() 
      } else {
        ##FCS2 :with 20 iteration
        start = Sys.time()
        #num_iter 수정
        num_iter = read.csv('~/seojin/EMB_iter.csv')[-1]$x
        fcs <- cbind(y[,i],X_mis[,2:(i+1)])%>%
          mice( m = m, maxit = num_iter[i]*2, seed = 123,print = FALSE) %>%
          mice::complete("all") }
      
      end = Sys.time() 
      
      tempFCS = sapply(fcs, function(x){coef(lm(x[,1]~as.matrix(x[,2:(i+1)])))[2]})
      tempse = sapply(fcs, function(x){summary(lm(x[,1]~as.matrix(x[,2:(i+1)])))$coefficients[2,2]})
      betaFCS[n] = mean(tempFCS)
      seFCS[n] = sqrt(mean(tempse^2)+var(tempFCS)*(1+1/m))
      ciFCS[n] = ifelse(betaFCS[n]-1.96*seFCS[n] <= beta[2] & betaFCS[n]+1.96*seFCS[n] >= beta[2],1,0)
      time[n] = as.numeric(end-start)
      
      tempbias[n] = betaFCS[n] - beta[2]
      temprmse[n] = sqrt(mean((tempFCS-beta[2])^2))
      
      if (n%%100==0) {cat('n=',n,'\n')}
    }
    comptime[i] = mean(time)
    biasFCS[i] = mean(tempbias)
    rmseFCS[i] = mean(temprmse)
    cicovFCS[i] = mean(ciFCS)
    cilenFCS[i] = 2*1.96*mean(seFCS)
    
  }
  return (data.frame(biasFCS = biasFCS, rmseFCS = rmseFCS, cicovFCS = cicovFCS, cilenFCS = cilenFCS, comptime = comptime))
}
fcs1=FCS(m = 10, t = t, mc = 1000, type = 1)
fcs2=FCS(m = 10, t = t, mc = 1000, type = 2)


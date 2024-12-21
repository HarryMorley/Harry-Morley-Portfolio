lib <- c('fBasics','data.table', 'xts','zoo', 'quantmod', 'PerformanceAnalytics', 
         'PortfolioAnalytics', 'GenSA', 'timeSeries','fPortfolio', 'quantmod', 
         'plyr', 'PortfolioAnalytics', 'ROI', 'GenSA','DEoptim', 'xtable', 
         "ROI.plugin.glpk", 'ROI.plugin.quadprog','portfolioBacktest','ggplot2')
loading.lib <- lapply(lib, require, character.only = TRUE)

library(xts); library(mvtnorm)
library(quantmod)
library(PerformanceAnalytics)
library(CVXR)
dataset <- readRDS("dataset.rds")
## # upload Fama-French factors

mydata <- read.csv("F-F_Research_Data_Factors_daily.CSV",sep=",", skip = 4)
mydata <- mydata[-nrow(mydata), ]  # remove last row
fama_lib <- xts(x = mydata[, c(2,3,4)], order.by = as.Date(paste(mydata[, 1]), "%Y%m%d"))
str(fama_lib)

SP500_index <- readRDS("SP500.rds")

begin_date <- "2018-01-01"   #"2015-01-01"
end_date <- "2019-12-31"   #"2019-12-20" #"2017-12-31"  #
period<-paste(begin_date,"/",end_date,sep="")

stockPrices<-dataset$adjusted[period]     #[1:1250]
stockPrices <- stockPrices[,1:9]
tclass(stockPrices) <- "Date"
X <- diff(log(stockPrices), na.pad = FALSE)
N <- ncol(X)  # number of stocks
T <- nrow(X)  # number of days

mu <- colMeans(X)
Sigma <- cov(X)


F_FamaFrench <- fama_lib[index(X)]/100

f_SP500 <- diff(log(SP500_index), na.pad = FALSE)


BBrMkt <- na.omit(as.xts(rowAvgs(dataset$BBr),order.by=index(dataset$BBr)))/100


##create the PNlog Market index
PNlogMkt <- na.approx(as.xts(rowAvgs(dataset$PNlog),order.by=index(dataset$PNlog)))

SentIndx <- PNlogMkt[index(X)] 


FFS <- na.omit(merge(fama_lib[index(X)],SentIndx))

##Sp500 + Sent index
SPS <- na.omit(merge(f_SP500[index(X)],SentIndx))

# split data into training and test data
T_trn <- round(0.5*T)
X_trn <- X[1:T_trn, ]
X_tst <- X[(T_trn+1):T, ]
F_FamaFrench_trn <- F_FamaFrench[1:T_trn, ]
F_FamaFrench_tst <- F_FamaFrench[(T_trn+1):T, ]
f_SP500_trn <- f_SP500[1:T_trn, ]
f_SP500_tst <- f_SP500[(T_trn+1):T, ]
SentIndx_trn <-SentIndx[1:T_trn,]
SentIndx_tst <-SentIndx[(T_trn+1):T,]
FFS_trn <-FFS[1:T_trn,]
FFS_tst <-FFS[(T_trn+1):T,]
SPS_trn <-SPS[1:T_trn,]
SPS_tst <-SPS[(T_trn+1):T,]

# Fama French 3 factor
F_ <- cbind(ones = 1, F_FamaFrench_trn)
Gamma <- t(solve(t(F_) %*% F_, t(F_) %*% X_trn))
colnames(Gamma) <- c("alpha", "beta1", "beta2", "beta3")
alpha <- Gamma[, 1]
B <- Gamma[, 2:4]
E <- xts(t(t(X_trn) - Gamma %*% t(F_)), index(X_trn))
PsiFF <- (1/(T_trn-4)) * t(E) %*% E
Sigma_FamaFrench <- B %*% cov(F_FamaFrench_trn) %*% t(B) + diag(diag(PsiFF))
mu_FF <- colMeans(Sigma_FamaFrench)

## SP500 + Sentiment 2-factor model
F_ <- cbind(ones = 1, SPS_trn)
Gamma <- t(solve(t(F_) %*% F_, t(F_) %*% X_trn))
colnames(Gamma) <- c("alpha", "beta1", "beta2")
alpha <- Gamma[,1]
B <- Gamma[, 2:3]
E <- xts(t(t(X_trn) - Gamma %*% t(F_)), index(X_trn))
PsiSPS <- (1/(T_trn-4)) * t(E) %*% E
Sigma_SPS <- B %*% cov(SPS_trn) %*% t(B) + diag(diag(PsiSPS))
mu_SS <- colMeans(Sigma_SPS)

## multiple robust solutions

portfolioMaxReturnRobustEllipsoid <- function(mu_hat, S, kappa = 0.1) {
  S12 <- chol(S)  # t(S12) %*% S12 = Sigma
  w <- Variable(length(mu_hat))
  prob <- Problem(Maximize( t(w) %*% mu_hat - kappa*p_norm(S12 %*% w,p=2) ), 
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}

#FF
w_GMRP_robust_FF0.63 <- portfolioMaxReturnRobustEllipsoid(mu,Sigma,kappa = 0.63)
names(w_GMRP_robust_FF0.63) <- colnames(X)
w_all_GMRP_robust_ellipsoid_FF0.63 <- cbind(w_GMRP_robust_FF0.63)

w_GMRP_robust_FF0.9 <- portfolioMaxReturnRobustEllipsoid(mu,Sigma,kappa = 0.9)
names(w_GMRP_robust_FF0.9) <- colnames(X)
w_all_GMRP_robust_ellipsoid_FF0.9 <- cbind(w_GMRP_robust_FF0.9)

w_GMRP_robust_FF0.3 <- portfolioMaxReturnRobustEllipsoid(mu,Sigma,kappa = 0.3)
names(w_GMRP_robust_FF0.3) <- colnames(X)
w_all_GMRP_robust_ellipsoid_FF0.3 <- cbind(w_GMRP_robust_FF0.3)

set.seed(357)
for (i in 1:6) {
  #X_noisy <- rmvnorm(n = T, mean = mu, sigma = Sigma)
  #mu_noisy <- colMeans(X_noisy)
  #Sigma_noisy <- cov(X_noisy)
  w_GMRP_robust_ellipsoid_noisyFF0.9 <- portfolioMaxReturnRobustEllipsoid(mu_FF, Sigma_FamaFrench, kappa = 0.9)
  w_GMRP_robust_ellipsoid_noisyFF0.3 <- portfolioMaxReturnRobustEllipsoid(mu_FF, Sigma_FamaFrench, kappa = 0.3)
  w_GMRP_robust_ellipsoid_noisyFF0.63 <- portfolioMaxReturnRobustEllipsoid(mu_FF, Sigma_FamaFrench, kappa = 0.63)
  
  w_all_GMRP_robust_ellipsoid_FF0.63 <- cbind(w_all_GMRP_robust_ellipsoid_FF0.63, w_GMRP_robust_ellipsoid_noisyFF0.63)
  
  w_all_GMRP_robust_ellipsoid_FF0.9 <- cbind(w_all_GMRP_robust_ellipsoid_FF0.9, w_GMRP_robust_ellipsoid_noisyFF0.9)
  
  w_all_GMRP_robust_ellipsoid_FF0.3 <- cbind(w_all_GMRP_robust_ellipsoid_FF0.3, w_GMRP_robust_ellipsoid_noisyFF0.3)
}

# plot to compare the allocations
barplot(t(w_all_GMRP_robust_ellipsoid_FF0.63), col = 1:7, #legend = colnames(w_all_GMRP_robust_ellipsoid), 
        beside = TRUE, args.legend = list(bg = "white"),
        main = paste("FF 3-F Robust Global maximum return portfolio allocation, k=0.63"), xlab = "stocks", ylab = "dollars")


chart.StackedBar(t(w_all_GMRP_robust_ellipsoid_FF0.63), 
                 main = paste("FF 3-F Robust Global maximum return portfolio allocation, k=0.63"), 
                 ylab = "w", space = 0, border = NA)

barplot(t(w_all_GMRP_robust_ellipsoid_FF0.9), col = 1:7, #legend = colnames(w_all_GMRP_robust_ellipsoid), 
        beside = TRUE, args.legend = list(bg = "white"),
        main = paste("FF 3-F Robust Global maximum return portfolio allocation, k=0.9"), xlab = "stocks", ylab = "dollars")


chart.StackedBar(t(w_all_GMRP_robust_ellipsoid_FF0.9), 
                 main = paste("FF 3-F Robust Global maximum return portfolio allocation, k=0.9"), 
                 ylab = "w", space = 0, border = NA)

barplot(t(w_all_GMRP_robust_ellipsoid_FF0.3), col = 1:7, #legend = colnames(w_all_GMRP_robust_ellipsoid), 
        beside = TRUE, args.legend = list(bg = "white"),
        main = paste("FF 3-F Robust Global maximum return portfolio allocation, k=0.3"), xlab = "stocks", ylab = "dollars")


chart.StackedBar(t(w_all_GMRP_robust_ellipsoid_FF0.3), 
                 main = paste("FF 3-F Robust Global maximum return portfolio allocation, k=0.3"), 
                 ylab = "w", space = 0, border = NA)


################ SS

w_GMRP_robust_SS0.63 <- portfolioMaxReturnRobustEllipsoid(mu,Sigma,kappa = 0.63)
names(w_GMRP_robust_SS0.63) <- colnames(X)
w_all_GMRP_robust_ellipsoid_SS0.63 <- cbind(w_GMRP_robust_SS0.63)

w_GMRP_robust_SS0.9 <- portfolioMaxReturnRobustEllipsoid(mu,Sigma,kappa = 0.9)
names(w_GMRP_robust_SS0.9) <- colnames(X)
w_all_GMRP_robust_ellipsoid_SS0.9 <- cbind(w_GMRP_robust_SS0.9)

w_GMRP_robust_SS0.3 <- portfolioMaxReturnRobustEllipsoid(mu,Sigma,kappa = 0.3)
names(w_GMRP_robust_SS0.3) <- colnames(X)
w_all_GMRP_robust_ellipsoid_SS0.3 <- cbind(w_GMRP_robust_SS0.3)

set.seed(357)
for (i in 1:6) {
  #X_noisy <- rmvnorm(n = T, mean = mu, sigma = Sigma)
  #mu_noisy <- colMeans(X_noisy)
  #Sigma_noisy <- cov(X_noisy)
  
  w_GMRP_robust_ellipsoid_noisySS0.63 <- portfolioMaxReturnRobustEllipsoid(mu_SS, Sigma_SPS, kappa = 0.63)
  w_all_GMRP_robust_ellipsoid_SS0.63 <- cbind(w_all_GMRP_robust_ellipsoid_SS0.63, w_GMRP_robust_ellipsoid_noisySS0.63)
  
  w_GMRP_robust_ellipsoid_noisySS0.9 <- portfolioMaxReturnRobustEllipsoid(mu_SS, Sigma_SPS, kappa = 0.9)
  w_all_GMRP_robust_ellipsoid_SS0.9 <- cbind(w_all_GMRP_robust_ellipsoid_SS0.9, w_GMRP_robust_ellipsoid_noisySS0.9)
  
  w_GMRP_robust_ellipsoid_noisySS0.3 <- portfolioMaxReturnRobustEllipsoid(mu_SS, Sigma_SPS, kappa = 0.3)
  w_all_GMRP_robust_ellipsoid_SS0.3 <- cbind(w_all_GMRP_robust_ellipsoid_SS0.3, w_GMRP_robust_ellipsoid_noisySS0.3)
}

# plot to compare the allocations
barplot(t(w_all_GMRP_robust_ellipsoid_SS0.63), col = 1:7, #legend = colnames(w_all_GMRP_robust_ellipsoid), 
        beside = TRUE, args.legend = list(bg = "white"),
        main = paste("Sentiment indicator Robust Global maximum return portfolio allocation, k=0.63"), xlab = "stocks", ylab = "dollars")


chart.StackedBar(t(w_all_GMRP_robust_ellipsoid_SS0.63), 
                 main = paste("Sentiment indicator Robust Global maximum return portfolio allocation, k=0.63"), 
                 ylab = "w", space = 0, border = NA)

barplot(t(w_all_GMRP_robust_ellipsoid_SS0.9), col = 1:7, #legend = colnames(w_all_GMRP_robust_ellipsoid), 
        beside = TRUE, args.legend = list(bg = "white"),
        main = paste("Sentiment indicator Robust Global maximum return portfolio allocation, k=0.9"), xlab = "stocks", ylab = "dollars")


chart.StackedBar(t(w_all_GMRP_robust_ellipsoid_SS0.9), 
                 main = paste("Sentiment indicator Robust Global maximum return portfolio allocation, k=0.9"), 
                 ylab = "w", space = 0, border = NA)

barplot(t(w_all_GMRP_robust_ellipsoid_SS0.3), col = 1:7, #legend = colnames(w_all_GMRP_robust_ellipsoid), 
        beside = TRUE, args.legend = list(bg = "white"),
        main = paste("Sentiment indicator Robust Global maximum return portfolio allocation, k=0.3"), xlab = "stocks", ylab = "dollars")


chart.StackedBar(t(w_all_GMRP_robust_ellipsoid_SS0.3), 
                 main = paste("Sentiment indicator Robust Global maximum return portfolio allocation, k=0.3"), 
                 ylab = "w", space = 0, border = NA)


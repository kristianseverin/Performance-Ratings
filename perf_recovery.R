### Parameter recovery for performance ratings

install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm, VGAM, reshape2, mnormt)

set.seed(43)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#------ create task environment -------------------
n_emp <- readRDS('/Users/kristian/Documents/Skole/7. semester/Decision Making/Exam/n_emp.Rds')
nteams <- length(n_emp)

# factor loadings (constant across all iterations)
# X <- array(NA, c(nteams, 4))
# for (n in 1:nteams) {
#   # # easy/simple sampling
#   # X[n,] <- runif(4,0,0.25) # sampling 4 factor loadings - all between 0 and 0.25
#   
#   # more complex sampling - sum of all 4 cannot exceed 1 (is that right?)
#   X_temp <- array(NA, 4)
#   Xn <- array(NA, 4)
#   for (i in 1:4) {
#     X_temp[i] <- runif(1,0,1)
#   }
#   X_sum <- sum(X_temp)
#   for (i in 1:4) {
#     Xn[i] <- ifelse(X_sum>1, (X_temp[i] - (X_sum-1)/4)-0.01, X_temp[i])
#   }
#   X[n,] <- Xn
#   
# }

# read X
X <- read_csv("/Users/kristian/Documents/Skole/7. semester/Decision Making/Exam/Factor_data.csv")
# choose sorted factor loadings
X <- X[, c("MR1", "MR2", "MR3", "MR4")]
# remove "China Market Group"
X <- X[-8,]
# make matrix
X <- as.matrix(X)

niterations <- 20 #2 #20 # 100

true_MR1 <- array(NA,c(niterations, nteams))
true_MR2 <- array(NA,c(niterations, nteams))
true_MR3 <- array(NA,c(niterations, nteams))
true_MR4 <- array(NA,c(niterations, nteams))

infer_MR1 <- array(NA,c(niterations, nteams))
infer_MR2 <- array(NA,c(niterations, nteams))
infer_MR3 <- array(NA,c(niterations, nteams))
infer_MR4 <- array(NA,c(niterations, nteams))

start_time = Sys.time()
for (i in 1:niterations) {
  
  # multivariate prior for the beta vector
  # # easy/simple sampling
  # betaX <- runif(4,0,0.25) # sampling 4 slopes - all between 0 and 0.25
  
  betaX_temp <- array(NA, 4)
  betaX <- array(NA, 4)
  for (i in 1:4) {
    betaX_temp[i] <- runif(1,0,1)
  }
  betaX_sum <- sum(betaX_temp)
  for (i in 1:4) {
    betaX[i] <- ifelse(betaX_sum>1, (betaX_temp[i] - (betaX_sum-1)/4)-0.01, betaX_temp[i])
  }

  beta0 <- runif(1,0,1)
  
  #------ National level regressions ---------------------------- 
  team_mu_probit <- array(NA,nteams)
  for (n in 1:nteams) {
    team_mu_probit[n] <- beta0 + (betaX[1]*(X[n,1])) + (betaX[2]*X[n,2]) + (betaX[3]*X[n,3]) + (betaX[4]*X[n,4])
  }
  
  #------- Reparameterising---------
  #reparamaterising beta prior for slope of preferencesin CC model
  team_mu <- array(NA,nteams)
  team_sigma <- array(NA,nteams)
  shape1_team <- array(NA,nteams)
  shape2_team <- array(NA,nteams)
  for (n in 1:nteams) { 
    # concentration (precision) of rate parameters (for beta priors)
    team_sigma[n] <- runif(1,1,100)
    #team_sigma[n] <- dunif(1,n_emp[n])  # attempts to incorporate sample size for the precision
    
    team_mu[n] <- probitlink(team_mu_probit[n], inverse = TRUE) # standardisation of rate estimate mean - probit link
    shape1_team[n] <- (team_mu[n]) * team_sigma[n]
    shape2_team[n] <- (1 - team_mu[n]) * team_sigma[n]  
  }
  
  source('/Users/kristian/Documents/Skole/7. semester/Decision Making/Exam/p_rating.R')
  perf_sims <- p_rating(nteams, n_emp, shape1_team, shape2_team)
  
  emp <- perf_sims$emp
  Y <- perf_sims$Y
  
  V <- solve(t(X)%*%X)
  
  data <- list("nteams", "n_emp", "emp","X","Y","V") 
  params <- c("beta0", "betaX") 
  
  # - run jags code
  perf.samples <- jags.parallel(data, inits=NULL, params,
                              model.file ="/Users/kristian/Documents/Skole/7. semester/Decision Making/Exam/performance.txt",
                              n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)
  
  true_MR1[i] <- betaX[1]
  true_MR2[i] <- betaX[2]
  true_MR3[i] <- betaX[3]
  true_MR4[i] <- betaX[4]
  
  # find maximum a posteriori
  Z <- perf.samples$BUGSoutput$sims.list  
  infer_MR1[i] <- MPD(Z$betaX[,1])
  infer_MR2[i] <- MPD(Z$betaX[,2])
  infer_MR3[i] <- MPD(Z$betaX[,3])
  infer_MR4[i] <- MPD(Z$betaX[,4])
  
}
end_time = Sys.time()
end_time - start_time


# plotting code courtesy of Lasse
source('216377/Module4/recov_plot.R')
pl1 <- recov_plot(true_MR1, infer_MR1, c("true MR1", "infer MR1"), 'smoothed linear fit')
pl2 <- recov_plot(true_MR2, infer_MR2, c("true MR2", "infer MR2"), 'smoothed linear fit')
pl3 <- recov_plot(true_MR3, infer_MR3, c("true MR3", "infer MR3"), 'smoothed linear fit')
pl4 <- recov_plot(true_MR4, infer_MR4, c("true MR4", "infer MR4"), 'smoothed linear fit')
ggarrange(pl1, pl2, pl3, pl4, pl5, pl6)

#traceplot(CC.samples, mfrow=c(3,2))

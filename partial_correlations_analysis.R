# set seed
set.seed(43)  # Answer to everything +1, also the worst possible hand in Meyer

# install and load packages
install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, ggplot2, glue, tidyverse, tidyr)

# setwd and read data
setwd("/Users/kristian/Documents/Skole/7. semester/Decision Making/Exam")
X <- read_csv("/Users/kristian/Documents/Skole/7. semester/Decision Making/Exam/Factor_data.csv")
Y <- read_csv("/Users/kristian/Documents/Skole/7. semester/Decision Making/Exam/all_sorted.csv")
n_emp_temp <- read_csv("/Users/kristian/Documents/Skole/7. semester/Decision Making/Exam/count_df.csv")
# remove "China Market Group"
n_emp_temp <- n_emp_temp[-8,]
# create vector 
n_emp <- as.vector(n_emp_temp$count)

# choose sorted factor loadings
X <- X[, c("MR1", "MR2", "MR3", "MR4")]
# make matrix
X <- as.matrix(X)

# choose columns of interest
Y <- Y[, c ("...1","L4", "PR")]

# make data wider by team_name
Y1 <- Y %>% 
  group_by(L4) %>% 
  mutate(rn = row_number()) %>% 
  pivot_wider(id_cols = -...1, names_from = L4, values_from = PR) %>% 
  select(-rn)
Y1 <- t(Y1)  # transpose to get 55X2263

# cut dataset in half due to it being doubble the size
Y1 <- Y1[,1:2263]

# make matrix
Y <- as.matrix(Y1)

# make vector
n_emp <- as.vector(n_emp)

# get dot product of factor loadings - necessary for JZS priors
V <- solve(t(X)%*%X)

# number of teams in the analysis
nteams <- 55



# function to run JAGS code and track parameters
partial.correlations <- function (X,Y,nteams, n_emp) {

   data <- list( "Y", "nteams","X","V", "n_emp") #data inputted into jags
   params <- c("beta0","betaX","prior_T") #parameters we'll track in jags
   
# - run jags code
   win.samples <- jags.parallel(data, inits=NULL, params,
                       model.file ="performance.txt",
                       n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster = 4)
   
   return(list(win.samples))
   
}

# run function
partial_correlations <- partial.correlations(X, Y, nteams, n_emp)

plot(density(partial_correlations$BUGSoutput$sims.list$beta0))






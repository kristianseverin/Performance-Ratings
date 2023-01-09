#install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library(TruncatedDistributions)

p_rating <- function(nteams, n_emp, shape1_team, shape2_team) {
  
  # arrays to populate for simulation
  emp <- array(NA,c(nteams,max(n_emp)))
  Y <- array(NA,c(nteams,max(n_emp)))
  
  #-------------------------------------------------------------------
  #-------------------  Individual level model -----------------------
  #-------------------------------------------------------------------
  
  #---- Group level factors ----------------------- 
  for (n in 1:nteams) {
    for (e in 1:n_emp[n]) {
      emp[n,e] <- rtbeta(1,shape1_team[n],shape2_team[n], a=0.001, b=0.999)
      Y[n,e] <- rbinom(1,1,emp[n,e]) 
    }         
  }
  
  result <- list(emp=emp,
                 Y=Y)
  
  return(result)
  
}

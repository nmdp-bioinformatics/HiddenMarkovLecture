## Install and load pacakges
if (!require("depmixS4", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("depmixS4", dependencies = TRUE)
  library("depmixS4", verbose=F, warn.conflicts =F)
}

## A simple example
data("speed")
set.seed(1)
# Speciy model parameters 
# response - responsive variable
# nstates - number of states
# trstart - strat transition matrix 
# instart - prior probablities (not specified in this example)
# respstart - parameters of the response models (not specified in this example)
mod <- depmix(response = rt ~ 1, data = speed, nstates = 2, trstart = runif(4))
# modle fitting
# Automatic generation of starting parameters using EM algorithm 
# in this case, random start values should not be generated 
fm <- fit(mod, emc=em.control(rand=FALSE))

summary(fm)

# Covariates on transition probabilities:
set.seed(1)
mod2 <- depmix(rt ~ 1, data = speed, nstates = 2, family = gaussian(), 
              transition = ~ scale(Pacc), instart = runif(2))
fm2 <- fit(mod2, verbose = FALSE, emc=em.control(rand=FALSE))

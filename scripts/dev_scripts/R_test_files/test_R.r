

library(statnet)
library(ergm)
library(sna)
library(coda)
library(network)
sessionInfo()
set.seed(0)

Nsample = 10
N = 10
T  =20
net <- network.initialize(N)
sampledMat_T_R =    array(0, dim=c(N,N,T,Nsample))

changeStats_T_R = list()
stats_T_R = list()
    for(t in 1:T){
        changeStats_t_R = list()
        stats_t_R = list()
        print(t)
        for(n in 1:Nsample){
            print(parDgpT[,t])
             net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parDgpT[,t],control = control.simulate.formula(MCMC.burnin = 100000))
             sampledMat_T_R[,,t,n] <- as.matrix.network( net)
             print(c(t,n))
             chStat_t <- ergmMPLE(formula_ergm)
             changeStats_t_R[[n]] <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
             stats_t_R[[n]] <- summary(formula_ergm)

             }
              changeStats_T_R[[t]] <-changeStats_t_R
              stats_T_R[[t]] <- stats_t_R
         }

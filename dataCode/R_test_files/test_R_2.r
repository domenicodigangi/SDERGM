library(statnet)
library(ergm)
library(sna)
library(coda)
library(network)
#library(VCERGM)
sessionInfo()

set.seed(0)

Ntrems = 2
#create an empty network, the formula defining ergm, sample the ensemble and store in R

load('~/Dropbox/Dynamic_Networks/data/congress_covoting_US/Rollcall_VCERGM.RData')
obsMat_T_R =     list()
estParSS_T_R =    list()
changeStats_T_R = list()
stats_T_R = list()
networks = Rollcall$networks # Networks
attr = Rollcall$attr # Political affiliation
T<-length(Rollcall$networks)
net <- network.initialize(5)

for(t in 1:T){
estParSS_t_R =    list()
changeStats_t_R = list()

net <-network(networks[[t]] , directed = FALSE)
#set.network.attribute(net, 'attr1', attr[[t]])
net %v% 'attr1' <-attr[[t]]
formula_ergm = net ~  edges +  gwesp(alpha = 0.5,fixed = TRUE)  + nodematch('attr1')#  #
#formula_ergm = net ~  triangle + kstar(2)  +   nodematch('attr1')# # edges +  triangle +  kstar(2) +   nodematch('attr1') #  gwesp(decay = 0.5,fixed = TRUE) + gwnsp(decay = 0.5,fixed = TRUE)  # triangle + kstar(2)  +   nodematch('attr1')# + nodematch('attr1') #  gwesp(decay = 0.5,fixed = TRUE) gwdegree(decay = 0.5,fixed = TRUE)# edges  +   kstar(2) #  + nodematch('attr1')  #
obsMat_T_R[[t]] <- as.matrix(net)
tmp <- ergm(formula_ergm,control = control.ergm(init=c(-6,1.5,2.5)))
estParSS_T_R[[t]] <- tmp[[1]]
print(t)
print(estParSS_T_R[[t]])
chStat_t <- ergmMPLE(formula_ergm)
changeStats_T_R[[t]] <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
stats_T_R[[t]] <- summary(formula_ergm)
}




















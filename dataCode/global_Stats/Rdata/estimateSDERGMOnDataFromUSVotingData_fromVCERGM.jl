
# sample sequences of ergms with different parameters' values from R package ergm
# and test the PseudoLikelihoodScoreDrivenERGM filter
using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase,CSV, RCall
using PyCall; pygui(); using PyPlot
# load the required packages in R
R"library(statnet)
    library(ergm)
    library(sna)
    library(coda)
    library(network)
    #library(VCERGM)
    data(Rollcall)
    sessionInfo()"

R"set.seed(0)"

Ntrems = 2
 #create an empty network, the formula defining ergm, sample the ensemble and store in R
 R"
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
       formula_ergm = net ~  edges +  gwesp(alpha = 0.5,fixed = TRUE) # + nodematch('attr1')#  #
       #formula_ergm = net ~  edges +  triangle +  nodematch('attr1')# # edges +  triangle +  kstar(2) +   nodematch('attr1') #  gwesp(decay = 0.5,fixed = TRUE) + gwnsp(decay = 0.5,fixed = TRUE)  # triangle + kstar(2)  +   nodematch('attr1')# + nodematch('attr1') #  gwesp(decay = 0.5,fixed = TRUE) gwdegree(decay = 0.5,fixed = TRUE)# edges  +   kstar(2) #  + nodematch('attr1')  #
      obsMat_T_R[[t]] <- as.matrix(net)
     tmp <- ergm(formula_ergm)
     estParSS_T_R[[t]] <- tmp[[1]]
     print(t)
     print(estParSS_T_R[[t]])
     chStat_t <- ergmMPLE(formula_ergm)
     changeStats_T_R[[t]] <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
     stats_T_R[[t]] <- summary(formula_ergm)
      }
  "
# LE MIE STIME NON FUNZIONANO SUI DATI USANDO COME STATISTICA I TRIANGOLI E LE K star, o cmq sono un po' noisy
# LE STIME DA VCERGM NON CONVERGONO INVECE CON STATISTICHE WELL BEHAVED.
#
 # import in julia
T = @rget(T) -1 #
 changeStats_T = @rget changeStats_T_R; tmp = Array{Array{Float64,2}}(T); for t=1:T tmp[t] =  (changeStats_T[t]);end;changeStats_T = tmp
 Nterms = length(changeStats_T[1][1,:])-2
 obsMat_T =  @rget(obsMat_T_R);tmp = [BitArray(obsMat_T[i]) for i=1:T]
 estParSS_T = @rget(estParSS_T_R);tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = estParSS_T[t]; end ; estParSS_T = tmp
 stats_T = @rget(stats_T_R); tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = stats_T[t];end;stats_T = tmp


# save

@save("/home/Domenico/Dropbox/Dynamic_Networks/data/congress_covoting_US/juliaEstimates.jld")

# #SE SIMULO E STIMO PARAMETRI COSTANTI UGUALI A QUELLI CHE STIMO SUI DATI COSA OTTENGO?
# staticPar,convFlag = estimate(model;UM = startUM,indTvPar = falses(Nterms))
# tmp = zeros(Nterms); for i=1:Nterms tmp[i] = staticPar[i][1]; end;parMatDgp_T = repmat(tmp,1,T)
# Nsample = 30
# @rput(parMatDgp_T)
# @rput(Nsample)
# R"
#
#  # sampledMat_T_R =    array(0, dim=c(N,N,T,Nsample))
#   estParSS_T_R =    list()
#   changeStats_T_R = list()
#   stats_T_R = list()
#    for(t in 1:T){
#        net <-network(networks[[t]] , directed = FALSE)
#        #set.network.attribute(net, 'attr1', attr[[t]])
#        net %v% 'attr1' <-attr[[t]]
#        formula_ergm = net ~  triangle +  kstar(2) +   nodematch('attr1')
#
#        estParSS_t_R =    list()
#        changeStats_t_R = list()
#        stats_t_R = list()
#        print(t)
#        for(n in 1:Nsample){
#            print(parMatDgp_T[,t])
#             net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parMatDgp_T[,t],control = control.simulate.formula(MCMC.burnin = 100000))
#             #sampledMat_T_R[,,t,n] <- as.matrix.network( net)
#
#             tmp <- ergm(formula_ergm,estimate = 'MPLE')
#             estParSS_t_R[[n]] <- tmp[[1]]
#             print(c(t,n))
#
#             print(12)
#             #print(estParSS_t_R[[n]])
#             chStat_t <- ergmMPLE(formula_ergm)
#             changeStats_t_R[[n]] <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
#             stats_t_R[[n]] <- summary(formula_ergm)
#
#             }
#              estParSS_T_R[[t]] <- estParSS_t_R
#              changeStats_T_R[[t]] <-changeStats_t_R
#              stats_T_R[[t]] <- stats_t_R
#         }"
#
#
#
#
# # import in julia
# sampledMat_T = BitArray( @rget(sampledMat_T_R))
# estParSS_T = @rget(estParSS_T_R);#tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = estParSS_T[t]; end ; estParSS_T = tmp
# changeStats_T = @rget changeStats_T_R;# tmp = Array{Array{Float64,2}}(T); for t=1:T tmp[t] =  changeStats_T[t];end;changeStats_T = tmp
# stats_T = @rget(stats_T_R); #tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = stats_T[t];end;stats_T = tmp
#
#
#
# estParSS_TPlot = Array{Matrix{Float64}}(Nterms)
# for i=1:Nterms
#     tmp = zeros(Nsample,T)
#     for t = 1:T
#         for n=1:Nsample
#             tmp[n,t] = estParSS_T[t][n][i]
#         end
#     end
#     estParSS_TPlot[i] = tmp
# end
#
# close()
#  for parInd = 1:Nterms
#     subplot(Nterms,1,parInd); #plot(1:T,gasFiltPar[parInd,:],"r")
#     for n=1:Nsample
#               plot(1:T,estParSS_TPlot[parInd][n,:],".g")
#     end
# end
#
#
#
#
#
# Nterms = 3
#
#
#
#
#
#
#
#
#
#
#













#

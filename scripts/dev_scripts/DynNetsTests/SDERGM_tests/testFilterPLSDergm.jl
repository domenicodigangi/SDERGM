
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
    sessionInfo()"

R"set.seed(0)"


Nsample = 100
T = 50
N = 50
Nterms = 2
# the matrix of parameters values for each t

parDgpT = zeros(Nterms,T)
 parDgpT[1,:] = randSteps(-2.5,-3,2,T)# -3# randSteps(0.05,0.5,2,T) #1.5#.00000000000000001
 parDgpT[2,:] = randSteps(0.05,0.5,2,T)#
 #parDgpT[2,:] = 0.25
 @rput T; @rput parDgpT;@rput N
 #create an empty network, the formula defining ergm, sample the ensemble and store in R
 R"
 net <- network.initialize(N)
 formula_ergm = net ~ edges + gwesp(decay = 0.25,fixed = TRUE)#gwnsp(decay = 0.25,fixed = TRUE)
   sampledMat_T_R =    array(0, dim=c(N,N,T))
   estParSS_T_R =    list()
   changeStats_T_R = list()
   stats_T_R = list()
    for(t in 1:T){
        print(t)
        print(parDgpT[,t])
         net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parDgpT[,t],control = control.simulate.formula(MCMC.burnin = 500000))
         sampledMat_T_R[,,t] <- as.matrix.network( net)
         tmp <- ergm(formula_ergm,estimate = 'MPLE')#)#
         estParSS_T_R[[t]] <- tmp[[1]]
         print(estParSS_T_R[[t]])
         chStat_t <- ergmMPLE(formula_ergm)
         changeStats_T_R[[t]] <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
         stats_T_R[[t]] <- summary(formula_ergm)
         }"
 # import in julia
 sampledMat_T = BitArray( @rget(sampledMat_T_R))
 estParSS_T = @rget(estParSS_T_R);tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = estParSS_T[t]; end ; estParSS_T = tmp
 changeStats_T = @rget changeStats_T_R; tmp = Array{Array{Float64,2}}(T); for t=1:T tmp[t] =  changeStats_T[t];end;changeStats_T = tmp
 stats_T = @rget(stats_T_R); tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = stats_T[t];end;stats_T = tmp
 close()
figure()
 subplot(2,2,1);  plot(1:T,sumSq(sumSq(sampledMat_T,1),1)./(N^2-N))
 subplot(2,2,2);  plot(1:T,stats_T[1,:]  )
 subplot(2,2,3); plot(1:T,parDgpT[1,:],"k",1:T,estParSS_T[1,:],"c")
 subplot(2,2,4); plot(1:T,parDgpT[2,:],"k",1:T,estParSS_T[2,:],"c")
 #cor(stats_T')

 estParStatic,convFlag,UM , ftot_0= estimate(GasNetModelDirBinGlobalPseudo(changeStats_T,testGasPar,trues(2),"");UM = [0;0],indTvPar = falses(Nterms))
 constParVec = zeros(sum(!indTvPar)); for i=1:sum(.!indTvPar) constParVec[i] = estParStatic[.!indTvPar][i][1]; end
 grad_T,s_T,I_T = gasScoreSeries(GasNetModelDirBinGlobalPseudo(changeStats_T,testGasPar,trues(2),""),constParVec)


 # run one auxiliary regression for each parameters
 using GLM
 for i=1:2
     g_T_t = Float64.(grad_T[i,2:end])
     s_T_tm1 = Float64.(s_T[i,1:end-1])
     g_T_tm1 = Float64.(grad_T[i,1:end-1])
     g_delta_t = Float64.(grad_T[(1:Nterms).!=i,2:end])
     I_T_t = Float64.(I_T[i,2:end])
     #ols = lm(@formula(x1 ~ x2 + x3 + x4 ),data)
     autocor(g_T_t,[1])
     y = ones(T-1)
     ols = lm([  g_delta_t' g_T_t  s_T_tm1.*g_T_t],y)
     yHat = predict(ols)
     yMean =1
     RSS = sum((yHat- yMean).^2)
     ESS = T - RSS
     pval =  ccdf(Chisq(1), ESS)
     scoreAcorr = autocor(Float64.(grad_T[i,:]),[1])[1]
     @show(i,pval,scoreAcorr)
 end

 figure()
 plot(s_T' -mean(s_T)); grid()

## test SDERGM
#  um = 0.5
#   B,A = 0.9, 0.02
#   testGasPar = [[-1.45],[um*(1-B),B ,A]]
#   indTvPar = trues(2);indTvPar[1] = false
#   testModel = DynNets.GasNetModelDirBinGlobalPseudo(changeStats_T,testGasPar,trues(1),"")
# #
#   gasFiltPar , pseudolike = DynNets.score_driven_filter(testModel,testGasPar[2],indTvPar;
#     obsT = changeStats_T,vConstPar = testGasPar[1])
#     close()
#     subplot(2,2,1);  plot(1:T,sumSq(sumSq(sampledMat_T,1),1)./(N^2-N))
#     subplot(2,2,2);  plot(1:T,stats_T[1,:]  )
#     subplot(2,2,3); plot(1:T,gasFiltPar[1,:],"k",1:T,estParSS_T[1,:],"c")
#     subplot(2,2,4);  plot(1:T,parDgpT[2,:],"k",1:T,estParSS_T[2,:],"c")
#     subplot(2,2,4); plot(1:T,gasFiltPar[2,:],"k",1:T,estParSS_T[2,:],"c")
#  pseudolike
# #
#
# Npoints = 100
#  testPseudoLike = zeros(Npoints)
#  indTvPar = falses(2);
#  parVals=linspace(-3,3,Npoints)
#  for n=1:Npoints
#      testGasPar = [[parVals[n]],[0.5]]
#     ~,testPseudoLike[n]  = score_driven_filter(testModel,zeros(1),indTvPar;
#        obsT = changeStats_T,vConstPar = [testGasPar[1]; testGasPar[2]])
#  end
#  plot(parVals,testPseudoLike)
um = 0.5
 B,A = 0.9, 0.02
 testGasPar = [[-1.45],[um*(1-B),B ,A]]
 model = GasNetModelDirBinGlobalPseudo(changeStats_T,testGasPar,trues(2),"")
 indTvPar = BitArray([true,true])
 estPar,convFlag,UM , ftot_0 = estimate(GasNetModelDirBinGlobalPseudo(changeStats_T,testGasPar,trues(2),"");UM = [0;0],indTvPar = indTvPar,
                                                                        indTargPar = indTvPar)#  falses(indTvPar))##
 gasParVec = zeros(sum(indTvPar)*3); for i=0:(sum(indTvPar)-1) gasParVec[1+3i : 3(i+1)] = estPar[indTvPar][i+1]; end
 constParVec = zeros(sum(!indTvPar)); for i=1:sum(.!indTvPar) constParVec[i] = estPar[.!indTvPar][i][1]; end


#gasParVec[1] = -0.101
 #gasParVec[3] = 0.01
gasFiltPar , pseudolike = score_driven_filter(GasNetModelDirBinGlobalPseudo(changeStats_T,testGasPar,trues(2),""),gasParVec,indTvPar;
                          vConstPar = constParVec,ftot_0= ftot_0)

 close()
    subplot(2,2,1);  plot(1:T,sumSq(sumSq(sampledMat_T,1),1)./(N^2-N))
    subplot(2,2,2);  plot(1:T,stats_T[1,:]  )
    subplot(2,2,3); plot(1:T,parDgpT[1,:],"k",1:T,gasFiltPar[1,:],"--k",1:T,estParSS_T[1,:],"c")
    subplot(2,2,4); plot(1:T,parDgpT[2,:],"k",1:T,gasFiltPar[2,:],"--k",1:T,estParSS_T[2,:],"c")


###
# close()
#  plt[:hist](estParSS_T[2,1:25])
#
#

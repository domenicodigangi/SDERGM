
# sample sequences of ergms with different parameters' values from R package ergm
# and test the PseudoLikelihoodScoreDrivenERGM filter
using Utilities,AReg,StaticNets,DynNets , JLD,MLBase,StatsBase,CSV, RCall
 using PyCall; pygui(); using PyPlot
 # load the required packages in R
 R"library(statnet)
    library(ergm)
    library(sna)
    library(coda)
    library(network)
    sessionInfo()"

 R"set.seed(0)"


 parDgpT =  dgp_paper_SDERGM(Nsample = 100,dgpType = "sin", T = 50, N = 50,
                              Nterms = 2, Nsteps1 = 1 ,Nsteps2 = 0)


# For each t, sample the ERGM with parameters corresponding to the DGP at time t
# and run a single snapeshot estimate
@rput T; @rput parDgpT;@rput N;@rput Nsample
 #create an empty network, the formula defining ergm, sample the ensemble and
 # store the sufficient statistics and change statistics in R
 R"
 net <- network.initialize(N)
  #formula_ergm = net ~ edges + gwesp(decay = 0.25,fixed = TRUE)
  formula_ergm = net ~ edges +  gwesp(decay = 0.25,fixed = TRUE)
   sampledMat_T_R =    array(0, dim=c(N,N,T,Nsample))
   estParSS_T_R =    list()
   changeStats_T_R = list()
   stats_T_R = list()
    for(t in 1:T){
        estParSS_t_R =    list()
        changeStats_t_R = list()
        stats_t_R = list()
        print(t)
        for(n in 1:Nsample){
            print(parDgpT[,t])
             net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parDgpT[,t],control = control.simulate.formula(MCMC.burnin = 100000))
             sampledMat_T_R[,,t,n] <- as.matrix.network( net)
             tmp <- ergm(formula_ergm)#,estimate = 'MPLE')#)#
             estParSS_t_R[[n]] <- tmp[[1]]
             print(c(t,n))

             print(estParSS_t_R[[n]])
             chStat_t <- ergmMPLE(formula_ergm)
             changeStats_t_R[[n]] <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
             stats_t_R[[n]] <- summary(formula_ergm)

             }
              estParSS_T_R[[t]] <- estParSS_t_R
              changeStats_T_R[[t]] <-changeStats_t_R
              stats_T_R[[t]] <- stats_t_R
         }"


 # import sampled networks in julia
 sampledMat_T = BitArray( @rget(sampledMat_T_R))
 estParSS_T = @rget(estParSS_T_R);#tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = estParSS_T[t]; end ; estParSS_T = tmp
 changeStats_T = @rget changeStats_T_R;# tmp = Array{Array{Float64,2}}(T); for t=1:T tmp[t] =  changeStats_T[t];end;changeStats_T = tmp
 stats_T = @rget(stats_T_R); #tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = stats_T[t];end;stats_T = tmp
 estParSS_T = Utilities.collapseArr3(estParSS_T)

# Save all Sampled Data
save_fold = "./data/estimatesTest/sdergmTest/R_MCMC_estimates/"
 @save(save_fold*"test_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE.jld",
              stats_T, changeStats_T,estParSS_T,sampledMat_T ,parDgpT,Nsample,T,N)


using Utilities,AReg,StaticNets,DynNets , JLD,MLBase,StatsBase,CSV, RCall
 using PyCall; pygui(); using PyPlot
 using JLD2,Utilities, GLM
 ## load R MCMC simulation and estimates and estimate sdergmTest

 if false
 Nsample = 100
 dgpType = "sin"
 T = 50
 N = 50
 Nterms = 2
 Nsteps1 ,Nsteps2 = 0,1
 load_fold = "./data/estimatesTest/sdergmTest/R_MCMC_estimates/"
 @load(load_fold*"test_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE.jld",
             stats_T, changeStats_T,estParSS_T,sampledMat_T ,parDgpT,Nsample)
 end

 onlyTest = false
 targetAllTv = true
 tvParFromTest = false; pValTh = 0.05
 gasParEst = fill(fill(Float64[],2), Nsample)
 startPointEst = fill(fill(0.0,2), Nsample)
 staticEst = fill(fill(0.0,2), Nsample)
 filtPar_T_Nsample = zeros(Nterms,T,Nsample)
 pVals_Nsample = zeros(Nterms,Nsample)
 scoreAutoc_Nsample = zeros(Nterms,Nsample)
 convFlag = falses(Nsample)
 tmp = Array{Array{Float64,2},2}(T,Nsample); for t=1:T,s=1:Nsample tmp[t,s] =  changeStats_T[t][s];end;#changeStats_T = tmp
 for n=1:Nsample
     @show(n)
      pVals_Nsample[:,n],tmpInfo,staticEst[n] =
         DynNets.pValStatic_SDERGM( DynNets.SdErgmDirBinGlobalPseudo(tmp[:,n],DynNets.fooGasPar,falses(Nterms),""))
      @show(tmpInfo)
     if tvParFromTest
         indTvPar =  pVals_Nsample[:,n] .< pValTh
     else
         indTvPar = BitArray([true,true])
     end
     if targetAllTv
         indTargPar = indTvPar
     else
          indTargPar =BitArray([false,false])
     end
    model = DynNets.SdErgmDirBinGlobalPseudo(tmp[:,n],DynNets.fooGasPar,indTvPar,"")

    if .!onlyTest
    estPar,convFlag[n],UM,startPoint = DynNets.estimate(model;UM =meanSq(estParSS_T[n,:,:],1),indTvPar = indTvPar, indTargPar = indTargPar)


    gasParEst[n] = estPar # store gas parameters
    startPointEst[n] = startPoint
    gasParVec = zeros(sum(indTvPar)*3); for i=0:(sum(indTvPar)-1) gasParVec[1+3i : 3(i+1)] = estPar[indTvPar][i+1]; end
    constParVec = zeros(sum(!indTvPar)); for i=1:sum(.!indTvPar) constParVec[i] = estPar[.!indTvPar][i][1]; end
    gasFiltPar , pseudolike = DynNets.score_driven_filter_or_dgp(model,gasParVec,indTvPar;
                              vConstPar = constParVec,ftot_0= startPoint )


    filtPar_T_Nsample[:,:,n] = gasFiltPar
    end
 end

 save_fold = "./data/estimatesTest/sdergmTest/"*
            "gas_MCMC_comparison_estimates/"
 @save(save_fold*"test_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE_target_$(targetAllTv).jld" ,
     stats_T, changeStats_T,estParSS_T,sampledMat_T ,parDgpT,
     Nsample,T,N,filtPar_T_Nsample,gasParEst,convFlag,pVals_Nsample,scoreAutoc_Nsample,staticEst)


if false

    sum(convFlag)
    figure()
     legTex = ["Constant ERGM";"SDERGM"; "Sequence of ERGM"; "DGP"]
     for n=1:Nsample
     indTvPar = BitArray([true,true])
     model = SdErgmDirBinGlobalPseudo(tmp[:,n],fooGasPar,indTvPar,"")
     estPar = gasParEst[n] # store gas parameters
     gasFiltPar  = filtPar_T_Nsample[:,:,n]
     parInd = 1
     subplot(1,2,1);plot(1:T,ones(T)*staticEst[n][parInd],"-b")
                    plot(1:T,gasFiltPar[parInd,:],"r")
                    plot(1:T,estParSS_T[n,:,parInd],".b")
                    plot(1:T,parDgpT[parInd,:],"k",linewidth=5)

     parInd = 2
     subplot(1,2,2);plot(1:T,ones(T)*staticEst[n][parInd],"-b")
                    plot(1:T,gasFiltPar[parInd,:],"r")
                    plot(1:T,estParSS_T[n,:,parInd],".b")
                    plot(1:T,parDgpT[parInd,:],"k",linewidth=5)
     end
     namePar1 = "Number of Links"
     namePar2 = "GWESP"
     subplot(1,2,1);    title(namePar1); legend(legTex)
     subplot(1,2,2);  title(namePar2 ); legend(legTex)

    #Plotta info relative ai tests
    figure()
         subplot(2,2,1); parInd = 1;   plt[:hist]((pVals_Nsample[parInd,:]), bins=logspace(minimum(log10(pVals_Nsample[1,:])),maximum(log10(pVals_Nsample[parInd,:])), 20))
         if maximum((pVals_Nsample[parInd,:])) >pValTh
              axvspan((pValTh),maximum((pVals_Nsample[parInd,:])),color = "r",alpha = 0.1);
         end
              xscale("log");
         subplot(2,2,2);  parInd = 2;   plt[:hist]((pVals_Nsample[parInd,:]), bins=logspace(minimum(log10(pVals_Nsample[1,:])),maximum(log10(pVals_Nsample[parInd,:])), 20))
         if maximum((pVals_Nsample[parInd,:])) >pValTh
              axvspan((pValTh),maximum((pVals_Nsample[parInd,:])),color = "r",alpha = 0.1);
         end
         subplot(2,2,3); plot((pVals_Nsample[1,:]),(pVals_Nsample[2,:]),".")
         axvline((pValTh));axhline((pValTh))
         xscale("log");
         xscale("log");yscale("log")
         xlabel("P-Value " * namePar1)
         ylabel("P-Value " * namePar2)
end

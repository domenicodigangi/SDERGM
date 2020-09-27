
# sample sequences of ergms with different parameters' values from R package ergm
# and test the PseudoLikelihoodScoreDrivenERGM filter
using HelperFunDom,AReg,StaNets,JLD,MLBase,StatsBase,CSV, RCall
using PyCall; pygui(); using PyPlot


using JLD,HelperFunDom, GLM
## load R MCMC simulation and estimates and estimate sdergmTest

Nsample = 100
 targetAllTv = true
 dgpType = "AR"
 T = 50
 N = 50
 Nterms = 2
 Nsteps1 ,Nsteps2 = 2,2
 load_fold = "./data/estimatesTest/sdergmTest/gas_MCMC_comparison_estimates/"
 @load(load_fold *"test_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE_target_$(targetAllTv).jld" ,
     stats_T, changeStats_T,estParSS_T, parMatDgp_T,
     Nsample,T,N,filtPar_T_Nsample,gasParEst,convFlag,pVals_Nsample,scoreAutoc_Nsample,staticEst)

pValTh = 0.05
 tmp = sum( pVals_Nsample.<pValTh,2)
 (-tmp[1] + tmp[2] + 100)/200
#compute the Average MSE
MSE_sample  = zeros(Nterms,3,Nsample)
 for n=1:Nsample
     if sum(isnan(filtPar_T_Nsample[:,:,n]))==0
    ssDiff = ssFiltPar = estParSS_T[n,:,:]' .-parMatDgp_T
    constDiff =   staticEst[n] .-parMatDgp_T
    gasDiff = filtPar_T_Nsample[:,:,n].-parMatDgp_T
    #MSE_sample[:,1,n] = mean(constDiff.^2,2)
    MSE_sample[:,2,n] = mean(ssDiff.^2,2)
    MSE_sample[:,3,n] = mean(gasDiff.^2,2)
    else
        pVals_Nsample[:,n] = 0.1
    end

 end
 avgMSE = meanSq(MSE_sample,3)
 display(avgMSE)



close()
 figure(figsize=(9,7))
 labSize = 26
 tickSize = 20
  titSize = 28
 legSize = 20
 legTex = ["SD-ERGM "; " ERGM "; "DGP"]
 legTex1 = ["Constant ERGM";"SDERGM <MSE> = $(avgMSE[2,1])"; "Sequence of ERGM = $(avgMSE[2,1])"; "DGP"]
 legTex2 = ["Constant ERGM";"SDERGM <MSE> = $(avgMSE[2,2])"; "Sequence of ERGM = $(avgMSE[2,1])"; "DGP"]
 for n=1:Nsample
 indTvPar = BitArray([true,true])
 gasFiltPar  = filtPar_T_Nsample[:,:,n]
 parInd = 1
 subplot(1,2,1);#plot(1:T,ones(T)*staticEst[n][parInd],"-b")
                plot(1:T,gasFiltPar[parInd,:],"r")
                plot(1:T,estParSS_T[n,:,parInd],".b",markersize = 2)
                plot(1:T,parMatDgp_T[parInd,:],"k",linewidth=5)

 parInd = 2
 subplot(1,2,2);#plot(1:T,ones(T)*staticEst[n][parInd],"-b")
                plot(1:T,gasFiltPar[parInd,:],"r")
                plot(1:T,estParSS_T[n,:,parInd],".b",markersize = 2)
                plot(1:T,parMatDgp_T[parInd,:],"k",linewidth=5)
 end
 namePar1 = "Number of Links"
 namePar2 = "GWESP"
 subplot(1,2,1);    title(namePar1,fontsize = titSize); legend(legTex,fontsize = legSize,framealpha = 0.95)
 ylim([-3.5,-1.8])
 xticks(fontsize = tickSize)
 yticks(fontsize = tickSize)
 subplot(1,2,2);  title(namePar2,fontsize = titSize); legend(legTex,fontsize = legSize,framealpha = 0.95)
 ylim([-0.5,1.2])

  xticks(fontsize = tickSize)
  yticks(fontsize = tickSize)
   tight_layout()

 #Plotta info relative ai tests
figure()
  pValTh = 0.05
      subplot(2,2,1); parInd = 1;   plt[:hist]((pVals_Nsample[parInd,:]), bins=logspace(minimum(log10(pVals_Nsample[1,:])),maximum(log10(pVals_Nsample[parInd,:])), 20))
      if maximum((pVals_Nsample[parInd,:])) >pValTh
           axvspan((pValTh),maximum((pVals_Nsample[parInd,:])),color = "r",alpha = 0.1);
      end
           xscale("log");
      subplot(2,2,2);  parInd = 2;   plt[:hist]((pVals_Nsample[parInd,:]), bins=logspace(minimum(log10(pVals_Nsample[2,:])),maximum(log10(pVals_Nsample[parInd,:])), 20))
      if maximum((pVals_Nsample[parInd,:])) >pValTh
           axvspan((pValTh),maximum((pVals_Nsample[parInd,:])),color = "r",alpha = 0.1);
      end
           xscale("log");
      subplot(2,2,3); plot((pVals_Nsample[1,:]),(pVals_Nsample[2,:]),".")
      axvline((pValTh));axhline((pValTh))
      xscale("log");yscale("log")
      xlabel("P-Value " * namePar1)
      ylabel("P-Value " * namePar2)
      tight_layout()

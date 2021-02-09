
# sample sequences of ergms with different parameters' values from R package ergm
# and test the PseudoLikelihoodScoreDrivenERGM filter
using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase,CSV, RCall
using PyCall; pygui(); using PyPlot


using JLD,Utilities, GLM
## load R MCMC simulation and estimates and estimate sdergmTest

Nsample = 50
 dgpType = "steps"
 T = 50
 N = 50
 Nterms = 2
 Nsteps1 ,Nsteps2 = 2,1
 load_fold = "./data/estimatesTest/sdergmTest/R_MCMC_estimates/"
 @load(load_fold*"test_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE.jld",
             stats_T, changeStats_T,estParSS_T,sampledMat_T ,parDgpT,Nsample)


onlyTest = false
 targetAllTv = true
 tvParFromTest =true; pValTh = 0.05
 gasParEst = fill(fill(Float64[],2), Nsample)
 startPointEst = fill(fill(0.0,2), Nsample)
 staticEst = fill(fill(0.0,2), Nsample)
 filtPar_T_Nsample = zeros(Nterms,T,Nsample)
 pVals_Nsample = zeros(Nterms,Nsample,2)
 scoreAutoc_Nsample = zeros(Nterms,Nsample)
 convFlag = falses(Nsample)
 tmp = Array{Array{Float64,2},2}(T,Nsample); for t=1:T,s=1:Nsample tmp[t,s] =  changeStats_T[t][s];end;#changeStats_T = tmp
 #for
     n=1#:Nsample
     @show(n)
      pVals_Nsample[:,n,1],tmpInfo,staticEst[n] =
         pValStatic_SDERGM( GasNetModelDirBinGlobalPseudo(tmp[:,n],fooGasPar,falses(Nterms),""))
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
    model = GasNetModelDirBinGlobalPseudo(tmp[:,n],fooGasPar,indTvPar,"")

    if .!onlyTest
    estPar,convFlag[n],UM,startPoint = estimate(model;UM =meanSq(estParSS_T[n,:,:],1),indTvPar = indTvPar, indTargPar = indTargPar)


    gasParEst[n] = estPar # store gas parameters
    startPointEst[n] = startPoint
    gasParVec = zeros(sum(indTvPar)*3); for i=0:(sum(indTvPar)-1) gasParVec[1+3i : 3(i+1)] = estPar[indTvPar][i+1]; end
    constParVec = zeros(sum(!indTvPar)); for i=1:sum(.!indTvPar) constParVec[i] = estPar[.!indTvPar][i][1]; end
    gasFiltPar , pseudolike = score_driven_filter(model,gasParVec,indTvPar;
                              vConstPar = constParVec,ftot_0= startPoint )


    filtPar_T_Nsample[:,:,n] = gasFiltPar
    pVals_Nsample[:,n,2],tmpInfo,~ =
       pValStatic_SDERGM( GasNetModelDirBinGlobalPseudo(tmp[:,n],fooGasPar,falses(Nterms),"");estPar = gasFiltPar)
    #end
 end
 close("all")

tmp = gasScoreSeries(model,gasFiltPar;matScal=false)
plot(tmp[1]')
autocor(Float64.(tmp[1])',[1])
pValStatic_SDERGM( GasNetModelDirBinGlobalPseudo(tmp[:,n],fooGasPar,falses(Nterms),""))

pValStatic_SDERGM( GasNetModelDirBinGlobalPseudo(tmp[:,n],fooGasPar,falses(Nterms),"");estPar = gasFiltPar)
  # Estimation
#Plotta info relative ai tests
for testInt =1:2
 figure()
 namePar1 = "Number of Links"
 namePar2 = "GWESP"
 subplot(1,2,1);    title(namePar1)
 subplot(1,2,2);  title(namePar2 );# legend(legTex)

 pVals_Nsample[isnan(pVals_Nsample)] = pValTh
  pValTh = 0.05
      subplot(2,2,1); parInd = 1;   plt[:hist]((pVals_Nsample[parInd,:,testInt]), bins=logspace(minimum(log10(pVals_Nsample[1,:,testInt])),maximum(log10(pVals_Nsample[parInd,:,testInt])), 20))
      if maximum((pVals_Nsample[parInd,:,testInt])) >pValTh
           axvspan((pValTh),maximum((pVals_Nsample[parInd,:,testInt])),color = "r",alpha = 0.1);
      end
           xscale("log");
      subplot(2,2,2);  parInd = 2;   plt[:hist]((pVals_Nsample[parInd,:,testInt]), bins=logspace(minimum(log10(pVals_Nsample[2,:,testInt])),maximum(log10(pVals_Nsample[parInd,:,testInt])), 20))
      if maximum((pVals_Nsample[parInd,:])) >pValTh
           axvspan((pValTh),maximum((pVals_Nsample[parInd,:,testInt])),color = "r",alpha = 0.1);
      end
           xscale("log");
      subplot(2,2,3); plot((pVals_Nsample[1,:,testInt]),(pVals_Nsample[2,:,testInt]),".")
      axvline((pValTh));axhline((pValTh))
      xscale("log");yscale("log")
      xlabel("P-Value " * namePar1)
      ylabel("P-Value " * namePar2)
end

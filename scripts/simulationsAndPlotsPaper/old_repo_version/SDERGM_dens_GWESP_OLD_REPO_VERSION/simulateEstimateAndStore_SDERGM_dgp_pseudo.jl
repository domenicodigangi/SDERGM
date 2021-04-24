
# sample sequences of ergms with different parameters' values from R package ergm
# and test the PseudoLikelihoodScoreDrivenERGM filter
using Utilities,AReg,StaticNets,DynNets , JLD,MLBase,StatsBase,CSV, RCall
 using PyCall; pygui(); using PyPlot
 using ErgmRcall


# load packages needed in R
R"rm(list = ls())
   library(statnet)
    library(ergm)
    library(sna)
    library(coda)
    library(network)
    sessionInfo()"

tot_Nsample = 200
N_each_chunk = 10
#divide the simulations in chuncks
for n_chunk = 1:Int(tot_Nsample/N_each_chunk)
T = 110
N = 100
UM = [-3 , 2.5]
B_Re  = 0.95
A_Re  = 0.05
indTvPar = trues(2)
indTargPar = trues(2)
Nsample = N_each_chunk
estPar_all = fill(fill(ones(2),2),Nsample)
obsT_all = fill(fill(ones(2,2),T),Nsample)

convFlag  = falses(Nsample)
n=1
while n <= Nsample
    # sample the dgp
    print(n)

    R"rm(list = ls())
       library(statnet)
        library(ergm)
        library(sna)
        library(coda)
        library(network)
        sessionInfo()"
    try
        obsT,statsT,fVecT = ErgmRcall.sdergmDGPpaper_Rcall(  DynNets.fooSdErgmDirBinGlobalPseudo, indTvPar ,T,N, UM = UM, B_Re = B_Re,A_Re = A_Re)
        model = DynNets.SdErgmDirBinGlobalPseudo(obsT,DynNets.fooGasPar,indTvPar,"")
        estPar,convFlag[n],UM_est,startPoint = DynNets.estimate(model;indTvPar = indTvPar, indTargPar = indTargPar)
        obsT_all[n] = obsT
        estPar_all[n] = estPar
        if !isnan(estPar[2][1])
            n = n+1
        end
    catch
        # clear all R variables
    end
end
# Save all Sampled Data
save_fold = "./data/estimatesTest/sdergmTest/distrib_static_pars/"
 @save(save_fold*"sdergm_dgp__$(N)_T_$(T)_Sample_$(Nsample)_UM_$(UM)_intTV_$(Int.(indTvPar[:]))_intTarg_$(Int.(indTargPar[:]))_N_each_chunk_$(N_each_chunk)_chunk_$(n_chunk).jld",
              obsT_all, estPar_all, UM,B_Re,A_Re,convFlag, Nsample,T,N, indTvPar, indTargPar )

end
#
#
# using Utilities,AReg,StaticNets,DynNets , JLD,MLBase,StatsBase,CSV, RCall
#  using PyCall; pygui(); using PyPlot
#  using JLD2,Utilities, GLM
#  ## load R MCMC simulation and estimates and estimate sdergmTest
#
#  if false
#  Nsample = 100
#  dgpType = "sin"
#  T = 50
#  N = 50
#  Nterms = 2
#  Nsteps1 ,Nsteps2 = 0,1
#  load_fold = "./data/estimatesTest/sdergmTest/R_MCMC_estimates/"
#  @load(load_fold*"test_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE.jld",
#              stats_T, changeStats_T,estParSS_T,sampledMat_T ,parDgpT,Nsample)
#  end
#
#  onlyTest = false
#  targetAllTv = true
#  tvParFromTest = false; pValTh = 0.05
#  gasParEst = fill(fill(Float64[],2), Nsample)
#  startPointEst = fill(fill(0.0,2), Nsample)
#  staticEst = fill(fill(0.0,2), Nsample)
#  filtPar_T_Nsample = zeros(Nterms,T,Nsample)
#  pVals_Nsample = zeros(Nterms,Nsample)
#  scoreAutoc_Nsample = zeros(Nterms,Nsample)
#  convFlag = falses(Nsample)
#  tmp = Array{Array{Float64,2},2}(T,Nsample); for t=1:T,s=1:Nsample tmp[t,s] =  changeStats_T[t][s];end;#changeStats_T = tmp
#  for n=1:Nsample
#      @show(n)
#       pVals_Nsample[:,n],tmpInfo,staticEst[n] =
#          DynNets.pValStatic_SDERGM( DynNets.SdErgmDirBinGlobalPseudo(tmp[:,n],DynNets.fooGasPar,falses(Nterms),""))
#       @show(tmpInfo)
#      if tvParFromTest
#          indTvPar =  pVals_Nsample[:,n] .< pValTh
#      else
#          indTvPar = BitArray([true,true])
#      end
#      if targetAllTv
#          indTargPar = indTvPar
#      else
#           indTargPar =BitArray([false,false])
#      end
#     model = DynNets.SdErgmDirBinGlobalPseudo(tmp[:,n],DynNets.fooGasPar,indTvPar,"")
#
#     if .!onlyTest
#     estPar,convFlag[n],UM,startPoint = DynNets.estimate(model;UM =meanSq(estParSS_T[n,:,:],1),indTvPar = indTvPar, indTargPar = indTargPar)
#
#
#     gasParEst[n] = estPar # store gas parameters
#     startPointEst[n] = startPoint
#     gasParVec = zeros(sum(indTvPar)*3); for i=0:(sum(indTvPar)-1) gasParVec[1+3i : 3(i+1)] = estPar[indTvPar][i+1]; end
#     constParVec = zeros(sum(!indTvPar)); for i=1:sum(.!indTvPar) constParVec[i] = estPar[.!indTvPar][i][1]; end
#     gasFiltPar , pseudolike = DynNets.score_driven_filter_or_dgp(model,gasParVec,indTvPar;
#                               vConstPar = constParVec,ftot_0= startPoint )
#
#
#     filtPar_T_Nsample[:,:,n] = gasFiltPar
#     end
#  end
#
#  save_fold = "./data/estimatesTest/sdergmTest/"*
#             "gas_MCMC_comparison_estimates/"
#  @save(save_fold*"test_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE_target_$(targetAllTv).jld" ,
#      stats_T, changeStats_T,estParSS_T,sampledMat_T ,parDgpT,
#      Nsample,T,N,filtPar_T_Nsample,gasParEst,convFlag,pVals_Nsample,scoreAutoc_Nsample,staticEst)
#
#
# if false
#
#     sum(convFlag)
#     figure()
#      legTex = ["Constant ERGM";"SDERGM"; "Sequence of ERGM"; "DGP"]
#      for n=1:Nsample
#      indTvPar = BitArray([true,true])
#      model = SdErgmDirBinGlobalPseudo(tmp[:,n],fooGasPar,indTvPar,"")
#      estPar = gasParEst[n] # store gas parameters
#      gasFiltPar  = filtPar_T_Nsample[:,:,n]
#      parInd = 1
#      subplot(1,2,1);plot(1:T,ones(T)*staticEst[n][parInd],"-b")
#                     plot(1:T,gasFiltPar[parInd,:],"r")
#                     plot(1:T,estParSS_T[n,:,parInd],".b")
#                     plot(1:T,parDgpT[parInd,:],"k",linewidth=5)
#
#      parInd = 2
#      subplot(1,2,2);plot(1:T,ones(T)*staticEst[n][parInd],"-b")
#                     plot(1:T,gasFiltPar[parInd,:],"r")
#                     plot(1:T,estParSS_T[n,:,parInd],".b")
#                     plot(1:T,parDgpT[parInd,:],"k",linewidth=5)
#      end
#      namePar1 = "Number of Links"
#      namePar2 = "GWESP"
#      subplot(1,2,1);    title(namePar1); legend(legTex)
#      subplot(1,2,2);  title(namePar2 ); legend(legTex)
#
#     #Plotta info relative ai tests
#     figure()
#          subplot(2,2,1); parInd = 1;   plt[:hist]((pVals_Nsample[parInd,:]), bins=logspace(minimum(log10(pVals_Nsample[1,:])),maximum(log10(pVals_Nsample[parInd,:])), 20))
#          if maximum((pVals_Nsample[parInd,:])) >pValTh
#               axvspan((pValTh),maximum((pVals_Nsample[parInd,:])),color = "r",alpha = 0.1);
#          end
#               xscale("log");
#          subplot(2,2,2);  parInd = 2;   plt[:hist]((pVals_Nsample[parInd,:]), bins=logspace(minimum(log10(pVals_Nsample[1,:])),maximum(log10(pVals_Nsample[parInd,:])), 20))
#          if maximum((pVals_Nsample[parInd,:])) >pValTh
#               axvspan((pValTh),maximum((pVals_Nsample[parInd,:])),color = "r",alpha = 0.1);
#          end
#               xscale("log");
#          subplot(2,2,3); plot((pVals_Nsample[1,:]),(pVals_Nsample[2,:]),".")
#          axvline((pValTh));axhline((pValTh))
#          xscale("log");yscale("log")
#          xlabel("P-Value " * namePar1)
#          ylabel("P-Value " * namePar2)
# end

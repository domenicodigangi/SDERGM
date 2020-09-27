#script that estimates the dir bin gas network model on emid data and evaluates GAS
# forecasting performances

using HelperFunDom,AReg,StaNets,JLD,MLBase,StatsBase,DynNets
using PyCall; pygui(:qt); using PyPlot
#estimate and save for half the datase (after LTRO) or whole data?

## Load dataj
using JLD
 halfPeriod = false
 fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
 loadFilePartialName = "Weekly_eMid_Data_from_"
 halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
 @load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks,
                YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)


#

#density plot
matY_T = YeMidWeekly_T[:,:,3:end]
 matA_T = matY_T.>0
    N = length(matA_T[1,:,1])
    T = size(matA_T)[3]
    Ttrain = 100#round(Int, T/2) #70 106 #
    threshVar = 0.00#5
    degsIO_T = [sumSq(matA_T,2);sumSq(matA_T,1)]
 allFitSS =  StaNets.estimate( StaNets.SnapSeqNetDirBin1(degsIO_T); identPost = false,identIter= true )
 sizeLab = 30
# density plot
# density_T = [sum(matA_T[:,:,t]) for t=1:T]./(N*(N-1))
# sizeLab = 30
#  close()
#  #subplot(2,1,1) ;
#  plot(density_T)
#  ylabel("Network Density",size = sizeLab)
#  #xlabel("Time",size = sizeLab)
#  grid()
#  plt[:axvline](x=Ttrain,color = "r",linestyle = "--")
#  S_T = [sum(matY_T[:,:,t]) for t=1:T]

#Estimate gas model on train sample
allFitConstTrain,~,~ =  StaNets.estimate( StaNets.NetModelDirBin1(meanSq(degsIO_T[:,1:Ttrain],2)) )
 modGasDirBin1_eMidTrain = DynNets.GasNetModelDirBin1(degsIO_T[:,1:Ttrain],"FISHER-DIAG")
 estTargDirBin1_eMidTrain,~ = DynNets.estimateTarg(modGasDirBin1_eMidTrain;SSest = allFitSS )
 gasParEstOnTrain = estTargDirBin1_eMidTrain
 modAllObs =  DynNets.GasNetModelDirBin1(degsIO_T,"FISHER-DIAG")
 GasforeFit,~ = DynNets.gasFilter(modAllObs,[gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])
 gasforeFit = Float64.(GasforeFit)

load_fold = "./data/multiStepForecast/"
 @load(load_fold* "aucForVariousNsampleAndNsteps_final.jld",
          aucGas, valsNsteps,valsNsample)

N_Nsteps,N_Nsample,Nforetype = size(aucGas)


#@load(save_fold* "aucForVariousNsampleAndNsteps_testMat.jld",
#       aucGas, valsNsteps,valsNsample)

sum(mean(expMat_T_means[:,:,2,1,:],3).>0.5)
close()
  tmp_means = [sum(expMat_T_means[:,:,2,1,i].>0.5) for i=1:N_Nsample]
  tmp_Y = [sum(expMat_T_Y[:,:,2,1,i].>0.5) for i=1:N_Nsample]
 subplot(2,1,1)
 indSteps = 1
 title("N steps =  $(valsNsteps[indSteps])")
 plot(valsNsample,tmp_means,".-")
 #plot(valsNsample,aucGas[indSteps,:,:],".-")
 subplot(2,1,2)
 indSam = 3
 title("N sample =  $(valsNsample[indSam])")
 plot(valsNsample,tmp_Y,".-")
 grid()
 #plot(tmp_Y)
 #plot(valsNsteps, aucGas[:,indSam,1:2],".-")

mean(tmp)
expMat_T[:,:,:,1] == expMat_T[:,:,:,4]

close()#
  remMat =squeeze(prod(.!matA_T,3),3)#
  sizeLeg = 15
  maxNsteps = 50
  legTex = []
  aucStat = zeros(maxNsteps)
    for n=1:maxNsteps
          shiftVal = n
          tmpTm1 = StaNets.nowCastEvalFitNet(   allFitSS[:,Ttrain+1:end],matA_T[:,:,Ttrain+1:end];mat2rem = remMat,shift = shiftVal)
          aucStat[n] = tmpTm1[3]
          legTex = [legTex; "t-$(shiftVal) Par  AUC = $(round(tmpTm1[3],3))"]
    end


  legend(legTex,fontsize = sizeLeg)
  title("T train = $(Ttrain) of $(T)  ",size = sizeLab)

  xlabel("False Positive Rate",size = sizeLab)
  ylabel("False Negative Rate",size = sizeLab)


close()

    legTex = "ERGM"
    Tfore = 16
    plot(1:Tfore,aucStat[1:Tfore,:],"r")
    for indSam=2
    plot(valsNsteps[1:Tfore],aucGas[1:Tfore,indSam,1])
     legTex = [legTex; "SDERGM " ]
    end
    legend(legTex,fontsize = sizeLeg)
    title("Multi Steps Ahead Forecast   ",size = sizeLab)
    xlabel("Number Of Steps Ahead",size = sizeLab)
    ylabel("Area Under The Curve",size = sizeLab)
    grid()




















#
# function forecastEvalGasNetDirBin1(obsNet_T::BitArray{3}, gasParEstOnTrain::Array{Array{Float64,1},1},Ttrain::Int,Nsteps::Int;Nsample::Int=100, thVarConst = 0.0005,mat2rem=falses(obsNet_T[:,:,1]) )
#     # this function evaluates the forecasting performances of the GAS model over the train sample.
#     # to do so uses gas parameters previously estimated and observations at time t of the test sample
#     # to forecast tv par at t+1, then using these estimates and t+1 observations goes to t+2 and so on
#
#     @show N2 = length(gasParEstOnTrain[1]);N = round(Int,N2/2)
#     T = length(obsNet_T[1,1,:])
#     @show Ttest = T-Ttrain
#     testObsNet_T = obsNet_T[:,:,Ttrain+1:end]
#     degsIO_T = [sumSq(obsNet_T,2); sumSq(obsNet_T,1)]
#
#
#
#     #disregard predictions for degrees that are constant in the training sample
#     # gas shoud manage but AR fitness no and I dont want to advantage gas
#     # to do so i need to count the number of links that are non constant (nnc)
#     Nc,isConIn,isConOut = StaNets.defineConstDegs(degsIO_T[:,1:Ttrain];thVarConst = thVarConst )
#     noDiagIndnnc = putZeroDiag(((.!isConIn).*(.!isConOut')).* (.!mat2rem)    )
#
#     @show Nlinksnnc =sum(noDiagIndnnc)
#     # forecast fitnesses using Gas parameters and observations
#     foreFit,~ = gasFilter( DynNets.GasNetModelDirBin1(degsIO_T),[gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])
#
#     TRoc = Ttest-1
#     #storage variables
#     foreVals = 100 * ones(Nsteps*Nlinksnnc*(TRoc))
#     realVals = trues( Nsteps*Nlinksnnc*(TRoc))
#
#     lastInd=1
#
#     for t=(1+Nsteps):TRoc
#         adjMat = testObsNet_T[:,:,t]   #[testObsNet_T[50:end,:,t];testObsNet_T[1:49,:,t]] #a shuffle to test that the code does not not work
#         adjMat_tm1 = testObsNet_T[:,:,t-1]
#         indZeroR = sum(adjMat_tm1,2).==0
#         indZeroC = sum(adjMat_tm1,1).==0
#         indZeroMat = indZeroR.*indZeroC
#     #    println(sum(indZeroMat)/(length(indZeroC)^2))
#
#         expMat = StaNets.expMatrix(StaNets.fooNetModelDirBin1,foreFit[:,Ttrain+1:end][:,t])
#
#          tmpAllMat,~ = multiSteps_updatedGasPar(modGasDirBin1_eMidTrain,N,foreFit[:,Ttrain+1:end][:,t-Nsteps],
#             gasParEstOnTrain[1],gasParEstOnTrain[2],gasParEstOnTrain[3],Nsample,Nsteps)
#
#          expMat=tmpAllMat[:,:,Nsteps+1]
#
#
#         #expMat[indZeroMat] = 0
#
#         foreVals[lastInd:lastInd+Nlinksnnc-1] = expMat[noDiagIndnnc]
#         realVals[lastInd:lastInd+Nlinksnnc-1] = adjMat[noDiagIndnnc]
#         lastInd += Nlinksnnc
#     end
#     return realVals,foreVals,foreFit
#  end

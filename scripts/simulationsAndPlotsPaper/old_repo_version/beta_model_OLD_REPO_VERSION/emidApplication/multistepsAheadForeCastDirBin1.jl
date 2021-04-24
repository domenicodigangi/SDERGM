#script that estimates the dir bin gas network model on emid data and evaluates GAS
# forecasting performances

using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase,DynNets
using PyCall; pygui(:qt); using PyPlot


## Load dataj
using JLD2
 halfPeriod = false
 fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
 loadFilePartialName = "Weekly_eMid_Data_from_"
 halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
 @load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks,
                YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)


#

#define sequence of matrices and degrees
matY_T = YeMidWeekly_T[:,:,3:end]
 matA_T = matY_T.>0
    N = length(matA_T[1,:,1])
    T = size(matA_T)[3]
    Ttrain = 100#round(Int, T/2) #70 106 #
    threshVar = 0.00#5
    degsIO_T = [sumSq(matA_T,2);sumSq(matA_T,1)]
 allFitSS =  StaticNets.estimate( StaticNets.SnapSeqNetDirBin1(degsIO_T); identPost = false,identIter= true )
 sizeLab = 30



 @load(fold_Path*"_rolling_SD_estimates_"*periodEndStr, gasFiltAndForeFitFromRollEst,allFitSS,gasEst_rolling)


initWorkers()
 using DynNets , StatsFuns
 @everywhere using DynNets, StatsFuns



#select links among largest banks in training sample
remMat =squeeze(prod(.!matA_T,3),3)#
  modAllObs = DynNets.SdErgmDirBin1(degsIO_T,"FISHER-DIAG")
  sizeLeg = 15
  legTex = []
  valsNsample = [5 10 20 30] #Vector(5:10:106)
  N_Nsample = length(valsNsample)
  maxNsteps = 25
  valsNsteps = Vector(1:maxNsteps) #[Vector(1:5); round(Int,linspace(6,maxNsteps,N_Nsteps-5))]
  N_Nsteps = length(valsNsteps)
  storeMat = false
  if storeMat
  expMat_T_means = SharedArray{Float64}(zeros(N,N,T,N_Nsteps,N_Nsample))
  expMat_T_Y = SharedArray{Float64}(zeros(N,N,T,N_Nsteps,N_Nsample))
  expMat_T_meanPar = SharedArray{Float64}(zeros(N,N,T,N_Nsteps,N_Nsample))
  else
      expMat_T_means = SharedArray{Float64}(zeros(N,N,T))
      expMat_T_Y = SharedArray{Float64}(zeros(N,N,T ))
      expMat_T_meanPar = SharedArray{Float64}(zeros(N,N,T))
  end
  aucGas = SharedArray{Float64}(zeros(N_Nsteps,N_Nsample,3))
 for n=1:N_Nsample
            @sync @parallel for k=1:N_Nsteps
              Nsteps = valsNsteps[k]
              Nsample = valsNsample[n]
              @show((Nsteps,Nsample))

             if storeMat
             @time expMat_T_means[:,:,:,k,n],forFit ,expMat_T_Y[:,:,:,k,n],expMat_T_meanPar[:,:,:,k,n], =
                   DynNets.multiStepsForecastExpMat_roll(modAllObs,matA_T,gasEst_rollin,gasFiltAndForeFitFromRollEst,Nsample,Nsteps,Ttrain)
               else
                   @time expMat_T_means,forFit ,expMat_T_Y,expMat_T_meanPar =
                             DynNets.multiStepsForecastExpMat_roll(modAllObs,matA_T,gasEst_rolling,gasFiltAndForeFitFromRollEst,Nsample,Nsteps,Ttrain)
             end
             #figure(11)
              tmpTm1 =  StaticNets.nowCastEvalFitNet(   zeros(2,2) ,matA_T[:,:,Ttrain+1:end];mat2rem = remMat,expMat_T=expMat_T_means[:,:,Ttrain+1:end],shift=0,plotFlag = false)

              aucGas[k,n,1] = tmpTm1[3]
              #figure(12)
              tmpTm1 =  StaticNets.nowCastEvalFitNet(   zeros(2,2) ,matA_T[:,:,Ttrain+1:end];mat2rem = remMat,expMat_T=expMat_T_Y[:,:,Ttrain+1:end],shift=0,plotFlag = false)

              aucGas[k,n,2] = tmpTm1[3]
             # figure(13)

              tmpTm1 =  StaticNets.nowCastEvalFitNet(   zeros(2,2) ,matA_T[:,:,Ttrain+1:end];mat2rem = remMat,expMat_T=expMat_T_meanPar[:,:,Ttrain+1:end],shift=0,plotFlag = false )
              aucGas[k,n,3] = tmpTm1[3]

        end
    end
 @sync @parallel for k=1:5
 end

 #close()
 #plot(valsNsteps,aucGas[:,1,1])


aucAR =  zeros(maxNsteps)
 ## multistep ahead forecast from AR1
 for Nsteps =1:maxNsteps
 ARForeFitFromRollEst = zeros(2N,T)
 foreFitAR_roll = zeros(2N,T)
 for t=Ttrain+1:T
    tEst = t-Ttrain
    tObs = t-Nsteps
    trainFit = allFitSS[:,tEst:tEst + Ttrain-1]
    Nconst,isConIn,isConOut = StaticNets.defineConstDegs(degsIO_T[:,tEst:tEst + Ttrain],thVarConst =0.005 )
    estAR1= StaticNets.estManyAR(trainFit;isCon = [isConIn;isConOut])
    ## Nsteps  ahead forecasts for each fitness for each time in Test sample
    for n=1:2N
        #compute forecasting for all links, diagonal will be disregarded after
        foreFitAR_roll[n,t] = StaticNets.oneStepForAR1(allFitSS[n,tObs],estAR1[n,:],steps = Nsteps)
    end
 end
  inds2forecast = (Ttrain+2:T)
  @show(Nsteps)
 tmpTm1 =  StaticNets.nowCastEvalFitNet(  foreFitAR_roll[:,inds2forecast] ,matA_T[:,:,inds2forecast];mat2rem = remMat,shift = 0,plotFlag=false)
 aucAR[Nsteps] = tmpTm1[3]
 end



close()#
  remMat =squeeze(prod(.!matA_T,3),3)#
  sizeLeg = 15
  maxNsteps = 50
  legTex = []
  aucStat = zeros(maxNsteps)
    for n=1:maxNsteps
          shiftVal = n
          tmpTm1 = StaticNets.nowCastEvalFitNet(   allFitSS[:,Ttrain+1:end],matA_T[:,:,Ttrain+1:end];mat2rem = remMat,shift = shiftVal,plotFlag = false)
          aucStat[n] = tmpTm1[3]
          legTex = [legTex; "t-$(shiftVal) Par  AUC = $(round(tmpTm1[3],3))"]
    end


save_fold = "./data/multiStepForecast/"
@save(save_fold* "aucForVariousNsampleAndNsteps_rolling_gas_estimates.jld",
              aucGas, valsNsteps,valsNsample,aucAR,aucStat)

  #
  # legend(legTex,fontsize = sizeLeg)
  # title("T train = $(Ttrain) of $(T)  ",size = sizeLab)
  #
  # xlabel("False Positive Rate",size = sizeLab)
  # ylabel("False Negative Rate",size = sizeLab)

#

if false
 close()
    legTex = []
    maxNsteps = 5
#    plot(1:maxNsteps,aucStat,"r")
    #legTex = "ERGM"
#    plot(1:maxNsteps,aucAR)
     #legTex = [legTex; "AR1" ]
    for indSam=1:N_Nsample
    plot(valsNsteps,aucGas[:,indSam,1])
     legTex = [legTex; "SDERGM avg mean mat N sample = $(valsNsample[indSam])" ]
    end
    legend(legTex,fontsize = sizeLeg)
    title("Multi Steps Ahead Forecast for ",size = sizeLab)
    xlabel("Number Of Steps Ahead",size = sizeLab)
    ylabel("Area Under The Curve",size = sizeLab)
    grid()
end



















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
#     Nc,isConIn,isConOut = StaticNets.defineConstDegs(degsIO_T[:,1:Ttrain];thVarConst = thVarConst )
#     noDiagIndnnc = putZeroDiag(((.!isConIn).*(.!isConOut')).* (.!mat2rem)    )
#
#     @show Nlinksnnc =sum(noDiagIndnnc)
#     # forecast fitnesses using Gas parameters and observations
#     foreFit,~ = score_driven_filter_or_dgp( DynNets.SdErgmDirBin1(degsIO_T),[gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])
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
#         expMat = StaticNets.expMatrix(StaticNets.fooErgmDirBin1,foreFit[:,Ttrain+1:end][:,t])
#
#          tmpAllMat,~ = multiSteps_predict_score_driven_par(modGasDirBin1_eMidTrain,N,foreFit[:,Ttrain+1:end][:,t-Nsteps],
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

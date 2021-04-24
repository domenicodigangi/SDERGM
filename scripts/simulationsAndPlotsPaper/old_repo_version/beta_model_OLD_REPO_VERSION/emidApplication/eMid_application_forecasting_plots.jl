#script that estimates the dir bin gas network model on emid data and evaluates GAS
# forecasting performances

using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase,DynNets
using PyCall; pygui(:qt); using PyPlot
#estimate and save for half the datase (after LTRO) or whole data?
halfPeriod = false

## Load dataj
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks,
                YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)


#density plot
matY_T = YeMidWeekly_T[:,:,3:end]
 matA_T = matY_T.>0
    N = length(matA_T[1,:,1])
    T = size(matA_T)[3]
    Ttrain = 100#round(Int, T/2) #70 106 #
    threshVar = 0.00#5
    degsIO_T = [sumSq(matA_T,2);sumSq(matA_T,1)]
 allFitSS =  StaticNets.estimate( StaticNets.SnapSeqNetDirBin1(degsIO_T); identPost = false,identIter= true )

density_T = [sum(matA_T[:,:,t]) for t=1:T]./(N*(N-1))

#Plot network density over time

close()
 sizeLab = 3
#  figure(figsize=(14,7))
   figure(figsize=(18,7))
  labSize = 26
  tickSize = 20
  legSize = 20
  #names = ["Edges" "GWESP" "Political Affiliation"]

  start_day = Dates.DateTime(2009, 6, 22)
  weeks = [start_day + Dates.Week(i) for i=1:T] #Dates.DateTime(2015, 2, 27),T);
  datestrings = Dates.format.(weeks, "u dd")
  plot(weeks,density_T)
  #plt[:axvline](x=weeks[Ttrain],color = "r",linestyle = "--")
  ylabel("Network Density",size =labSize)
  #plot(1:T,ones(T).*vConstPar[parInd],"-b")
  xticks(fontsize = tickSize)
  yticks(fontsize = tickSize)
  legTex = [ "Density "; "Separation Between Train and Test samples"]
  legend(legTex,fontsize = legSize,framealpha = 0.95)
  tight_layout()
  grid()
  S_T = [sum(matY_T[:,:,t]) for t=1:T]

  #text(weeks[30],0.012,"Train Sample",fontsize = labSize) #
  #text(weeks[Ttrain + 20],0.027,"Test Sample",fontsize = labSize) #


#Estimate gas model on train sample
allFitConstTrain,~,~ =  StaticNets.estimate( StaticNets.ErgmDirBin1(meanSq(degsIO_T[:,1:Ttrain],2)) )
 modGasDirBin1_eMidTrain = DynNets.SdErgmDirBin1(degsIO_T[:,1:Ttrain],"FISHER-DIAG")
 estTargDirBin1_eMidTrain,~ = DynNets.estimateTarg(modGasDirBin1_eMidTrain;SSest = allFitSS )
 gasParEstOnTrain = estTargDirBin1_eMidTrain
 GasforeFitTrain,~ = DynNets.score_driven_filter_or_dgp( DynNets.SdErgmDirBin1(degsIO_T),[gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])
 gasforeFitTrain = Float64.(GasforeFitTrain)


@load(fold_Path*"_rolling_SD_estimates_"*periodEndStr,gasFiltAndForeFitFromRollEst,allFitSS,gasEst_rolling)



### estimate AR1 on rolling window of cross sectional estimates and use those for 1 step ahead forecasts
Nsteps = 1
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


#foreEval = DynNets.forecastEvalGasNetDirBin1(matA_T,gasParEstOnTrain,Ttrain;thVarConst = threshVar,mat2rem = remMat)
close()
 figure(figsize=(9,7))
    labSize = 32
    tickSize = 24
    legSize = 24

 #select links among largest banks in training sample
 remMat =squeeze(prod(.!matA_T,3),3)
 inds2forecast = (Ttrain+2:T)
    tmpTm1 =  StaticNets.nowCastEvalFitNet(  gasFiltAndForeFitFromRollEst[:,inds2forecast] ,matA_T[:,:,inds2forecast];mat2rem = remMat,shift = 0)
    legTex = ["SD-ERGM   t-1 AUC = $(round(tmpTm1[3],3))"]
    tmpTm1 =  StaticNets.nowCastEvalFitNet(  foreFitAR_roll[:,inds2forecast] ,matA_T[:,:,inds2forecast];mat2rem = remMat,shift = 0)
    legTex = [legTex;"AR1     t-1 AUC = $(round(tmpTm1[3],3))"]
    # ROC for predictions from mean over Ttrain estimates
    #tmpTm1 = StaticNets.nowCastEvalFitNet(   allFitSS[:,inds2forecast] ,matA_T[:,:,inds2forecast];mat2rem = remMat)
    #legTex = [legTex; "ERGM Contemporaneous   AUC = $(round(tmpTm1[3],3))"]
    # tmpTm1 =  StaticNets.nowCastEvalFitNet(  gasforeFit[:,1:T] ,matA_T,Ttrain;mat2rem = remMat,shift = -1)
    # legTex = [legTex;"GAS t AUC = $(round(tmpTm1[3],3))"]
    shiftVal = 1
    tmpTm1 = StaticNets.nowCastEvalFitNet(   allFitSS[:,inds2forecast] ,matA_T[:,:,inds2forecast];mat2rem = remMat,shift = shiftVal)
    legTex = [legTex; "ERGM t-$(shiftVal)   AUC = $(round(tmpTm1[3],3))"]

    legend(legTex,fontsize = legSize)
    title(" ",size = labSize)
    #grid()
    xlabel("False Positive Rate",size = labSize)
    ylabel("False Negative Rate",size = labSize)
    xticks(fontsize = tickSize)
    yticks(fontsize = tickSize)
    xlim([0,0.6])
    tight_layout()


aucGas[1,1,1] = tmpTm1[3]
size(aucGas)
#
error()
save_fold = "./data/multiStepForecast/"
@load(save_fold* "aucForVariousNsampleAndNsteps_rolling_gas_estimates.jld",
       aucGas, valsNsteps,valsNsample,aucStat,aucAR)
@load(save_fold* "aucForVariousNsampleAndNsteps_gas_estimates.jld",
      aucGas, valsNsteps,valsNsample,aucStat,aucAR)




  # Multistep ahead forecast da rolling AR1 estimates


close()
  figure(figsize=(9,7))
  maxNsteps = 8
    plot(1:maxNsteps,aucStat[1:maxNsteps],"r")
    legTex = "ERGM"
    plot(1:maxNsteps,aucAR[1:maxNsteps])
     legTex = [legTex; "AR1" ]
   for i=4
    plot(1:maxNsteps,aucGas[1:maxNsteps,i,1])
     legTex = [legTex; "SD-ERGM " ]
   end
    legend(legTex,fontsize = legSize)
    title(" ",size = labSize)
    xlabel("Number Of Steps Ahead",size = labSize)
    ylabel("Area Under The Curve",size = labSize)
    grid()
    xticks(fontsize = tickSize)
    yticks(fontsize = tickSize)
    tight_layout()

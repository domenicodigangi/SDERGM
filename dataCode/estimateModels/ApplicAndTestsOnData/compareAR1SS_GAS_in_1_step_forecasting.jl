#script that estimates the dir bin gas network model on emid data and evaluates GAS
# forecasting performances

using HelperFunDom,AReg,StaNets,JLD,MLBase,StatsBase,DynNets
using PyCall; pygui(:qt); using PyPlot
#estimate and save for half the datase (after LTRO) or whole data?
halfPeriod = false

## Load dataj
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks,
                YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)


#

## Load  piero's DATA
using MAT
pieroEmidEst = matread("/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/matlabFiles/data\ Piero/forecastEMID/fitness_timeseries.mat")
pieroEmidData = matread("/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/matlabFiles/data\ Piero/forecastEMID/eMID_weeklyaggregated_binaryDIRECT_postLTRO.mat")

#density plot
matY_T = YeMidWeekly_T[:,:,3:end]
 matA_T = matY_T.>0
    N = length(matA_T[1,:,1])
    T = size(matA_T)[3]
    Ttrain = 100#round(Int, T/2) #70 106 #
    threshVar = 0.00#5
    degsIO_T = [sumSq(matA_T,2);sumSq(matA_T,1)]
 allFitSS =  StaNets.estimate( StaNets.SnapSeqNetDirBin1(degsIO_T); identPost = false,identIter= true )

density_T = [sum(matA_T[:,:,t]) for t=1:T]./(N*(N-1))
sizeLab = 30
 close()
 subplot(2,1,1) ;plot(density_T)
 ylabel("Network Density",size = sizeLab)
 #xlabel("Time",size = sizeLab)
 grid()
 plt[:axvline](x=Ttrain,color = "r",linestyle = "--")
 S_T = [sum(matY_T[:,:,t]) for t=1:T]
 subplot(2,1,2);plot(S_T)
  ylabel("Total Volume",size = sizeLab)
  xlabel("Time",size = sizeLab)
  grid()
 plt[:axvline](x=Ttrain,color = "r",linestyle = "--")

#Estimate gas model on train sample
allFitConstTrain,~,~ =  StaNets.estimate( StaNets.NetModelDirBin1(meanSq(degsIO_T[:,1:Ttrain],2)) )
 modGasDirBin1_eMidTrain = DynNets.GasNetModelDirBin1(degsIO_T[:,1:Ttrain],"FISHER-DIAG")
 estTargDirBin1_eMidTrain,~ = DynNets.estimateTarg(modGasDirBin1_eMidTrain;SSest = allFitSS )
 gasParEstOnTrain = estTargDirBin1_eMidTrain
 GasforeFit,~ = DynNets.gasFilter( DynNets.GasNetModelDirBin1(degsIO_T),[gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])
 gasforeFit = Float64.(GasforeFit)


#foreEval = DynNets.forecastEvalGasNetDirBin1(matA_T,gasParEstOnTrain,Ttrain;thVarConst = threshVar,mat2rem = remMat)
close()
 #select links among largest banks in training sample
 remMat =falses(N,N)#squeeze(prod(.!matA_T,3),3)#
 sizeLeg = 15
 tmpTm1 =  StaNets.nowCastEvalFitNet(  gasforeFit[:,1:T] ,matA_T,Ttrain;mat2rem = remMat,shift = 0)
    legTex = ["GAS t-1 AUC = $(round(tmpTm1[3],3))"]
    # ROC for predictions from mean over Ttrain estimates
    tmpTm1 = StaNets.nowCastEvalFitNet(   allFitSS[:,1:end] ,matA_T,Ttrain;mat2rem = remMat)
    legTex = [legTex; "t Par  AUC = $(round(tmpTm1[3],3))"]
    # tmpTm1 =  StaNets.nowCastEvalFitNet(  gasforeFit[:,1:T] ,matA_T,Ttrain;mat2rem = remMat,shift = -1)
    # legTex = [legTex;"GAS t AUC = $(round(tmpTm1[3],3))"]
    shiftVal = 1
    tmpTm1 = StaNets.nowCastEvalFitNet(   allFitSS[:,1:end] ,matA_T,Ttrain;mat2rem = remMat,shift = shiftVal)
    legTex = [legTex; "t-$(shiftVal) Par  AUC = $(round(tmpTm1[3],3))"]

    tmpTm1 = StaNets.nowCastEvalFitNet( repmat( allFitConstTrain,1,T),matA_T,Ttrain;mat2rem = remMat)
    legTex = [legTex; "Const Par  AUC = $(round(tmpTm1[3],3))"]

    tmpTm1 = StaNets.forecastEvalAR1Net( matA_T,Ttrain,mat2rem = remMat)
    legTex = [legTex; "AR1 SS = $(round(tmpTm1[3],3))"]
    legend(legTex,fontsize = sizeLeg)
    title("T train = $(Ttrain) of $(T)  ",size = sizeLab)

    xlabel("False Positive Rate",size = sizeLab)
    ylabel("False Negative Rate",size = sizeLab)



# Piero estimates of the fitnesses
N2,T = size(pieroEmidEst["thetaEST0"]);N =round(Int,N2/2);TtrainP = T - 50
 fitP1_T = -pieroEmidEst["thetaEST0x"]# [pieroEmidEst["thetaEST0"][1:N,:];pieroEmidEst["thetaEST0"][1:N,:]]
 fitP2_T = -pieroEmidEst["thetaEST0"]
 matAP_T = pieroEmidData["A"].>0
 shiftVal = 1
 mat2rem = falses(N,N)# squeeze(prod(.!matAP_T,3),3)
 tmpTm1 = StaNets.nowCastEvalFitNet(fitP2_T,matAP_T,TtrainP,shift = shiftVal,mat2rem=mat2rem)
 legTex = ["theta0 t - $(shiftVal)   AUC = $(round(tmpTm1[3],3))"]
 tmpTm1 = StaNets.nowCastEvalFitNet(fitP1_T,matAP_T,TtrainP,shift = shiftVal,mat2rem=mat2rem)
 legTex = [legTex ; "theta0X t - $(shiftVal)   AUC = $(round(tmpTm1[3],3))"]
 legend(legTex)


N = length(matAP_T[1,:,1])
    T = size(matAP_T)[3]
    Ttrain = 106#round(Int, T/2) #70 106 #
    threshVar = 0.00#5
    degsIOP_T = [sumSq(matAP_T,2);sumSq(matAP_T,1)]

 allFitSS =   StaNets.estimate( StaNets.SnapSeqNetDirBin1(degsIOP_T); identPost = false,identIter= true,targetErr = 1e-6 )
 allFitConstTrain,~,~ =  StaNets.estimate( StaNets.NetModelDirBin1(meanSq(degsIOP_T[:,1:Ttrain],2)) )
 modGasDirBin1_eMidTrain = DynNets.GasNetModelDirBin1(degsIOP_T[:,1:Ttrain],"FISHER-DIAG")
 estTargDirBin1_eMidTrain,~ = DynNets.estimateTarg(modGasDirBin1_eMidTrain;SSest = allFitSS )
 gasParEstOnTrain = estTargDirBin1_eMidTrain
 GasforeFit,~ = DynNets.gasFilter( DynNets.GasNetModelDirBin1(degsIOP_T),
                    [gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])

 #Nconst,isConIn,isConOut = defineConstDegs(degsIOP_T[:,1:Ttrain];thVarConst = thVarConst )

 #foreEval = DynNets.forecastEvalGasNetDirBin1(matAP_T,gasParEstOnTrain,Ttrain;thVarConst = threshVar,mat2rem = remMat)
 close()
    #select links among largest banks in training sample
    remMat = mat2rem# falses(N,N)#;for i=1:N remMat[i,i] = true; end
    tmpGas =  StaNets.nowCastEvalFitNet(   Float64.(GasforeFit) ,matAP_T,Ttrain;mat2rem = remMat,shift=0)
     legTex = ["GAS  AUC = $(round(tmpGas[3],3))"]
     legend(legTex)
     #tmpAR = StaNets.forecastEvalAR1Net(matAP_T,Ttrain;thVarConst = threshVar,mat2rem = remMat)
     #legTex = [legTex ; "AR1  AUC = $(round(tmpAR[3],3))"]
     legend(legTex)
     title("T train = $(Ttrain) of $(T)  ")
       tmpTm1 = StaNets.nowCastEvalFitNet(   allFitSS[:,1:end] ,matAP_T,Ttrain;mat2rem = remMat)
     legTex = [legTex; "t Par  AUC = $(round(tmpTm1[3],3))"]
     # tmpTm1 = StaNets.nowCastEvalFitNet(  [ allFitSS[:,1]  allFitSS[:,1:end-1] ] ,matAP_T,Ttrain;mat2rem = remMat)
     # legTex = [legTex; "t-1 Par  AUC = $(round(tmpTm1[3],3))"]
     shiftVal = 1
     tmpTm1 = StaNets.nowCastEvalFitNet(    allFitSS[:,1:end],matAP_T,Ttrain;mat2rem = remMat,shift=shiftVal)
     legTex = [legTex; "t-$(shiftVal) Par  AUC = $(round(tmpTm1[3],3))"]

     tmpTm1 = StaNets.nowCastEvalFitNet( repmat( allFitConstTrain,1,T),matAP_T,Ttrain;mat2rem = remMat)
     legTex = [legTex; "Const Par  AUC = $(round(tmpTm1[3],3))"]

     tmpTm1 = StaNets.forecastEvalAR1Net( matAP_T,Ttrain;mat2rem = remMat)
     legTex = [legTex; "AR1 SS = $(round(tmpTm1[3],3))"]
     legend(legTex)
     title("T train = $(Ttrain) of $(T)  ")
     grid()
     xlabel("False Positive Rate")
     ylabel("False Negative Rate")


using RCall

function Rroc(score,response)
    R"library('pROC')"
        @rput(score)
        @rput(response)
        R" plot(roc(response, score , direction='<'), print.auc=TRUE)
        par(new=TRUE)
        "
end

Rroc(tmpGas[5], tmpGas[4])
Rroc(tmpTm1[5], tmpTm1[4])



#

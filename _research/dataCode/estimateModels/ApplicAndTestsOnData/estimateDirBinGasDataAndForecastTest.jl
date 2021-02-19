#script that estimates the dir bin gas network model on emid data and evaluates GAS
# forecasting performances

using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase
using PyCall; pygui(:qt); using PyPlot
using DynNets
#estimate and save for half the datase (after LTRO) or whole data?
halfPeriod = true

## Load data
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod ? periodEndStr =  "2012_03_12_to_2015_02_27.jld" : periodEndStr =  "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)

#
# ## Binary  estimates
# divide the data in Train and test samples and estimate the gas model on the train
N2,T = size(degsIO_T);N = round(Int,N2/2)
Ttrain = round(Int,T/2)
Ttest = T - Ttrain
modGasDirBin1_eMid = DynNets.GasNetModelDirBin1(degsIO_T)
modGasDirBin1_eMidTrain = DynNets.GasNetModelDirBin1(degsIO_T[:,1:Ttrain])
SSestForTargTrain = DynNets.estimateSnapSeq(modGasDirBin1_eMidTrain)
# ## test targeted estimate

estTargDirBin1_eMidTrain,~ = DynNets.estimateTarg(modGasDirBin1_eMidTrain;SSest = SSestForTargTrain)
#estTargDirBin1_eMid,~ = DynNets.estimateTarg(modGasDirBin1_eMid)

# ## save estimates
#@save(fold_Path*"EstDirBin1Gas"*loadFilePartialName*periodEndStr,estTargDirBin1_eMidTrain,modGasDirBin1_eMidTrain,SSestForTargTrain)
##load the estimates instead of re estimating every time
#@load(fold_Path*"EstDirBin1Gas"*loadFilePartialName*periodEndStr,estTargDirBin1_eMidTrain,SSestForTargTrain)

## Forecast functions


N2,T = size(degsIO_T);N = round(Int,N2/2)

t=1

ftotIO_t_Start = identify!(modGasDirBin1_eMid,gasPar[1].*(1-gasPar[2]))
tvPartp1 = predict_score_driven_par( modGasDirBin1_eMid,N,degsIO_T[:,t],
                         copy(ftotIO_t_Start),zeros(N2,N2),trues(N2),
                         identify!(modGasDirBin1_eMid,gasPar[1]),gasPar[2],gasPar[3])[1]


## --------------Arrivato qui al test della funzione che valuta il foreasting
gasParEstOnTrain = estTargDirBin1_eMidTrain
testObsNet_T = AeMidWeekly_T[:,:,Ttrain+1:end]

foreEval = forecastEvalGasNetDirBin1(gasParEstOnTrain,testObsNet_T)
rocCurve(foreEval[1],foreEval[2])

#check that forecasted fitnesses and degrees are reasonable
foreFit,~ = score_driven_filter_or_dgp( DynNets.GasNetModelDirBin1(degsIO_T),[gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])
plot(foreFit')
expDegsIO_Ttest = zeros(N2,Ttest)
for t=1:Ttest
    expDegsIO_Ttest[:,t] = [sumSq(StaticNets.expMatrix(StaticNets.fooNetModelDirBin1,foreEval[3][:,t]),2);sumSq(StaticNets.expMatrix(StaticNets.fooNetModelDirBin1,foreEval[3][:,t]),1)]
end
plot(expDegsIO_Ttest' )



## Wrap everithing up in a function









#













#

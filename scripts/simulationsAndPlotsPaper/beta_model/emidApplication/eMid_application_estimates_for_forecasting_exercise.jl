#script that estimates the dir bin gas network model on emid data and evaluates GAS
# forecasting performances

using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase,DynNets
using PyCall; pygui(:qt); using PyPlot
#estimate and save for half the datase (after LTRO) or whole data?
halfPeriod = false

## Load dataj
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod ? periodEndStr =  "2012_03_12_to_2015_02_27.jld" : periodEndStr =  "2009_06_22_to_2015_02_27.jld"
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
 allFitSS =  StaticNets.estimate( StaticNets.SnapSeqNetDirBin1(degsIO_T);
                        identPost = false,identIter= true )


gasEst_rolling = Array{Vector{Vector{Float64}},1}(T-Ttrain)
#the following variable stores the filtered values of the Parameters
# using at each step the observation t and static gas parameters estimated
# using observations from t-101:t-1
# te values between 1:Ttrain are filtered using the contstant Parameters estimated
# in the first Ttrain observations
gasFiltAndForeFitFromRollEst = zeros(2N,T)
#Estimate gas model on rolling windows of length Ttrain
modGasDirBin1_eMidTrain = DynNets.GasNetModelDirBin1(degsIO_T[:,1:Ttrain],"FISHER-DIAG")
gasParEstOnTrain_roll_t,~ = DynNets.estimateTarg(modGasDirBin1_eMidTrain;SSest = allFitSS[:,1:Ttrain] )
gasEst_rolling[1] = gasParEstOnTrain_roll_t
GasforeFit,~ = DynNets.gasFilter( DynNets.GasNetModelDirBin1(degsIO_T[:,1:Ttrain+1]),
                                    [gasParEstOnTrain_roll_t[1];gasParEstOnTrain_roll_t[2];gasParEstOnTrain_roll_t[3]])
gasforeFit = Float64.(GasforeFit)
gasFiltAndForeFitFromRollEst[:,1:Ttrain+1] = gasforeFit[:,1:Ttrain+1]
#

for t = 2:T-Ttrain
    @show(t)
 modGasDirBin1_eMidTrain = DynNets.GasNetModelDirBin1(degsIO_T[:,t:Ttrain+t],"FISHER-DIAG")
 gasParEstOnTrain_roll_t,~ = DynNets.estimateTarg(modGasDirBin1_eMidTrain;SSest = allFitSS[:,1+t:Ttrain+t] )
 gasEst_rolling[t] = gasParEstOnTrain_roll_t
 GasforeFit,~ = DynNets.gasFilter( DynNets.GasNetModelDirBin1(degsIO_T[:,t:Ttrain+t]),
                                     [gasParEstOnTrain_roll_t[1];gasParEstOnTrain_roll_t[2];gasParEstOnTrain_roll_t[3]])
 gasforeFit = Float64.(GasforeFit)
 # add the last value from this sequence of filtered ones (for which the rolling gas estimates are used) to the sequence of filtered parameters
 gasFiltAndForeFitFromRollEst[:,Ttrain+t] = gasforeFit[:,Ttrain+1]
 end

@save(fold_Path*"_rolling_SD_estimates_"*periodEndStr, gasFiltAndForeFitFromRollEst,allFitSS,gasEst_rolling)

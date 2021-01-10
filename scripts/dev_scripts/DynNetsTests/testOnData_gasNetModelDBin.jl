using Utilities, DynNets,JLD

halfPeriod = false
## Load data
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)

N2,T = size(degsIO_T);N = round(Int,N2/2)

## Binary  estimates tests
eMidMod = DynNets.GasNetModelDirBin1(degsIO_T)
SSest = estimateSnapSeq(eMidMod)
## test targeted estimate
@time estimateTarg(eMidMod;SSest = SSest)

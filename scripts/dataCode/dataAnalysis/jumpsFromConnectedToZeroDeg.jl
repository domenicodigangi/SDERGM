using HelperFunDom,StaNets,JLD,StatsFuns,StatsBase
using PyCall; pygui(:qt); using PyPlot


## Load data
#whole period data file
loadFileName = "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/Weekly_eMid_Data_from_2009_06_22_to_2015_02_27.jld"
@load(loadFileName, AeMidWeekly_Active_T,banksIDs,inactiveBanks, YeMidWeekly_Active_T,weekInd,datesONeMid ,degsIO_T,strIO_T)

##
inDegs=sumSq(AeMidWeekly_Active_T,2)
outDegs = sumSq(AeMidWeekly_Active_T,1)
N,T = size(inDegs)
NtimesT = prod(size(inDegs))

# An AR fitness model might have troubles explaining lots of jumps from highly active to inactive
# To check for it let's #Identify jumps from active to Inactive and relative entity
inactiveIn = (inDegs.==0)
degBeforeInactiveIn = inDegs[:,1:end-1][inactiveIn[:,2:end]]
freqIn = [(0:maximum(degBeforeInactiveIn)) counts(degBeforeInactiveIn.+1,maximum(degBeforeInactiveIn.+1))]
inactiveOut = (outDegs.==0)
degBeforeInactiveOut = outDegs[:,1:end-1][inactiveOut[:,2:end]]
freqOut = [(0:maximum(degBeforeInactiveOut)) counts(degBeforeInactiveOut.+1,maximum(degBeforeInactiveOut.+1))]
#Fractions of zero degrees observations
(sum(inactiveIn[:,2:end])/(NtimesT),sum(inactiveOut[:,2:end])/(NtimesT))
# Fraction  of jumps from a given degree(greater than zero) to zero
plot(freqOut[2:end,1],freqOut[2:end,2]/NtimesT,"r. ");plot(freqIn[2:end,1],freqIn[2:end,2]/NtimesT,"b. ");grid()
# plot fraction of jumps to zero starting from a degree equal or lower than that on x axis
plot(cumsum(freqIn[:,2]./sum(freqIn[:,2])),"r. "); plot(cumsum(freqOut[:,2]./sum(freqOut[:,2])),"b+ "); grid( )
# In summary the number of zero degree observations is quite high but
# those going to zero from a large degree are ~5% for jumps larger than 2 and  <2% for 5

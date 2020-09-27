
# script that wants to numerically test chatteris diaconis (misspelled with 99% prob)
# for the estimates of beta, fitness,ergm (many names..) in the DirBin1 case
using HelperFunDom,AReg,StaNets,JLD,MLBase,StatsBase,KernelDensity,DynNets
using PyCall; pygui(); using PyPlot

Nsample = 100


# load my data
halfPeriod = false
## Load dataj
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)

matY_T = YeMidWeekly_T
strIO_T = [sumSq(matY_T,2);sumSq(matY_T,1)]
N,~,T = size(matY_T)

i = 2
nnzStrIO = [ strIO_T[i,strIO_T[i,:].>0] for i =1:2N]
nnzLenghtIO = [length(nnzStrIO[i]) for i=1:2N]
indNnzStrIO = nnzLenghtIO.> 1
acorrNnzStrIO = [autocor(nnzStrIO[indNnzStrIO][i],[1])[1] for i=1:sum(indNnzStrIO)]

plot(nnzLenghtIO,".")
plot(strIO_T[i,strIO_T[i,:].>0])
plot(acorrNnzStrIO)

fracsI = Float64[]
 fracsO = Float64[]
 for t=1:T
    Y = YeMidWeekly_T[:,:,t]
    YfracIn = Y./sum(Y,2)
    YfracOut = Y./sum(Y,1)

    tmpO = YfracIn[.!isnan(YfracIn)][YfracIn[.!isnan(YfracIn)].>0]
    tmpI = YfracOut[.!isnan(YfracOut)][YfracOut[.!isnan(YfracOut)].>0]
    fracsI = [fracsI; tmpI]
    fracsO = [fracsO; tmpO]
end

xGrid = linspace(0,1,10000)
 plt[:hist](fracsI,100)

uniI = unique(fracsI )
countUniI = [sum(uniI[i] .== fracsI) for i=1:length(uniI)]
plot(uniI,countUniI,".")
tabI = sortrows([countUniI uniI ],rev = true)
sum(countUniI[countUniI.>10])./sum(countUniI)
semilogx(tabI[:,1])
prevDegEqFlag = (degsIO_T[:,1:end-1] .== degsIO_T[:,2:end])
prevDegEqFlagI,prevDegEqFlagO = prevDegEqFlag[1:N,:], prevDegEqFlag[1+1N:end,:]
tmp = (sum(prevDegEqFlagI,1) - N).*(sum(prevDegEqFlagO,1) - N)/(N*(N-1))
PyPlot.plot(tmp')

steps = 1
 prevDegDiff = abs.(degsIO_T[:,1:end-steps] - degsIO_T[:,steps+1:end])
 sum(prevDegDiff.<=1)
 degDiffMax = 20
 plt[:hist](prevDegDiff[prevDegDiff.<=degDiffMax],degDiffMax,density = true);grid()
 title("Absolute Values of One Step Degree Variation ")
 ylabel("Fraction of Occurences")
# plot(1:degDiffMax,cumsum(counts(prevDegDiff[:]))[1:degDiffMax]/sum(counts(prevDegDiff[:])))


prevStrDiff = (strIO_T[:,1:end-1] - strIO_T[:,2:end])[:]

prevStrDiffNnz = zeros(1); for i=1:2N  prevStrDiffNnz = vcat(prevStrDiffNnz , diff(strIO_T[i,:][strIO_T[i,:].!=0]) ); end; deleteat!(prevStrDiffNnz,1)
posDiffNnz = prevStrDiffNnz[prevStrDiffNnz.>0]
negDiffNnz = prevStrDiffNnz[prevStrDiffNnz.<0]

 diff(strIO_T[i,:][strIO_T[i,:].!=0])./strIO_T[i,:][strIO_T[i,:].!=0][2:end]

strIO_T

i =12
log(strIO_T[i,:][strIO_T[i,:].!=0][2:end]./strIO_T[i,:][strIO_T[i,:].!=0][1:end-1] )

prevStrRelDiffNnz = zeros(1); for i=1:2N  prevStrRelDiffNnz = vcat(prevStrRelDiffNnz , log(strIO_T[i,:][strIO_T[i,:].!=0][2:end]./strIO_T[i,:][strIO_T[i,:].!=0][1:end-1] )); end;
deleteat!(prevStrRelDiffNnz,1)
plot(prevStrRelDiffNnz[:])
using AverageShiftedHistograms,Plots
plotly()
Plots.plot(ash(prevStrRelDiffNnz),leg=false,size = (1200,700))
xlabel!("log( str_tm1/str_t)")
grid()
posDiff = prevStrDiff[prevStrDiff.>0]
negDiff = prevStrDiff[prevStrDiff.<0]
zeroDiff = prevStrDiff[prevStrDiff.==0]
posecdf = ecdf(posDiff)
negecdf = ecdf(abs.(negDiff))
posecdfNnz = ecdf(posDiffNnz)
negecdfNnz = ecdf(negDiffNnz)

Npoints = 100
tot = length(prevStrDiff)
x = exp(linspace(0,log(maximum(posDiffNnz)),Npoints))
 loglog(x,(1-posecdf(x))*(length(posDiffNnz)/tot))
 grid(which = :both ,alpha = 0.5)
 x = exp(linspace(0,log(maximum(abs.(negDiffNnz))),Npoints))
  loglog(x,(1-negecdf(x))*(length(negDiffNnz)/tot))
  grid(which = :both ,alpha = 0.5)

degDiffMax = 20
 plt[:hist](prevDegDiff[prevDegDiff.<=degDiffMax],degDiffMax,density = true);grid()
 title("Absolute Values of One Step Degree Variation ")
 ylabel("Fraction of Occurences")

prevDegEqFlagI,prevDegEqFlagO = prevDegEqFlag[1:N,:], prevDegEqFlag[1+1N:end,:]

tmp = (sum(prevDegEqFlagI,1) - N).*(sum(prevDegEqFlagO,1) - N)/(N*(N-1))
plot(tmp')

##
using HelperFunDom
plot(mean(degsIO_T,2),autocor(degsIO_T',[1],demean=true)',".")

nnzDegs_T  = degsIO_T[sumSq(degsIO_T.>0,2).>T/10,:]

plot(sum(prevDegDiff.>0,2),autocor(degsIO_T',[20 ],demean=true)',".")
 grid()
autocor(randi(0,30,297),[1])
#
lagsVec = [1]
maxDeg = 20
 T = 297
 Nsample = 100
 sampleAutocorr = zeros(Nsample,2)
 for i=1:Nsample
    sampleAutocorr[i,:] = [autocor([randi(maxDeg) for i=1:T ],lagsVec) autocor(maxDeg*rand(T) ,lagsVec)]
 end
 plot(sampleAutocorr)


tmp = ones(100)
 tmp[1 ] = 12
 autocor(tmp,lagsVec)






#

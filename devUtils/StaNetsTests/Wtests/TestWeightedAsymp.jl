
# script that wants to numerically test chatteris diaconis (misspelled with 99% prob)
# for the estimates of beta, fitness,ergm (many names..) in the DirBin1 case
using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase#,DynNets
using PyCall; pygui(:qt); using PyPlot

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

t=1
(sum(strIO_T[1:N,t]) - sum(strIO_T[N+1:end,t]) )< 10^Base.eps()? S =sum(strIO_T[N+1:end,t]) :error()
SperLink  = S/(N*N-N)
tmp = -0.5log(SperLink)

1/exp(2tmp)
1/exp(tmp)
plot(strIO_T')
lstrI = log.(strIO_T[1:N,t])
lstrO = log.(strIO_T[1+N:end,t])
1./exp(-lstrI)
expMat = expMatrix(fooNetModelDirW1Afixed,-[lstrI;lstrO],trues(N,N))

prevDegEqFlag = (degsIO_T[:,1:end-1] .== degsIO_T[:,2:end])
prevDegEqFlagI,prevDegEqFlagO = prevDegEqFlag[1:N,:], prevDegEqFlag[1+1N:end,:]
tmp = (sum(prevDegEqFlagI,1) - N).*(sum(prevDegEqFlagO,1) - N)/(N*(N-1))
plot(tmp')

prevDegDiff = abs.(degsIO_T[:,1:end-1] - degsIO_T[:,2:end])
sum(prevDegDiff.<=1)
degDiffMax = 20
plt[:hist](prevDegDiff[prevDegDiff.<=degDiffMax],degDiffMax,density = true);grid()
 title("Absolute Values of One Step Degree Variation ")
 ylabel("Fraction of Occurences")
# plot(1:degDiffMax,cumsum(counts(prevDegDiff[:]))[1:degDiffMax]/sum(counts(prevDegDiff[:])))

prevDegEqFlagI,prevDegEqFlagO = prevDegEqFlag[1:N,:], prevDegEqFlag[1+1N:end,:]

tmp = (sum(prevDegEqFlagI,1) - N).*(sum(prevDegEqFlagO,1) - N)/(N*(N-1))
plot(tmp')
allFitSS =  StaticNetsW.estimate( StaticNetsW.SnapSeqNetDirW1(strIO_T./10000000); identPost = false,identIter= true )
dgpPar = meanSq(allFitSS,2 )
expDgpMat = StaticNets.expMatrix(fooNetModelDirBin1,dgpPar)
dgpDegs = [sumSq(expDgpMat,2); sumSq(expDgpMat,1)]
N = length(matA_T[1,:,1])
Nsample = 1000
L = maximum(abs.(dgpPar))
C = L
errBnd = C*sqrt(log(N)/N )
prob = 1 - C/(N^2)
dgpMod = NetModelDirBin1(zeros(dgpDegs),dgpPar,Int.(1:N))
findmax(dgpDegs)
#sample the dgp model
sam = sampl(dgpMod,Nsample )

est = zeros(2N,Nsample)
degsIO_Sample = [sumSq(sam,2);sumSq(sam,1)]
estPar = zeros(2N,Nsample)
estDegs = zeros(2N,Nsample)
#estimate
for i=1:Nsample
    println(i)
    estPar[:,i]= estimate(fooNetModelDirBin1;degIO =degsIO_Sample[:,i])[1]
    expEstMat = StaticNets.expMatrix(fooNetModelDirBin1,estPar[:,i])
    estDegs[:,i] = [sumSq(expEstMat,2); sumSq(expEstMat,1)]

end

estPar[.!isfinite(estPar)] = -5

diffdat = estDegs .- dgpDegs # abs.(estPar .- dgpPar)
plt[:hist](diffdat[:]  ,100,density=true )
 title("Sample  Static Model, Run MLE and Conpare Degrees")
 grid(which = "both" , linewidth=1,alpha = 0.5)
 xlim(-3, 3)
 xlabel("deg_dgp - deg_MLE")

sampMeanDegs = meanSq(degsIO_Sample,2)
ind1= findmin(sampMeanDegs)[2]
ind2 = findmax(sampMeanDegs)[2]
diffdat = abs.(estPar .- dgpPar)./abs.(dgpPar)
plt[:hist]([diffdat[ind1,:] diffdat[ind2,:]]  ,50,density=true )
 title("DGP fit = mean sing snap estimates  ")
 legend(["Sample Mean Deg = $(round(sampMeanDegs[ind1],2)) " ; "Sample Mean Deg = $(round(sampMeanDegs[ind2],2)) "])
  xlabel("θ_est - θ_dgp  ")
 grid(which = "both" , linewidth=1,alpha = 0.5)

degsIO_Sample[1,:]
degsIO_Sample[ind2,:]


funProb(x::Array{T,1} where T<: Real ) = 1./(1+exp.(-x))
derfunProb(x::Array{T,1} where T<: Real )  =  exp.(-x)./((1+exp.(-x)).^2)
x = Vector(-5:0.001:5)
othPar = -3:1.5:3
for i = 1:length(othPar) plot(x, funProb(x.+othPar[i] )); end
 legTex = ["θ_2 =  $(othPar[i])" for i = 1:length(othPar)]
 title( "1/(1+ exp( - θ - θ_2 ) )")
 legend(legTex)
 xlabel("θ")
 grid()










#

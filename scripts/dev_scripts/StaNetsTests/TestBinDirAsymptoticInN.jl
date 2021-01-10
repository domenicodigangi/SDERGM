
# script that wants to numerically test chatteris diaconis (misspelled with 99% prob)
# for the estimates of beta, fitness,ergm (many names..) in the DirBin1 case
using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase,DynNets
using PyCall; pygui(:qt); using PyPlot

Nsample = 100
sizes = [50 100  200] #400 800 1200]
Nsizes = length(sizes)
estParSiz = Array{Array{Float64,2},1}(Nsizes)
dgpParSiz = Array{Array{Float64,1},1}(Nsizes)

for k=1:Nsizes
    ##Define the Dgp model
    N = sizes[k]
    deltaN = 5
    maxDeg = sqrt(N)
    tmpdegsIO = round.([ linspace(deltaN,maxDeg,N);Vector(linspace(deltaN,maxDeg,N))])
    dgpPar,dgpDegs, = estimate(NetModelDirBin1(tmpdegsIO))
    #remove infinities
    dgpPar[.!isfinite(dgpPar)] = -3

    dgpParSiz[k] = dgpPar

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
    #estimate
    for i=1:Nsample
        println(i)
        estPar[:,i]= estimate(fooNetModelDirBin1;degIO =degsIO_Sample[:,i])[1]
    end

    estPar[.!isfinite(estPar)] = -3
    estParSiz[k] = estPar
end

k = 1

diff =[ abs.(estParSiz[k] .- dgpParSiz[k]) for k=1:Nsizes]
plt[:hist]([diff[i][:] for i = 1:Nsizes] ,15,density=true )
 title(" A First Test for Asymptotic Consistency    ")
 legend(["N = $(sizes[k])  Bnd = C(L) $(round(sqrt(log(sizes[k])/sizes[k]),3))" for k=1:Nsizes])
  xlabel("θ_est - θ_dgp  ")
 grid(which = "both" , linewidth=1,alpha = 0.5)
##
Mod = NetModelDirBin1(deg,zeros(N),groupsInds)
groupsNames = unique(groupsInds)
NG = length(groupsNames)
tmpDic = countmap(groupsInds)
#number of nodes per group
NGMemb = sortrows([[i for i in keys(tmpDic)] [  tmpDic[i] for i in keys(tmpDic) ] ])
# agrregate degree of each group
degGroups = zeros(2NG) ; for i =1:N degGroups[groupsInds[i]] += deg[i]; degGroups[NG + groupsInds[i]] += deg[N + i];  end
L = sum(degGroups,1)[1]
L/(N^2-N)
groupsPar = ones(2NG)
firstOrderCond(Mod;degIO = deg,parGroupsIO = groupsPar, groupsInds = groupsInds )
expMatrix(Mod,groupsPar)

targetErrValStaticNets = 0.005
estPar,estIt,estMod = estimate(Mod)

##
sam = sampl(estMod,1000)

meanMatSam = squeeze(mean(sam,3),3)
degsSam = expValStatsFromMat(estMod,meanMatSam)
degGroupsSam = zeros(2NG) ; for i =1:N degGroupsSam[groupsInds[i]] += degsSam[i]; degGroupsSam[NG + groupsInds[i]] += degsSam[ N + i] end
relErrSam = maximum(abs.(( degGroupsSam - degGroups )./degGroups )[degGroups.!=0])

## Test of estimates efficiency on the data

## Repeat on piero's DATA
using MAT
pieroEmidEst = matread("/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/matlabFiles/data\ Piero/forecastEMID/fitness_timeseries.mat")
pieroEmidData = matread("/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/matlabFiles/data\ Piero/forecastEMID/eMID_weeklyaggregated_binaryDIRECT_postLTRO.mat")

pieroA_T = Bool.(pieroEmidData["A"])

# load my data
halfPeriod = false
## Load dataj
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)

matA_T = AeMidWeekly_T# pieroA_T
degsIO_T = [sumSq(matA_T,2);sumSq(matA_T,1)]
T = size(matA_T)[3]
allFitSS =  StaticNets.estimate( StaticNets.SnapSeqNetDirBin1(degsIO_T); identPost = false,identIter= true )
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

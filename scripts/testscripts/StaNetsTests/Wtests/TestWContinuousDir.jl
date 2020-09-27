N = 10
NG = N#5


maxW = 2000
minW = 10
tmpStrIO = round.([ Vector(linspace(minW,maxW,N));Vector(linspace(minW,maxW ,N))])
N, StrI,StrO = splitVec(tmpStrIO)
sum(StrI) == sum(StrO)?():error()
Amat = trues(N,N);for n=1:N Amat[n,n] =false end
Amat[StrI.==0,:] = false
Amat[:,StrO.==0] = false
dgpStrIO  =   DegSeq2graphDegSeq(fooNetModelDirW1Afixed,tmpStrIO;Amat=Amat)
N, StrI,StrO = splitVec(dgpStrIO)
sum(StrI) == sum(StrO)?():error()
Amat = trues(N,N);for n=1:N Amat[n,n] =false end
Amat[StrI.==0,:] = false
Amat[:,StrO.==0] = false



Mod =  NetModelDirW1Afixed(dgpStrIO,Amat)
S = 0.5*sum(dgpStrIO,1)[1]
Par = ones(2N) * 0.5 /(S/(N^2-N))
expMatTest=expMatrix(Mod,Par,Amat)
firstOrderCond(Mod; strIO =dgpStrIO, parNodesIO = Par,Amat=Amat)

(dgpStrIO .== 0 )
##
estPar,estIt,estMod = estimate(Mod)
StrO .==0
degsO = sumSq(Amat,1)
(degsO.==0)
## primi test
uBndPar = bndPar2uBndPar(Mod,estPar)
sum(abs.(estPar - uBndPar2bndPar(Mod,uBndPar)))


##
using Plots
plotly()
plot(estPar[abs.(estPar).<100])


##
N,parI,parO = splitVec(estPar)
minVal = minimum(parI)
parI = parI - minVal
parO = parO + minVal

sum(parI) - sum(parO)

(sum(parI.<0) + sum(parO.<0)) == 0 ? ():error()


##
expMatEst= expMatrix(Mod,estPar)

sam = sampl(estMod,10000)

tmp = sumSq(sam[N,:,:] .- expMatEst[N,:],1)./sqrt(sumSq( expMatEst[N,:],1).^2)
histogram( tmp )
meanMatSam = squeeze(mean(sam,3),3)
strSam = expValStatsFromMat(estMod,meanMatSam)
strGroupsSam = zeros(2NG) ; for i =1:N strGroupsSam[groupsInds[i]] += strSam[i]; strGroupsSam[NG + groupsInds[i]] += strSam[ N + i] end
relErrSam = maximum(abs.(( strGroupsSam - strGroups )./strGroups )[strGroups.!=0])
absErrZeros = maximum(abs.(( strGroupsSam - strGroups ) )[strGroups.==0])

## Test sui dati
using JLD
 load_fold =   "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
# file_name = "Weekly_eMid_Estimates.jld"
# load_path = load_fold*file_name#
# @load(load_path,estSS,estGas,estSS_1GW,estGas_1GW,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)
#
 N,~,T = size(YeMidWeekly_T)
# ##
# estPar_T = zeros(T,2N)
#
# for t = 1:T
#     Y = YeMidWeekly_T[:,:,t]
#     strIOeMid_t = expValStatsFromMat(fooNetModelDirW1,Y)
#     eMidMod = NetModelDirW1(strIOeMid_t)
#
#     estPar,estIt,estMod = estimate(eMidMod)
#     estPar_T[t,:] = estPar
#     N,parI,parO = splitVec(estPar)
#     minVal = minimum(parI)
#     parI = parI - minVal
#     parO = parO + minVal
#
#     println(t)
#      (sum(parI.<0) + sum(parO.<0)) == 0 ? ():error()
# end
#
# estContW_SS_T = estPar_T
# file_name = "Weekly_eMid_Estimates_ContinuousWeights.jld"
# load_path = load_fold*file_name#
# @save(load_path,estContW_SS_T,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)

file_name = "Weekly_eMid_Estimates_ContinuousWeights.jld"
load_path = load_fold*file_name#
@load(load_path,estContW_SS_T,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)
 N,~,T = size(YeMidWeekly_T)

##
estParPos_T=zeros(estContW_SS_T)
for t=1:T
    estPar = estContW_SS_T[t,:]
    N,parI,parO = splitVec(estPar)
    minVal = minimum(parI)
    parI = parI - minVal
    parO = parO + minVal
    estParPos_T[t,:] = [parI;parO]
end

estParPos_T


##
t=1
par = estParPos_T[t,:]
N,parI,parO = splitVec(par)

minimum(parI)
minimum(parO)
sum(parI) - sum(parO)

##
using Plots
plotly()
pdegsI = plot(degsIO_T[:,1:N],title = "  $(round.(squeeze(mean(degsIO_T[:,1:N],1),1),1))")
pparO = plot(estParPos_T[:,N+1:2N])#,title = round.(squeeze(mean(Fitness_T[:,N+1:2N],1),1),1))  $(round.(UMdgp[N+1:2N],1))")
pdegsO = plot(degsIO_T[:,N+1:2N],title =  "  $(round.(squeeze(mean(degsIO_T[:,N+1:2N],1),1),1))")

plot(pdegsI,pdegsO,layout = (2,1),legend=:none,size = (1200,600))


## test IPF function
YallT = meanSq(YeMidWeekly_T,3)
strIOallT = expValStatsFromMat(fooNetModelDirW1,YallT)
#zeros on the diagonal
A =  ones(N,N); for i=1:N A[i,i] = 0 end  #Int.(YallT .> 0) #
nnzInds = strIO.!=0
~,nnzIndsI,nnzIndsO = splitVec(nnzInds)

@time  ipfMat,~ = estimateIPFMat(fooNetModelDirW1,strIOallT,A)
eMidMod = NetModelDirW1(strIOallT)
@time  estPar,estIt,estMod = estimate(eMidMod)
gammaMat = expMatrix(estMod,estPar)

gammaMat = sortMat(gammaMat )
ipfMat = sortMat(ipfMat)
gammaMatA = (gammaMat.>0)
ipfMatA = ipfMat.>0
sum(gammaMatA) - sum(ipfMatA)
sum(.!ipfMatA.& .!gammaMatA)
plotMat = log.(gammaMat./ipfMat)
 contour(plotMat,fill=true,yflip = true)
##

contour(sortMat(log.(gammaMat)),fill=true,yflip = true)
function sortMat(M::Array{<:Real,2})
    strI = sum(M,2)
    N = length(strI)
    M1 = [strI Vector(1:N) M ]

    tmp = sortrows(M1,rev=true)
    inds = Int.(tmp[:,2])
    Msorted = tmp[:,3:end]
    #M2 = [strI Msorted' ]
    #Msorted = sortrows(M2,rev=true)[:,2:end]'
    Msorted = Msorted[:,inds]

    return Msorted
end

sortMat(gammaMat[1:4,1:4,])
##






#






















#

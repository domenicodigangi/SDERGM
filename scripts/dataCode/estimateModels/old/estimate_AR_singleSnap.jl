using HelperFunDom,StaNets,JLD
using PyCall; pygui(:qt); using PyPlot


## Load data
#whole period data file
loadFileName = "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/Weekly_eMid_Data_from_2009_06_22_to_2015_02_27.jld"
@load(loadFileName, AeMidWeekly_Active_T,banksIDs,inactiveBanks, YeMidWeekly_Active_T,weekInd,datesONeMid ,degsIO_T,strIO_T)


## Methods definitions

struct SnapSeqNetDirBin1
    obsT:: Array{<:Real,2}  # In degrees and out Degrees for each t . The last
                            # dimension is always time. 2NxT
    Par::Array{<:Real,2} # In and Out Parameters [InW_T;outW_T] for each t
end
SnapSeqNetDirBin1(obs::Array{<:Real,2}) = SnapSeqNetDirBin1(obs,zeros(obs))


function estimate(Model::SnapSeqNetDirBin1; degsIO_T::Array{<:Real,2}=Model.obsT,targetErr::Real=1e-12)
   N2,T = size(degsIO_T);N = round(Int,N2/2)
   UnParT = zeros(2N,T)
   prog = 0
   for t = 1:T
       degs_t = degsIO_T[:,t]
       UnParT[:,t]  = identify!(fooNetModelDirBin1,estimate(fooNetModelDirBin1; degIO = degs_t , targetErr =  targetErr)[1])
       round(t/T,2)>prog ? (prog=round(t/T,2);println(prog) ):()
   end
 return UnParT
end

function bigNegParFun(Model::NetModelDirBin1,N::Int; smallProb = 1e-1)
    #define a number big enough to play the role of Inf for
    # purposes of sampling N(N-1) bernoulli rvs with prob 1/(1+exp(  bigNumb))
    # prob of sampling degree > 0 (N-1)/(1+exp(bigNumb)) < 1e-6
    bigNumb = - log((N-1)/smallProb - 1)
    return bigNumb
 end

## Application to Data
snapModData = SnapSeqNetDirBin1(degsIO_T)
singSnapEst_T = estimate(snapModData)
N2,T = size(degsIO_T);N = round(Int,N2/2)

## Stima local likelihood, la ottengo facendo i locally weighted degrees
Nalpha = 30 #Number of different alpha values for EWMA
degsIO_EWMA_T_Alpha = zeros(Float64,2N,T,Nalpha)
alphaVals = linspace(0.4+1/Nalpha,1,Nalpha)
# compute the EWMA degrees for different values of alpha
for alphaInd=1:Nalpha
    alpha = alphaVals[alphaInd]
    alphaVec =ones(1) *alpha
    for t=1:T
        wVec = alphaVec./sum(alphaVec)
        alphaVec*=alpha; push!(alphaVec,alpha)
        #println(wVec)
        t>1? degsIO_EWMA_T_Alpha[:,t,alphaInd] = sumSq(wVec'.*degsIO_T[:,1:t],2):()# (degIO_EWMA_T[:,t] = degsIO_T[1,:])
    end
end
# estimate the single snapshots on the EWMA degrees
singSnapEst_T_Alpha = zeros(2N,T,Nalpha)
for n=1:Nalpha
    singSnapEst_T_Alpha[:,:,n] = estimate( SnapSeqNetDirBin1(degsIO_EWMA_T_Alpha[:,:,n]))
    println(n)
end

plot(singSnapEst_T[1:N,:]')
#starting in a uniform point do I always end up is some particular identification?
x = (sumSq(singSnapEst_T[1:N,:],1)-sumSq(singSnapEst_T[1+N:2N,:],1))
y = Float64.(sumSq(degsIO_T[:,:],1))
plot(x,y)
cor(x,y)

##
degIO = degsIO_T[:,1]
parUn,it = estimate(fooNetModelDirBin1; degIO = degIO , targetErr =  1e-15)























#

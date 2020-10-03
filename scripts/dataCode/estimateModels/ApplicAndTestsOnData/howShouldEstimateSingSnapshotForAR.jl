using HelperFunDom,StaticNets,JLD,MLBase
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
function estimate(Model::SnapSeqNetDirBin1; degsIO_T::Array{<:Real,2}=Model.obsT,
                        targetErr::Real=1e-12,identPost=false,identIter=false,zeroDegFit="small")
   N2,T = size(degsIO_T);N = round(Int,N2/2)
   UnParT = zeros(2N,T)
   prog = 0
   for t = 1:T
       degs_t = degsIO_T[:,t]
       UnPar_t , Niter  = estimate(fooNetModelDirBin1; degIO = degs_t , targetErr =  targetErr,identIter=identIter)
       identPost?UnParT[:,t] = identify!(fooNetModelDirBin1,UnPar_t): UnParT[:,t] = UnPar_t
       round(t/T,2)>prog ? (prog=round(t/T,2);println((prog,Niter)) ):()
   end
   if uppercase(zeroDegFit) == "INF"
       # zero degrees have -Inf fitnesses
   elseif uppercase(zeroDegFit) == "10AS21"
       # such that the difference  fit(2) - fit(1) == fit(1)-fit(0) , where
       # fit(n) means the fitness associated with a node of deg n.
       # when more than one choice is possible (the fitness depends on both degrees)
       # than choose fit(2) and fit(1) that maximize the difference

   elseif uppercase(zeroDegFit) == "SMALL"
       # such that the probability of observing a link is very small
       UnParT[degsIO_T.==0] = -5
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
singSnapEstNoIdent_T = estimate(snapModData)
singSnapEstIdent_T = estimate(snapModData; identPost = true,identIter=false)
singSnapEstIdentIter_T = estimate(snapModData; identPost = false,identIter=true)
N2,T = size(degsIO_T);N = round(Int,N2/2)

# plot(singSnapEstNoIdent_T[1:30,:]')
# plot(singSnapEstIdent_T[1:3,:]')
# plot(singSnapEstIdentIter_T[1:3,:]')
## check that equal degrees are associated to equal fitnesses.
using StatsBase
singSnapEst_T = singSnapEstIdent_T
fit,degs = [singSnapEst_T[1:N,1] singSnapEst_T[1+N:2N,1]],[degsIO_T[1:N,1] degsIO_T[1+N:2N,1]]
[unique(degs,1) unique(fit,1)]
# actually is a pair (InDeg,OutDeg) that is associated with a pair (InFit,OutFit). That is because
# there cannot be single loops. Hence the value of a fitness has to take into account
# that links cannot come from the node itself. So all things equal a node with
# (1,1) should have an higher  InFit than a (1,0) because it needs to attract more
# from the other ones
##
D = maximum(degsIO_T)
storeFitDistrPerDeg =  Matrix{Vector{Float64}}(D,T)
meansFit = Matrix{Float64}(D,T)
stdFit = Matrix{Float64}(D,T)
testSingSnapEst_T = singSnapEst_T
t=1;d=1
for t=1:T
    degsIOt = degsIO_T[:,t]
    degsUniIOt = sort(unique(degsIOt))
    degsUniIOt = degsUniIOt[2:end]
    println(length(degsUniIOt))
    for d=1:length(degsUniIOt)
        deg = degsUniIOt[d]
        indsDeg = degsIOt[1:N].== deg
        storeFitDistrPerDeg[d,t] = testSingSnapEst_T[1:N,t][indsDeg]
        meansFit[d,t] = mean(testSingSnapEst_T[1:N,t][indsDeg])
        stdFit[d,t] = std(testSingSnapEst_T[1:N,t][indsDeg])
    end

end

deg = 1:4
plot(meansFit[deg,:]',"-");plot(meansFit[deg,:]' + stdFit[deg,:]',"--");plot(meansFit[deg,:]' - stdFit[deg,:]',"--")
tmp = 1./sumSq(degsIO_T,1)
cor(meansFit[2,:],tmp)
mean(meansFit[4,:])
#plot(tmp)
# the mean value of the fitnesses associated to a given degree is highly correlated
# with the inverse of total number of links


## lets estimate some AR(1)

using AReg

function estManyAR(obs::Array{<:Real,2};p::Int=1)
    #obs is a matrix with a time series for each row (time is last dimension)
    D,T = size(obs)
    isConst = ([length(unique(obs[i,:])) for i=1:D ] .==1)
    parEst = zeros(D,p+1)
    arp = AReg.ARp(ones(p+1),.1,[])
    for d=1:D
        if ~isConst[d]
            parEst[d,:] = AReg.fitARp(obs[d,:],arp).estimates
        end
        #println(d)
    end
    return parEst,isConst
end
function oneStepForAR1(obs_t::Float64,AR1Par::Vector{<:Real})
    #given w and β AR1 parameters, compute the one step ahead forecast from obs_t
    obstp1 = AR1Par[1] + AR1Par[2] * obs_t
    return obstp1
end
#Estimate AR1 on the training sample
function forecastEvalAR1(allFit::Array{Float64,2})
    N2,T = size(allFit)
    Ttrain = round(Int,T/2)
    Ttest = T-Ttrain
    trainFit = allFit[:,1:Ttrain]
    testFit = allFit[:,1:Ttest]
    estAR1,isCon = estManyAR(trainFit)
    ## One step ahead forecasts for each fitness for each time in Test sample
    foreFit = zeros(testFit[:,1:end])
    for n=1:2N,t=1:Ttest-1
        foreFit[n,t+1] = oneStepForAR1(testFit[n,t],estAR1[n,:])
    end

    foreVals = 100 * ones(N*(N-1)*(Ttest-1))
    realVals = trues(N*(N-1)*(Ttest-1))
    noDiagInd = .!Bool.(eye(N))
    lastInd=1
    NmaxLinks = N*(N-1)
    Nobs = NmaxLinks*T
    for t=2:Ttest
        adjMat =  AeMidWeekly_Active_T[:,:,Ttrain+1:end][:,:,t]
        expMat = expMatrix(fooNetModelDirBin1,foreFit[:,t])
        foreVals[lastInd:lastInd+NmaxLinks-1] = expMat[noDiagInd]
        realVals[lastInd:lastInd+NmaxLinks-1] = adjMat[noDiagInd]

        lastInd += NmaxLinks
    end
    return realVals,foreVals
end
##

forecastEvalAR1(allFit)
function rocCurve(realVals::BitArray{1},foreVals::Vector{Float64}; Th::Vector{Float64}=Vector(0:0.00001:1))
    rocOut = roc(realVals,foreVals,Th)
    Nth = length(Th)
    tpr = [rocOut[i].tp/rocOut[i].p for i=1:Nth]
    fpr = [rocOut[i].fp/rocOut[i].n for i=1:Nth]
    auc = sum(-diff(fpr).*tpr[1:end-1])
    plot(fpr,tpr)
    xlim([0,1])
    ylim([0,1])
    title("AUC =  $(round(auc,3))")
    grid()
    return tpr, fpr,auc
end
rocCurve(allFit::Array{Float64,2}) = (tmp = forecastEvalAR1(allFit); rocCurve(tmp[1],tmp[2]))
rocCurve(realVals,foreVals)


rocCurve(singSnapEstIdent_T)
rocCurve(singSnapEstNoIdent_T)
rocCurve(singSnapEstIdentIter_T)






















#

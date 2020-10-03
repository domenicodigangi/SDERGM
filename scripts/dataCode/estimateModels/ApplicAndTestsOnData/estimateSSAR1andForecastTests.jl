using HelperFunDom,AReg,StaticNets,JLD,MLBase,StatsBase
using PyCall; pygui(:qt); using PyPlot


## Load data
#half period
loadFileName = "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/Weekly_eMid_Data_from_2012_03_12_to_2015_02_27.jld"
#whole period
#loadFileName = "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/Weekly_eMid_Data_from_2009_06_22_to_2015_02_27.jld"
@load(loadFileName, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)


## Methods definitions

## Application to Data
N2,T = size(degsIO_T);N = round(Int,N2/2)
#estimate different versions
snapModData = SnapSeqNetDirBin1(degsIO_T)
singSnapEstNoIdent_T = estimate(snapModData)
singSnapEstIdent_T = estimate(snapModData; identPost = true,identIter=false)
singSnapEstIdentIter_T = estimate(snapModData; identPost = false,identIter=true)
# plot(singSnapEstNoIdent_T[1:30,:]')
# plot(singSnapEstIdent_T[1:3,:]')
# plot(singSnapEstIdentIter_T[1:3,:]')

## check that equal degrees are associated to equal fitnesses.
    # singSnapEst_T = singSnapEstIdent_T
    # fit,degs = [singSnapEst_T[1:N,1] singSnapEst_T[1+N:2N,1]],[degsIO_T[1:N,1] degsIO_T[1+N:2N,1]]
    # [unique(degs,1) unique(fit,1)]
    # # actually is a pair (InDeg,OutDeg) that is associated with a pair (InFit,OutFit). That is because
    # # there cannot be single loops. Hence the value of a fitness has to take into account
    # # that links cannot come from the node itself. So all things equal a node with
    # # (1,1) should have an higher  InFit than a (1,0) because it needs to attract more
    # # from the other ones
    # ##
    # D = maximum(degsIO_T)
    # storeFitDistrPerDeg =  Matrix{Vector{Float64}}(D,T)
    # meansFit = Matrix{Float64}(D,T)
    # stdFit = Matrix{Float64}(D,T)
    # testSingSnapEst_T = singSnapEst_T
    # t=1;d=1
    # for t=1:T
    #     degsIOt = degsIO_T[:,t]
    #     degsUniIOt = sort(unique(degsIOt))
    #     degsUniIOt = degsUniIOt[2:end]
    #     println(length(degsUniIOt))
    #     for d=1:length(degsUniIOt)
    #         deg = degsUniIOt[d]
    #         indsDeg = degsIOt[1:N].== deg
    #         storeFitDistrPerDeg[d,t] = testSingSnapEst_T[1:N,t][indsDeg]
    #         meansFit[d,t] = mean(testSingSnapEst_T[1:N,t][indsDeg])
    #         stdFit[d,t] = std(testSingSnapEst_T[1:N,t][indsDeg])
    #     end
    #
    # end
    #
    # deg = 1:4
    # plot(meansFit[deg,:]',"-");plot(meansFit[deg,:]' + stdFit[deg,:]',"--");plot(meansFit[deg,:]' - stdFit[deg,:]',"--")
    # tmp = 1./sumSq(degsIO_T,1)
    # cor(meansFit[2,:],tmp)
    # mean(meansFit[4,:])
    # #plot(tmp)
    # # the mean value of the fitnesses associated to a given degree is highly correlated
    # # with the inverse of total number of links


## lets estimate some AR(1)
forecastEvalAR1Net(singSnapEstIdent_T,AeMidWeekly_T)
#Compute rocCurve for each method
rocCurve(singSnapEstIdent_T,AeMidWeekly_T)
rocCurve(singSnapEstNoIdent_T,AeMidWeekly_T)
rocCurve(singSnapEstIdentIter_T,AeMidWeekly_T)

## Test that the naive model is Biased in the estimate of B (as Piero mentioned)
#estimate one AR1 for each fitness
using AReg,RCall
testSSFitSeq = singSnapEstIdent_T
(estAR,isCon) = estManyAR(testSSFitSeq)
# check results with R estimates
testSSFitSeqNnc =  testSSFitSeq[.!isCon,:]
# run estimates in R
@rput testSSFitSeqNnc
R"estAR_Rnnc = matrix(,nrow = dim(testSSFitSeqNnc)[1],ncol = 2)"
R"estARvar_Rnnc = matrix(,nrow = dim(testSSFitSeqNnc)[1],ncol = 1)"
R"for(n in 1:dim(testSSFitSeqNnc)[1]){ x = ar(testSSFitSeqNnc[n,],order.max = 1,demean = TRUE); estAR_Rnnc[n,] = c(x$x.mean,x$ar); estARvar_Rnnc[n] = c(x$var.pred) }"
@rget estAR_Rnnc estARvar_Rnnc
estAR_Rnnc[:,1] = estAR_Rnnc[:,1].*(1-estAR_Rnnc[:,2]) # the R function used computes the mean i use w
estAR_R = zeros(estAR); estAR_R[.!isCon,:] = estAR_Rnnc
estARvar_R = zeros(N2); estARvar_R[.!isCon,:] = estARvar_Rnnc

diff = (estAR_R .-estAR)[.!isCon,:]
plot(diff[:,1])
#plt[:hist](diff)[:,2] # Histogram
plot(estAR[:,2])
indRprob = estAR_R[:,2].<0
plot(degsIO_T[indRprob,:])
sum(degsIO_T[indRprob,:].>0,2)
estAR_R[indRprob,2] = 0
estAR[indRprob,2] = 0
estAR[indRprob,1] = -25
indJprob = estAR[:,2].<0
(sum(indRprob), sum(indJprob))
estAR[indRprob & (.!indJprob),:]
rocCurve(testSSFitSeq;estAndInd = (estAR_R,isCon))
rocCurve(testSSFitSeq;estAndInd = (estAR,isCon))
# My AR estimates (OLS) are in line with the Yule Walker ones from R. There are a few cases that
# give problems to both. due to almost constant degrees and, as consequence, fit
# some persistency parameters are estimated negative. I choose to put such cases
# to B = 0  for the moment, because this helps the out of sample performances of
# the naive model (keep calling it naive but we'll see who wins...)

## Simulate and estimate the AR1 on single snapshots to test for bias in the estimates
testSSFitSeq = singSnapEstIdent_T
D,T = size(testSSFitSeq)
i=1
using AReg
#estimate simulated data
testSSFitReEst =  estimate(snapModData; identPost = true,identIter=false)


#
function simulFit(Model::SnapSeqNetDirBin1,dgpARparAndVar::Array{<:Real,2},T::Int)
    N2 = length(dgpARparAndVar[:,1])
    simFitT = zeros(N2,T)
    for i=1:N2 simFitT[i,:]= AReg.simulateAR1(dgpARparAndVar[i,1:2],dgpARparAndVar[i,3];T=T) end
    return simFitT
 end

reEstDgpARparAndVar = [estAR estARvar_R]; reEstDgpARparAndVar[:,3] = 1
Ntest = 10
storeReEstAR = zeros(N2,3,Ntest)
for j=1:Ntest
    #use as AR parameters those estimated on the data
    #simulate AR paths for fit
    fitSimulT = simulFit(fooSnapSeqNetDirBin1,reEstDgpARparAndVar,T)
    #plot(fitSimulT')
    #given fit_T simulate one network
    sampledA_T = sampl(fooSnapSeqNetDirBin1,fitSimulT)
    samplDegIO_T = [sumSq(sampledA_T,2); sumSq(sampledA_T,1)]
    singSnapReEstIdent_T = estimate(fooSnapSeqNetDirBin1;degsIO_T = samplDegIO_T,  identPost = true,identIter=true)

    #estimate(fooNetModelDirBin1;degIO = samplDegIO_T[:,47],targetErr = 1e-4)#,identIter=identIter)
    #estimate the B
    (reEstAR,isCon) = estManyAR(singSnapReEstIdent_T)

    storeReEstAR[:,:,j] = [reEstAR isCon]
end
##

plt[:hist](diff) # Histogram
plot(diff)
diffMat = storeReEstAR[:,2,:] .-reEstDgpARparAndVar[:,2]

plot(reEstDgpARparAndVar[:,2],mean(diffMat,2),".")



# check Bias in my AR! estimates. There is not!!
Nsim = 10000
    arp = AReg.ARp(ones(1+1),.1,[])
    tmp = zeros(Nsim,2)
    parSim = [-1.5,0.85]
    for t=1:Nsim
        tmp[t,:] = AReg.fitARp(AReg.simulateAR1(parSim,20.0;T=T),arp).estimates .- parSim
    end
    plt[:hist](tmp[:,2],50)

















#

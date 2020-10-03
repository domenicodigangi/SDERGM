using JLD, Plots,DynNets,HelperFunDom,StatPlots,StatsBase
plotly()

load_fold =   "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
file_name = "Weekly_eMid_Estimates.jld"
load_path = load_fold*file_name#
@load(load_path,estSS,estGas,estSS_1GW,estGas_1GW,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)

##
N,degsI_T,degsO_T = splitMat(degsIO_T)
pdegsI = plot(degsI_T,title = "  $(round.(squeeze(mean(degsIO_T[:,1:N],1),1),1))")
pdegsO = plot(degsO_T,title =  "  $(round.(squeeze(mean(degsIO_T[:,N+1:2N],1),1),1))")
plot(pdegsI,pdegsO,layout = (2,1),legend=:none,size = (1200,600))


##

function plotFilter(filtPar_T,textTitle)
    ~,filtParI_T,filtParO_T = splitMat(filtPar_T)
    indzI = sumSq(degsI_T,1).==0
    indzO = sumSq(degsO_T,1).==0
    pparI = plot(filtParI_T[:,.!indzI],title = textTitle)

    pparO =  plot(filtParO_T[:,.!indzO])

    #plot(pdegsI,pdegsO,layout = (2,1),legend=:none,size = (1200,600))
    plot(pparI,pparO,layout = (2,1),legend=:none,size = (1200,600))
end


function plotSingSnapSeq(filtPar_T)
    ~,filtParI_T,filtParO_T = splitMat(filtPar_T)
    filtParI_T[degsI_T.==N] = maximum(filtParI_T[degsI_T.!=N])
    filtParI_T[degsI_T.==0] = minimum(filtParI_T[degsI_T.!=0])
    filtParO_T[degsO_T.==N] = maximum(filtParO_T[degsO_T.!=N])
    filtParO_T[degsO_T.==0] = minimum(filtParO_T[degsO_T.!=0])


    indzI = sumSq(degsI_T,1).==0
    indzO = sumSq(degsO_T,1).==0
    pparI = plot(filtParI_T[:,.!indzI])
    pdegsI = plot(degsI_T,title = "  $(round.(squeeze(mean(degsIO_T[:,1:N],1),1),1))")
    pparO =  plot(filtParO_T[:,.!indzO])
    pdegsO = plot(degsO_T,title =  "  $(round.(squeeze(mean(degsIO_T[:,N+1:2N],1),1),1))")
    #plot(pdegsI,pdegsO,layout = (2,1),legend=:none,size = (1200,600))
    plot(pparI,pparO,layout = (2,1),legend=:none,size = (1200,600))
end

## Filter and Plot
estModNpar = DynNets.GasNetModelDirBin1(degsIO_T,estGas[1],[Vector(1:N),ones(Int,N),ones(Int,N)])
filterParNpar = DynNets.gasFilter(estModNpar)[1]
plotFilter(filterParNpar, "N W Par")

estMod1par = DynNets.GasNetModelDirBin1(degsIO_T,estGas_1GW[1],[ones(Int,N),ones(Int,1),ones(Int,1)])
filterPar1par = DynNets.gasFilter(estMod1par)[1]
plotFilter(filterPar1par,"1 W Par")
filtPar_T,like = DynNets.gasFilter(estMod1par)
plotSingSnapSeq(estSS)


## Compare filtered paths

bias1_N = meanSq(filterParNpar .- filterPar1par,1)
squareErr1_N = sqrt(meanSq((filterParNpar .- filterPar1par).^2,1))

plotly()
scatter(meanSq(degsIO_T,1), log.(squareErr1_N) )


## Compare link probabilities
T = length(degsIO_T[:,1])
expMat(parIO) = StaticNets.expMatrix(DynNets.StaModType(estModNpar),parIO)
expM_1_T = zeros(T,N,N)
expM_N_T = zeros(T,N,N)
expM_SS_T = zeros(T,N,N)
for t = 1:T
    expM_1_T[t,:,:] = expMat(filterPar1par[t,:])
    expM_N_T[t,:,:] = expMat(filterParNpar[t,:])
    expM_SS_T[t,:,:] = expMat(estSS[t,:])
end
expLinks1_T = sumSq(sumSq(expM_1_T,2),2)
expLinksN_T = sumSq(sumSq(expM_N_T,2),2)
expLinksSS_T = sumSq(sumSq(expM_SS_T,2),2)
# a matrix that ensures the same total number of links bu uniformly distributed
expMrand = ones(expM_N_T)
for t = 1:T expMrand[t,:,:] = expLinksN_T[t]/(N^2-N) .*putZeroDiag(expMrand[t,:,:]) end
sumSq(sumSq(expMrand,2),2)

klDist(A::Array{<:Real,2},B::Array{<:Real,2}) = sum( (A.*log.(A./B))[(B.!=0).&(A.!=0)])
eucDist(A::Array{<:Real,2},B::Array{<:Real,2}) = sqrt(mean( (A.-B ).^2))
distancesKL = zeros(T,3)
distancesE = zeros(T,3)
for t=1:T
    distancesKL[t,:] = [   klDist(expM_N_T[t,:,:], expM_1_T[t,:,:]);
                    klDist(expM_N_T[t,:,:], expM_SS_T[t,:,:]);
                    klDist(expM_N_T[t,:,:], expMrand[t,:,:])]
    distancesE[t,:] = [   eucDist(expM_N_T[t,:,:], expM_1_T[t,:,:]);
                    eucDist(expM_N_T[t,:,:], expM_SS_T[t,:,:]);
                    eucDist(expM_N_T[t,:,:], expMrand[t,:,:])]
end
plTotLinks = plot([expLinks1_T expLinksN_T expLinksSS_T],label = ["1 W Par" "N W Par" "Sing Snap"],title=" Expected Number of Total Links  ",size = (1200,600))
pldistance = plot(distancesE,label = ["1 W Par"  "Sing Snap" "Rand Mat"],title="Distance from   Filtered Mean Matrix With N W Par ",size = (1200,600))
##



N= 1000
@time x = 0.3.*ones(N,N)




##






##







##

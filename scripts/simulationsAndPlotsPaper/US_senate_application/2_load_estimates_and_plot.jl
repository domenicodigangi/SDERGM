
using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"
DrWatson.greet()

using Test
using ScoreDrivenERGM
import ScoreDrivenERGM:StaticNets,DynNets, Utilities
using LinearAlgebra
using SharedArrays
using PyPlot
using Statistics
using StatsBase
using DataStructures
using JLD2
using RCall

using PyCall



model = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges +  gwesp(decay = 0.25,fixed = TRUE)", true), indTvPar = [false, true], scoreScalingType="FISH_D")

JLD2.@load( datadir("US_congr_covoting", "estimates") * "\\SS_est_MCMC.jld2",estParSS_T)#, obsMat_T, changeStats_T, stats_T)

JLD2.@load( datadir("US_congr_covoting", "estimates")*"\\SD_est_$(DynNets.name(model))", res_conf, res_est, estParStatic)




#%% Adapt plot to fit paper format
begin 
fig = figure(figsize=(18,7))
names = ["Number of Links" "GWESP"]
labSize = 26
tickSize = 20
legSize = 20
SDFiltPar_T = res_est.fVecT_filt
confBands1 = res_conf.confBandsFiltPar
T = size(SDFiltPar_T)[2]
Nterms = DynNets.number_ergm_par(model)
for parInd = 1:Nterms
    x = collect(1:T).+40
    subplot(1,Nterms,parInd)
    plot(x, SDFiltPar_T[parInd,:],"r")
    plot(x, estParSS_T[parInd,:],".b")
    plot(x, ones(T).*estParStatic[parInd],"--b")

    b= 1
    nBands=1
    fill_between(x, confBands1[parInd, :, b,1], y2 =confBands1[parInd,:,b, 2],color =(0.9, 0.2 , 0.2, 0.1), alpha = 0.2*b/nBands  )#, color='b', alpha=.1)
    
    title(names[parInd],fontsize = labSize+4)
    ylabel("\$ \\theta_$(parInd) \$",fontsize = labSize)
    xlabel("Congress Number",fontsize = labSize)
    xticks(fontsize = tickSize)
    yticks(fontsize = tickSize)
    legTex = [ "SD-ERGM "; "Cross Sect. ERGM"; "Constant ERGM" ]
    legend(legTex,fontsize = legSize,framealpha = 0.95)
end
    
tight_layout()

end
##

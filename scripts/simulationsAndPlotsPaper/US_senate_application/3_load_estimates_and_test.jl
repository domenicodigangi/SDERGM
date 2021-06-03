#########
#Created Date: Tuesday June 1st 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Thursday June 3rd 2021 5:18:30 pm
#Modified By:  Domenico Di Gangi
#-----
#Description:
#-----
########



using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"
DrWatson.greet()

using ScoreDrivenERGM
import ScoreDrivenERGM:StaticNets,DynNets, Utilities
using PyPlot
using StatsBase
using JLD2

using RCall


R"""library(stats)
    library(aTSA)
    library(forecast)
    """


model = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges +  gwesp(decay = 0.25,fixed = TRUE)", true), indTvPar = [false, true], scoreScalingType="FISH_D")

JLD2.@load( datadir("US_congr_covoting", "estimates") * "\\SS_est_MCMC.jld2",estParSS_T)#, obsMat_T, changeStats_T, stats_T)

JLD2.@load( datadir("US_congr_covoting", "estimates")*"\\SD_est_$(DynNets.name(model))", res_conf, res_est, estParStatic)

SDFiltPar_T = res_est.fVecT_filt
confBands1 = res_conf.confBandsFiltPar

@rput SDFiltPar_T





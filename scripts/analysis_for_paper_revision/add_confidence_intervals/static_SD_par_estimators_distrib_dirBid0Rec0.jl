"""
What is the distribution of static parameters of the SD filter,  when used as DGP? is it well approximated by the MLE normal?
"""

#region import and models
using ScoreDrivenExponentialRandomGraphs

import ScoreDrivenExponentialRandomGraphs:StaticNets, DynNets

import ScoreDrivenExponentialRandomGraphs.DynNets:GasNetModel,GasNetModelDirBin0Rec0, sample_dgp, statsFromMat, array2VecGasPar, unrestrict_all_par, conf_bands_par_uncertainty, avg_grad_and_hess_obj_SD_filter_time_seq, conf_bands_par_uncertainty, number_ergm_par, filter_and_conf_bands, conf_bands_coverage, estimate, mle_distrib_filtered_par, plot_filtered_and_conf_bands
using ScoreDrivenExponentialRandomGraphs.Utilities

using PyPlot
pygui(true)


using ForwardDiff
using StatsBase
using LinearAlgebra
using Distributions
using Statistics

using JLD

model_mle = DynNets.GasNetModelDirBin0Rec0_mle()
model_pmle = DynNets.GasNetModelDirBin0Rec0_pmle()
indTvPar = trues(2)

#endregion

# work in progress, arrivato qui------------------------------------
DynNets.starting_point_optim(model_mle)


using Statistics
using Distributions

X = rand(MvNormal(zeros(3),Diagonal([1.0,3.0,5.0])), 1000)

cov(X')
var(X[2,:])


#region filter checks 
begin
T=300
N=100
quantilesVals = [[0.975, 0.025]]
parDgpT = DynNets.dgp_misspecified(model_mle, "AR", N, T;  minValAlpha = 0.2, maxValAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false, B =0.95, sigma = 0.001)

# quick visual checks
# DynNets.sample_est_mle_pmle(model_mle, parDgpT, N, 1; plotFlag = true)
A_T_dgp = sample_dgp(model_mle, parDgpT,N)


res_mle = filter_and_conf_bands(model_mle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = parDgpT)

res_pmle = filter_and_conf_bands(model_pmle, A_T_dgp, quantilesVals; parDgpT = parDgpT, plotFlag=true)


end
#endregion

res_mle[8]
res_pmle[8]


res_mle[8]
res_pmle[8]

#region distribution of filtered parameters checks 
begin
T=100
N=50
model = model_pmle
quantilesVals = [[0.975, 0.025] ]

#parDgpT = DynNets.dgp_misspecified(model_mle, "SIN", N, T;  minValAlpha = 0.2, maxValAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false, B =0.990, sigma = 0.005)
# quick visual checks
# DynNets.sample_est_mle_pmle(model_mle, parDgpT, N, 1; plotFlag = true)
#A_T_dgp = sample_dgp(model_mle, parDgpT,N)
obsT = [statsFromMat(model, A_T_dgp[:,:,t]) for t in 1:T ]

estSdResPar, conv_flag, UM_mple, ftot_0 = estimate(model, obsT;indTvPar=indTvPar, indTargPar=falses(2))

vEstSdResPar = array2VecGasPar(model, estSdResPar, indTvPar)

fVecT_filt , target_fun_val_T, sVecT_filt = DynNets.score_driven_filter(model,  vEstSdResPar, indTvPar; obsT = obsT, ftot_0 = ftot_0)


distribFilteredSD, filtCovHatSample, errFlag = mle_distrib_filtered_par(model, obsT, indTvPar, ftot_0, vEstSdResPar)



confBandsParFilt, confBandsParGauss = DynNets.conf_bands_buccheri(model, obsT, indTvPar, fVecT_filt, distribFilteredSD, filtCovHatSample, quantilesVals)

confBandsParSimul = DynNets.conf_bands_par_uncertainty_blasques(model, obsT, fVecT_filt, distribFilteredSD, quantilesVals)




plot_filtered_and_conf_bands(model, N, fVecT_filt, confBandsParFilt;  parDgpT=parDgpT, confBands2=confBandsParSimul)



deltas = distribFilteredSD[:,1,end-10,] 

fig1, ax = subplots(2,1)
for p in 1:2
    x = 1:T
  
    ax[p].hist(deltas)
  
    ax[p].grid()  
end

end


#endregion

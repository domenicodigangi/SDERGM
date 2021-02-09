"""
What is the distribution of static parameters of the SD filter,  when used as DGP? is it well approximated by the MLE normal?
"""

#region import and models
begin

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

end
#endregion

#region  
begin
T=100
N=100
nSample = 100
model = model_pmle
alpha = 0.2
beta = mean(DynNets.beta_min_max_from_alpha_min(0.1, N))
UM = [alpha, beta]
vUnPar, ~ = DynNets.starting_point_optim(model, indTvPar, UM)
vResParDgp = DynNets.restrict_all_par(model, indTvPar, vUnPar)

BDgp = 0.96
ADgp = 0.2
vResParDgp[2:3:end] .= BDgp
vResParDgp[3:3:end] .= ADgp

#start optimization on the correct values
vParUnOptim_0, ARe_min = DynNets.starting_point_optim(model, indTvPar, UM; indTargPar = falses(number_ergm_par(model)))

vParOptim_0 = Real.(DynNets.restrict_all_par(model, indTvPar, vParUnOptim_0))



nStaticSdPar = length(vResParDgp)
vResParEstDistrib = zeros(length(vResParDgp), nSample)
vUnParEstDistrib = zeros(length(vResParDgp), nSample)
HessSum = zeros(nStaticSdPar, nStaticSdPar, nSample)
OpGradSum = zeros(nStaticSdPar, nStaticSdPar, nSample)

quantilesVals = [[0.975, 0.025] ]
end

using Profile
# sample SD dgp
fVecTDgp, A_T_dgp, ~ = DynNets.score_driven_filter( model_mle,N,  vResParDgp, indTvPar; dgpNT = (N,T))
    
res_mle_pmle = DynNets.filter_and_conf_bands(model_pmle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = fVecTDgp, parUncMethod = "WHITE-MLE")



begin
T=300
N=300
nSample = 100
model = model_pmle
alpha = 0.2
beta = mean(DynNets.beta_min_max_from_alpha_min(0.1, N))
UM = [alpha, beta]
vUnPar, ~ = DynNets.starting_point_optim(model, indTvPar, UM)
vResParDgp = DynNets.restrict_all_par(model, indTvPar, vUnPar)

BDgp = 0.96
ADgp = 0.2
vResParDgp[2:3:end] .= BDgp
vResParDgp[3:3:end] .= ADgp

#start optimization on the correct values
vParUnOptim_0, ARe_min = DynNets.starting_point_optim(model, indTvPar, UM; indTargPar = falses(number_ergm_par(model)))

vParOptim_0 = Real.(DynNets.restrict_all_par(model, indTvPar, vParUnOptim_0))



nStaticSdPar = length(vResParDgp)
vResParEstDistrib = zeros(length(vResParDgp), nSample)
vUnParEstDistrib = zeros(length(vResParDgp), nSample)
HessSum = zeros(nStaticSdPar, nStaticSdPar, nSample)
OpGradSum = zeros(nStaticSdPar, nStaticSdPar, nSample)

quantilesVals = [[0.975, 0.025] ]

for n=1:2
    # sample SD dgp
    global fVecTDgp, A_T_dgp, ~ = DynNets.score_driven_filter( model_mle,N,  vResParDgp, indTvPar; dgpNT = (N,T))


    #global res_mle_mle = DynNets.filter_and_conf_bands(model_mle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = fVecTDgp, parUncMethod = "WHITE-MLE")
    
    global res_mle_pmle = DynNets.filter_and_conf_bands(model_pmle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = fVecTDgp, parUncMethod = "WHITE-MLE")

   # global res_mle_boot = DynNets.filter_and_conf_bands(model, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = fVecTDgp, parUncMethod = "PAR-BOOTSTRAP-SAMPLE")

    # global res_mle_boot = DynNets.filter_and_conf_bands(model, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = fVecTDgp, parUncMethod = "PAR-BOOTSTRAP-COVMAT")

end
end



for n=1:nSample
    # sample SD dgp
    fVecTDgp, A_T_dgp, ~ = DynNets.score_driven_filter( model, N, vResParDgp, indTvPar; dgpNT = (N,T))

    
    obsT = [statsFromMat(model, A_T_dgp[:,:,t]) for t in 1:T ]


    arrayAllParHat, conv_flag,UM, ftot_0 = estimate(model, N, obsT; indTvPar=indTvPar, indTargPar=falses(2), vParOptim_0=vParOptim_0)

    
    vResParEstDistrib[:,n] = DynNets.array2VecGasPar(model, arrayAllParHat, indTvPar)

    OpGradSum[:,:,n], HessSum[:,:,n] = DynNets.A0_B0_est_for_white_cov_mat_obj_SD_filter_time_seq(model, N, obsT, vResParEstDistrib[:,n], indTvPar, ftot_0)

    # estCovWhite[:,:,n], errorFlag, OpGradSum[:,:,n], HessSum[:,:,n] = DynNets.white_estimate_cov_mat_static_sd_par(model, obsT, indTvPar, ftot_0, vResParEstDistrib[:,n]; returnAllMats=true)

end

# end
#endregion
n=1
begin
figure()
fVecTDgp, A_T_dgp, ~ = DynNets.score_driven_filter( model, vResParDgp, indTvPar; dgpNT = (N,T))

obsT = [statsFromMat(model, A_T_dgp[:,:,t]) for t in 1:T ]

arrayAllParHat, conv_flag,UM, ftot_0 = estimate(model, obsT; indTvPar=indTvPar, indTargPar=falses(2), vParOptim_0=vParOptim_0)


vResParEstDistrib[:,n] = DynNets.array2VecGasPar(model, arrayAllParHat, indTvPar)

fVecTFilt, A_T_dgp, ~ = DynNets.score_driven_filter( model, vResParEstDistrib[:,n], indTvPar, obsT=obsT, ftot_0=ftot_0)

plot(fVecTFilt[2,:])
plot(fVecTDgp[2,:])
end


begin
indPar = 3
fig, ax = subplots()
ax.hist(vResParEstDistrib[indPar,:], 30)
ylim = ax.get_ylim()
ax.vlines(vResParDgp[indPar], ylim[1], ylim[2], color = "r" )
ax.vlines(vResPar_0[indPar], ylim[1], ylim[2], color = "r" )
#region distribution of filtered parameters checks 
end


#region distribution of white estimators
# 1. C'é un bias di campione finito. Lo stimatore, MLE o White, della covarianza degli stimatori dei parametri é ok? 
begin
estCovWhite = zeros(nStaticSdPar, nStaticSdPar, nSample)
estCovWhitePos = zeros(nStaticSdPar, nStaticSdPar, nSample)
errFlag = falses(nSample)
for n=1:nSample
    opGrad = OpGradSum[:,:,n]
    hess = HessSum[:,:,n]

        parCovHat = pinv(hess) * opGrad * pinv(hess)
        
        parCovHatPosDef, minEigenVal = make_pos_def(parCovHat)

        estCovWhite[:,:,n] = parCovHat
        estCovWhitePos[:,:,n] = parCovHatPosDef
end

parBootCov = cov(vResParEstDistrib')
parNames = ["w_theta", "B_theta", "A_theta", "w_eta", "B_eta", "A_eta"]

fig, ax = subplots(3,2)
for i = 1:6
    j=i
    ax[i].hist(estCovWhite[i,j,:], range=quantile(estCovWhite[i,j,:], [0.01, 0.99]))
    ylim = ax[i].get_ylim()
    ax[i].vlines(parBootCov[i,j], ylim[1], ylim[2], color = "r" )
    ax[i].set_title(parNames[i])
end
suptitle("White Estimators of Diagonal Elements of Covariance Matrix \n Par Bootstrap Estimate in Red. nSample = $nSample ")

tight_layout()
end
#endregion

# Domande
# 1. C'é un bias di campione finito. Lo stimatore, MLE o White, della covarianza degli stimatori dei parametri é ok? 
# 2. Visto che posso costruire la distribuzione dei parametri con parametric bootstrap, perché non usare quella?
# 3. Che relazione c'é tra la distribuzione da parametric bootstrap e la mia idea di usare non parametric bootstrap campionando le osservazioni (che portano in dote la storia fino a quel punto, ma considerando la likelihood della sola osservazione)? 


#endregion

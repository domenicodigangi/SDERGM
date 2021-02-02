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

NON CATTURO LA PARAMETER UNCERTAINTY!!

#region  
begin
T=200
N=50
nSample = 3
model = model_mle
alpha = 0.2
beta = mean(DynNets.beta_min_max_from_alpha_min(0.1, N))
UM = [alpha, beta]
vUnPar, ~ = DynNets.starting_point_optim(model, indTvPar, UM)
vResParDgp = DynNets.restrict_all_par(model, indTvPar, vUnPar)

BDgp = 0.96
ADgp = 0.005
vResParDgp[2:3:end] .= BDgp
vResParDgp[3:3:end] .= ADgp

#start optimization on the correct values
vParOptim_0, ARe_min = DynNets.starting_point_optim(model, indTvPar, UM; indTargPar = falses(number_ergm_par(model)))

vResPar_0 = deepcopy(vResParDgp) 
vResPar_0[3:3:end] .=vResPar_0[3:3:end]./10 
vParOptim_0 = Real.(DynNets.unrestrict_all_par(model, indTvPar, vResPar_0))


nStaticSdPar = length(vResParDgp)
vResParEstDistrib = zeros(length(vResParDgp), nSample)
vUnParEstDistrib = zeros(length(vResParDgp), nSample)
HessSum = zeros(nStaticSdPar, nStaticSdPar, nSample)
OpGradSum = zeros(nStaticSdPar, nStaticSdPar, nSample)

quantilesVals = [[0.975, 0.025] ]

for n=1:2
    # sample SD dgp
    fVecTDgp, A_T_dgp, ~ = DynNets.score_driven_filter( model, vResParDgp, indTvPar; dgpNT = (N,T))


    res_mle = DynNets.filter_and_conf_bands(model, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = fVecTDgp)

end
end
for n=1:nSample
    # sample SD dgp
    fVecTDgp, A_T_dgp, ~ = DynNets.score_driven_filter( model, vResParDgp, indTvPar; dgpNT = (N,T))

    
    obsT = [statsFromMat(model, A_T_dgp[:,:,t]) for t in 1:T ]

    arrayAllParHat, conv_flag,UM, ftot_0 = estimate(model, obsT; indTvPar=indTvPar, indTargPar=falses(2), vParOptim_0=vParOptim_0)

    
    vResParEstDistrib[:,n] = DynNets.array2VecGasPar(model, arrayAllParHat, indTvPar)

    vecUnParAll = DynNets.unrestrict_all_par(model, indTvPar, vResParEstDistrib[:,n])

    vUnParEstDistrib[:,n] = vecUnParAll
    
    OpGradSum[:,:,n], HessSum[:,:,n] = DynNets.A0_B0_est_for_white_cov_mat_obj_SD_filter_time_seq(model, obsT, vecUnParAll, indTvPar, ftot_0)

    # estCovWhite[:,:,n], errorFlag, OpGradSum[:,:,n], HessSum[:,:,n] = DynNets.white_estimate_cov_mat_static_sd_par(model, obsT, indTvPar, ftot_0, vResParEstDistrib[:,n]; returnAllMats=true)

end

# end
#endregion


begin
indPar = 3
fig, ax = subplots()
ax.hist(vResParEstDistrib[indPar,:])
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

for n=1:nSample
    opGrad = OpGradSum[:,:,n]
    hess = HessSum[:,:,n]

    parCovHat = pinv(hess) * opGrad * pinv(hess)
    
    parCovHatPosDef, minEigenVal = make_pos_def(parCovHat)

    estCovWhite[:,:,n] = parCovHat
    estCovWhitePos[:,:,n] = parCovHatPosDef
end

end

parBootCov = cov(vResParEstDistrib')

begin

for i =1:6
    j=i
    figure()
    hist(estCovWhite[i,j,:])
    hist(estCovWhitePos[i,j,:])
    ax = gca()
    ylim = ax.get_ylim()
    vlines(parBootCov[i,j], ylim[1], ylim[2], color = "r" )
    title("$i")
end

end
#endregion

# Domande
# 1. C'é un bias di campione finito. Lo stimatore, MLE o White, della covarianza degli stimatori dei parametri é ok? 
# 2. Visto che posso costruire la distribuzione dei parametri con parametric bootstrap, perché non usare quella?
# 3. Che relazione c'é tra la distribuzione da parametric bootstrap e la mia idea di usare non parametric bootstrap campionando le osservazioni (che portano in dote la storia fino a quel punto, ma considerando la likelihood della sola osservazione)? 


#endregion

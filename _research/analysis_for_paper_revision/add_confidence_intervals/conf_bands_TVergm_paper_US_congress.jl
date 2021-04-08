"""
Estimate a SDERGM on US concress covoting network, filter TV parameters and add confidence intervals
"""

#region import and models
using ScoreDrivenERGM

import ScoreDrivenERGM:StaticNets, DynNets

import ScoreDrivenERGM.DynNets:GasNetModel,GasNetModelDirBin0Rec0, sample_mats_sequence, stats_from_mat, array2VecGasPar, unrestrict_all_par, conf_bands_par_uncertainty, avg_grad_and_hess_obj_SD_filter_time_seq, conf_bands_par_uncertainty, number_ergm_par, estimate_filter_and_conf_bands, conf_bands_coverage, estimate
using ScoreDrivenERGM.Utilities

using PyPlot
pygui(true)


using ForwardDiff
using StatsBase
using LinearAlgebra
using Distributions
using Statistics

using JLD2

model_mle = DynNets.GasNetModelDirBin0Rec0_mle()
model_pmle = DynNets.GasNetModelDirBin0Rec0_pmle()
indTvPar = trues(2)

#endregion

#region quick checks 
begin
T=200
N=300
quantilesVals = [0.975, 0.95, 0.05, 0.025]
parDgpT = DynNets.sample_time_var_par_from_dgp(model_mle, "sin", N, T;  minAlpha = 0.1, maxAlpha = 0.2, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false)
# quick visual checks
# DynNets.sample_est_mle_pmle(model_mle, parDgpT, N, 1; plotFlag = true)
A_T_dgp = sample_mats_sequence(model_mle, parDgpT,N)
res_mle = estimate_filter_and_conf_bands(model_mle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = parDgpT)

res_pmle = estimate_filter_and_conf_bands(model_pmle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = parDgpT)

end


#endregion


using StaticNets, DynNets
using PyPlot
pygui(true)
using JLD2

import Utilities:meanSq, maxLargeVal

@load("./data/congress_covoting_US/juliaEstimates.jld", estParSS_T, obsMat_T,changeStats_T, stats_T)


#estParSS_T = reshape(estParSS_T, 1, 2, 73)

# Estimate SD-PML-ERGM of the paper

indTvPar = BitArray([false, true])
indTargPar =indTvPar #BitArray([true,true])

obsT = changeStats_T #[Real.(changeStats_t) for changeStats_t in changeStats_T]

model = DynNets.GasNetModelDirBinGlobalPseudo(obsT, DynNets.fooGasPar, indTvPar,"")

estPar,convFlag,UM,startPoint,hess_opt = DynNets.estimate(model;UM =meanSq(estParSS_T[:,:],1), indTvPar = indTvPar, indTargPar = indTargPar,hess_opt_flag = true)



targetAllTv = false
 gasParSampled  = Array{Array{Float64,2}}(undef, Nsample)
 flag_cov_issue = falses(Nsample)
 gasParEst= fill(fill(Float64[],2), Nsample)
 startPointEst = fill(fill(0.0,2), Nsample)
 staticEst = fill(fill(0.0,2), Nsample)
 filtPar_T_Nsample = zeros(Nterms,T,Nsample)
 pVals_Nsample = zeros(Nterms,Nsample)
 scoreAutoc_Nsample = zeros(Nterms,Nsample)
 convFlag = falses(Nsample)


T_train = size(obsT)[1]

restrCase = false
if restrCase
    # restricted parameters
    f_T(x) = DynNets.logLike_T(model::DynNets.GasNetModelDirBinGlobalPseudo,
                                obsT, x, model.indTvPar)
    allpar0 = DynNets.array2VecGasPar(model, estPar, model.indTvPar)
    vec_of_f_t =[x -> DynNets.logLike_t(model::DynNets.GasNetModelDirBinGlobalPseudo,
                                        obsT[1:t], x, model.indTvPar)
                 for t in 1:T_train]
else
    #Not restricted parameters
    f_T(x) = DynNets.logLike_T(model::DynNets.GasNetModelDirBinGlobalPseudo,
                                obsT,DynNets.restrictGasPar(model, x, model.indTvPar),
                                model.indTvPar)
    allpar0 = DynNets.unRestrictGasPar(model,
                                        DynNets.array2VecGasPar(model, estPar, model.indTvPar),
                                        indTvPar )
    vec_of_f_t =[x -> DynNets.logLike_t(model::DynNets.GasNetModelDirBinGlobalPseudo,
                                        obsT[1:t],
                                        DynNets.restrictGasPar(model, x, model.indTvPar),
                                         model.indTvPar)
                                        for t in 1:T_train]
end

f_T(allpar0)
covA0 =  ForwardDiff.hessian(f_T, allpar0)./T_train
# run once to precompile the functions
[f(allpar0) for f in vec_of_f_t]
vec_g_t = [ForwardDiff.gradient(f, allpar0) for f in vec_of_f_t[2:end]]
global covB0 = zeros(size(covA0))
for g_t in vec_g_t
    global covB0 += g_t * g_t'
end
global covB0 = covB0./ T_train
diagCorrect = 1e-3
covMatHat = pinv(covA0) * covB0 * pinv(covA0)
covMatHat = Symmetric(covMatHat)
covMatHat_corr = (covMatHat .+ Diagonal(diagCorrect.*ones(size(covMatHat)[1])))
parDistr = MvNormal(allpar0, diag(covMatHat_corr))
samp = rand(parDistr, 6000)
indsToKeep =    (0 .< samp[end-1, :] .< 1) .& (samp[end, :] .>0) .&
                (0 .< samp[2, :] .< 1) .& (samp[3, :] .>0)
sum(indsToKeep)
hist(samp[end,:])
N_samp_par = 1000
sampled_gas_par = samp[:, indsToKeep][:, 1:N_samp_par]

gasParSampled[n] = sampled_gas_par
    gasParEst[n] = estPar # store gas parameters
    startPointEst[n] = startPoint
    gasParVec = zeros(sum(indTvPar)*3); for i=0:(sum(indTvPar)-1) gasParVec[1+3i : 3(i+1)] = estPar[indTvPar][i+1]; end
    constParVec = zeros(sum(.!indTvPar)); for i=1:sum(.!indTvPar) constParVec[i] = estPar[.!indTvPar][i][1]; end
    gasFiltPar , pseudolike = DynNets.score_driven_filter_or_dgp(model,gasParVec,indTvPar;
                              vConstPar = constParVec)#,ftot_0= startPoint )


    filtPar_T_Nsample[:,:,n] = gasFiltPar


plot(gasFiltPar)

##  Load Samles of many SD static parameters, filter the data once for each set
 # of static pars  and plot
using Utilities,AReg,StaticNets,DynNets , JLD,MLBase,StatsBase,CSV, RCall
 using PyCall; pygui(); using PyPlot
 using JLD2, GLM
 Nsample = 1
 dgpType = "sin"
 T = 50
 N = 50
 Nconf_bands = 500
 Nterms = 2
 Nsteps1 ,Nsteps2 = 0,1
 load_fold = "./data/estimatesTest/sdergmTest/conf_bands/"
 @load(load_fold*"conf_bands_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE_target_$(targetAllTv)_N_conf_$(Nconf_bands)_robust.jld" ,
     stats_T, changeStats_T,estParSS_T,sampledMat_T ,parDgpT,
     Nsample,T,N,filtPar_T_Nsample,gasParEst,convFlag,pVals_Nsample,scoreAutoc_Nsample,staticEst,gasParSampled)

gasFiltPar_conf_all = fill(zeros(N_samp_par,Nterms,T),Nsample)
 tmp = Array{Array{Real,2},2}(undef,T,Nsample); for t=1:T,s=1:Nsample tmp[t,s] =  Real.(changeStats_T[t][s]);end;#changeStats_T = tmp
 # Iterate across all SD static pars vectors
 for n=1
     @show(n)
    indTvPar = BitArray([true,true])
    model = DynNets.GasNetModelDirBinGlobalPseudo(tmp[:,n],fooGasPar,indTvPar,"")

    N_samp_par = size(gasParSampled[n])[2]
    # Works only when all parameters are time varying
    gasFiltPar_conf = zeros(N_samp_par,Nterms,T)

    for i= 1:N_samp_par
        @show(i)
        samPar = gasParSampled[n][:,i]
        gasParVec = samPar
        tmp, ~ = DynNets.score_driven_filter_or_dgp(model,gasParVec,indTvPar )
            gasFiltPar_conf[i,:,:] = tmp
    end
    gasFiltPar_conf_all[n] = gasFiltPar_conf
 end
using Statistics
# Compute confidence bands
quant_vals = [0.95, 0.05]
 conf_band = zeros(length(quant_vals),T,2)
 for t=1:T
     for k=1:Nterms
      filt_t = gasFiltPar_conf_all[n][:,k,t]
      conf_band[:,t,k] = Statistics.quantile(filt_t[.!isnan.(filt_t)],quant_vals)
    end
 end

# Plot DGP filtered Path and confidence bands
figure()
     legTex = ["DGP";"SDERGM"; "95%"]
     n=1
     indTvPar = BitArray([true,true])
     model = DynNets.GasNetModelDirBinGlobalPseudo(tmp[:,n],fooGasPar,indTvPar,"")
     estPar = gasParEst[n] # store gas parameters
     gasFiltPar  = filtPar_T_Nsample[:,:,n]
     parInd = 1
     subplot(1,2,1);
     plot(1:T,parDgpT[parInd,:],"k",linewidth=4)
     plot(1:T,gasFiltPar[parInd,:],"r")
     plt.fill_between(1:T, conf_band[1,:,parInd], y2 =conf_band[2,:,parInd],color =(0.9, 0.2 , 0.2, 0.1)  )#, color='b', alpha=.1)
     #N_bands = length(gasFiltPar_conf_all[n])
     # for i=1:N_bands
     #     plot(1:T,gasFiltPar_conf_all[n][i][parInd,:],"r",alpha = 0.02)
     # end
    parInd = 2
    subplot(1,2,2);
    plot(1:T,parDgpT[parInd,:],"k",linewidth=4)
    plot(1:T,gasFiltPar[parInd,:],"r")
     plt.fill_between(1:T, conf_band[1,:,parInd], y2 =conf_band[2,:,parInd],color =(0.9, 0.2 , 0.2, 0.1)  )#, color='b', alpha=.1)
    namePar1 = "Number of Links"
    namePar2 = "GWESP"
    subplot(1,2,1);    title(namePar1); legend(legTex)
    subplot(1,2,2);  title(namePar2 ); legend(legTex)

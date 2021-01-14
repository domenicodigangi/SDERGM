"""
Simulations to estimate coverage of confidence bands with Blasques methods and Buccheri's enhancement
"""


using ScoreDrivenExponentialRandomGraphs

import ScoreDrivenExponentialRandomGraphs:StaticNets, DynNets

import ScoreDrivenExponentialRandomGraphs.DynNets:GasNetModel,GasNetModelDirBin0Rec0

using PyPlot
pygui(true)



#region sample, estimate and filter misspecified dgp

model_mle = DynNets.fooGasNetModelDirBin0Rec0_mle
model_pmle = DynNets.fooGasNetModelDirBin0Rec0_pmle
N=30
T=200
Nsample = 3
indTvPar = trues(2)
# quick visual check on the distribution of  filters for the dgp
begin
global parMatDgp_T = DynNets.dgp_misspecified(model_mle, "sin", N, T;  minValAlpha = 0.15, maxValAlpha = 0.3)

DynNets.sample_est_mle_pmle(model_mle, parMatDgp_T, N, Nsample; plotFlag = true)
end

# sample once and  estimate SD
A_T_dgp = DynNets.sample_dgp(model_mle, parMatDgp_T,N)
stats_T_dgp = [DynNets.statsFromMat(model_mle, A_T_dgp[:,:,t]) for t in 1:T ]
change_stats_T_dgp = DynNets.change_stats(model_pmle, A_T_dgp)
obsT = stats_T_dgp
model = model_mle

estSdResPar, conv_flag, UM_mple, ftot_0_mple = DynNets.estimate(model; indTvPar=indTvPar, obsT = obsT)
vEstSdResPar = DynNets.array2VecGasPar(model, estSdResPar, indTvPar)

fVecT_filt , target_fun_val_T, sVecT_filt = DynNets.score_driven_filter(model,  vEstSdResPar, indTvPar; obsT = change_stats_T_dgp, ftot_0 = ftot_0_mple)

#endregion




#region Full Auto Diff approach to Estimate VarCovar matrix for static parameters of SD model

fVecT_filt , target_fun_val_T, sVecT_filt = DynNets.score_driven_filter( model,  vEstSdResPar, indTvPar; obsT = change_stats_T_dgp, ftot_0 = ftot_0_mple)


number_ergm_par(model::T where T <:GasNetModelDirBin0Rec0) = 2



"""
Given the flag of constant parameters, a starting value for their unconditional means (their constant value, for those constant), return a starting point for the optimization
"""
function starting_point_optim(model::T where T <:GasNetModel, indTvPar, UM; indTargPar =  falses(100))
    
    nTvPar = sum(indTvPar)
    NTargPar = sum(indTargPar)
    nErgmPar = length(indTvPar)
    
    # #set the starting points for the optimizations
    B0_Re  = 0.98; B0_Un = log(B0_Re ./ (1 .- B0_Re ))
    ARe_min =0.000001
    A0_Re  = 0.0005 ; A0_Un = log(A0_Re  .-  ARe_min)
    
    # starting values for the vector of parameters that have to be optimized
    vParOptim_0 = zeros(nErgmPar + nTvPar*2 - NTargPar)
    last = 0
    for i=1:nErgmPar
        if indTvPar[i]
            if indTargPar[i]
                vParOptim_0[last+1:last+2] = [ B0_Un; A0_Un]
                last+=2
            else
                vParOptim_0[last+1:last+3] = [UM[i]*(1 .- B0_Re) ; B0_Un; A0_Un]
                last+=3
            end
        else
            vParOptim_0[last+1] = UM[i]
            last+=1
        end
    end
    return vParOptim_0
end


"""
Given vecAllPar divide it into a vector of Score Driven parameters and one of costant parameters
"""
function divide_SD_par_from_const(model::T where T <:GasNetModel, indTvPar,  vecAllPar::Array{<:Real,1})

    nTvPar = sum(indTvPar)
    nErgmPar = length(indTvPar)

    vecSDParAll = zeros(Real,3nTvPar )
    vecConstPar = zeros(Real,nErgmPar-nTvPar)

    lastInputInd = 0
    lastIndSD = 0
    lastConstInd = 0
    #extract the vector of gas parameters, addimng w from targeting when needed
    for i=1:nErgmPar
        if indTvPar[i] 
            vecSDParAll[lastIndSD+1] = vecAllPar[lastInputInd + 1]
            vecSDParAll[lastIndSD+2] = vecAllPar[lastInputInd + 2]
            vecSDParAll[lastIndSD+3] = vecAllPar[lastInputInd + 3]
            lastInputInd +=3
            lastIndSD +=3
        else
            vecConstPar[lastConstInd+1] = vecAllPar[lastInputInd  + 1]
            lastInputInd +=1
            lastConstInd +=1
        end
    end
    return vecSDParAll, vecConstPar
end

function merge_SD_par_and_const(model::T where T <:GasNetModel, indTvPar,  vecSDPar::Array{<:Real,1}, vecConstPar)

    nTvPar = sum(indTvPar)
    nConstPar = sum(.!indTvPar)
    nConstPar == length(vConstPar) ? () : error()

    nErgmPar = length(indTvPar)
    nAllPar = 3*nTvPar + nConstPar

    vecAllPar = zeros(Real, nAllPar)

    lastIndAll = 0
    lastIndSD = 0
    lastIndConst = 0
    for i=1:nErgmPar
        if indTvPar[i] 
            
            vecAllPar[lastIndAll+1] = vecSDPar[lastIndSD + 1]
            vecAllPar[lastIndAll+2] = vecSDPar[lastIndSD + 2]
            vecAllPar[lastIndAll+3] = vecSDPar[lastIndSD + 3]

            lastIndAll +=3
            lastIndSD +=3
        else
            vecAllPar[lastIndAll+1] = vecConstPar[lastIndConst + 1]
                        
            lastInputInd +=1
            lastConstInd +=1
        end
    end
    return vecAllPar
end




"""
Restrict the  Score Driven parameters  to appropriate link functions to ensure that they remain in the region where the SD dynamics is well specified (basically 0<=B<1  A>=0)
"""
function restrict_SD_static_par(model::T where T <:GasNetModel, vecUnSDPar::Array{<:Real,1})

    nSDPar = length(vecUnSDPar)
    nTvPar = nSDPar/3

    nTvPar == 3 ? () : error()

    arrayOfVecsReSd = [ [vecUnSDPar[i], link_R_in_0_1(vecUnSDPar[i+1]), link_R_in_R_pos(vecUnSDPar[i+2]) ] for i in 1:3:nSDPar]

    vecReSDPar = reduce(vcat, arrayOfVecsReSd)

    return vecReSDPar
end


"""
Restrict the  Score Driven parameters  to appropriate link functions to ensure that they remain in the region where the SD dynamics is well specified (basically 0<=B<1  A>=0)
"""
function unrestrict_SD_static_par(model::T where T <:GasNetModel, vecReSDPar::Array{<:Real,1})

    nSDPar = length(vecReSDPar)
    nTvPar = nSDPar/3

    nTvPar == 3 ? () : error()

    arrayOfVecsUnSd = [ [vecReSDPar[i], inv_link_R_in_0_1(vecReSDPar[i+1]), inv_link_R_in_R_pos(vecReSDPar[i+2]) ] for i in 1:3:nSDPar]

    vecUnSDPar = reduce(vcat, arrayOfVecsUnSd)

    return vecUnSDPar
end


"""
"""
function target_unc_mean(UM, indTargPar)

end



using ForwardDiff

function avg_grad_and_hess_obj_SD_filter_time_seq(model, obsT, vecUnParAll, indTvPar, ftot_0)

    function obj_fun_T(xUn)

        vecSDParUn, vecConstPar = divide_SD_par_from_const(model, indTvPar, xUn)

        vecSDParRe = restrict_SD_static_par(model, vecSDParUn)

        oneInADterms  = (StaticNets.maxLargeVal + vecReSDPar[1])/StaticNets.maxLargeVal

        ~, target_fun_val_T, ~ = score_driven_filter( model,  vecSDParRe, indTvPar; obsT = obsT, vConstPar =  vecConstPar, ftot_0 = ftot_0 .* oneInADterms)
    
        return target_fun_val_T
    end

    grad = ForwardDiff.gradient(obj_fun_T, vecUnParAll)
    
    hess =  ForwardDiff.hessian(obj_fun_T, vecUnParAll)

    return grad, hess
end

objective_function_SD_filter_time_sequence(model, obsT, indTvPar, vecReSDPar, vecConstPar, ftot_0)




#endregion



#region not completely auto diff Estimate VarCovar matrix for static parameters of SD model

import ScoreDrivenExponentialRandomGraphs.DynNets:target_function_t, target_function_t_grad, target_function_t_hess

using StatsBase

function grad_and_hess_wrt_ergm_par(model, obsT, fVecT_filt )
    T = size(fVecT_filt)[2]
    grad_tot_T = zeros(2, T)
    hess_tot_T = zeros(2, 2, T)

    for t in 1:T
        obs_t = obsT[t]
        ftot_t = fVecT_filt[:,t]
        grad_tot_T[:,t] = target_function_t_grad(model, obs_t, ftot_t)
        hess_tot_T[:,:,t] = target_function_t_hess(model, obs_t, ftot_t)
    end

    return grad_tot_T, hess_tot_T
end




#endregion



#region Sample MV normal and get confidence bands


#endregion


#region Compute coverage

#endregion 


#region old explorations



#### Estimate SD static parameters
using Utilities, AReg, JLD, MLBase, StatsBase, CSV, RCall, Distributions,
      LinearAlgebra, ForwardDiff, GLM

using Revise
using StaticNets, DynNets
 using PyCall; pygui(:qt); using PyPlot

 Nconf_bands = 500
 Nterms = 2
@load("./data/congress_covoting_US/juliaEstimates.jld", estParSS_T, obsMat_T,
        changeStats_T, stats_T)
T = length(changeStats_T)
Nsample = 1
tmp = Array{Array{Real,2},2}(undef,T,Nsample)
 for t=1:T,s=1:Nsample tmp[t,s] =  Real.(changeStats_T[t]); end

estParSS_T = reshape(estParSS_T, 1, 2, 73)

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


n=1
indTvPar = BitArray([false, true])
indTargPar =indTvPar #BitArray([true,true])
obsT = tmp[:,n]
model = DynNets.GasNetModelDirBinGlobalPseudo(obsT, DynNets.fooGasPar, indTvPar,"")
estPar,convFlag[n],UM,startPoint,hess_opt = DynNets.estimate(model;UM =meanSq(estParSS_T[n,:,:],1),
                          indTvPar = indTvPar, indTargPar = indTargPar,hess_opt_flag = true)
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
    gasFiltPar , pseudolike = DynNets.score_driven_filter(model,gasParVec,indTvPar;
                              vConstPar = constParVec)#,ftot_0= startPoint )


    filtPar_T_Nsample[:,:,n] = gasFiltPar


plot(gasFiltPar)

##  Load Samles of many SD static parameters, filter the data once for each set
 # of static pars  and plot
using Utilities,AReg,StaticNets,DynNets , JLD,MLBase,StatsBase,CSV, RCall
 using PyCall; pygui(); using PyPlot
 using JLD, GLM
 Nsample = 1
 dgpType = "sin"
 T = 50
 N = 50
 Nconf_bands = 500
 Nterms = 2
 Nsteps1 ,Nsteps2 = 0,1
 load_fold = "./data/estimatesTest/sdergmTest/conf_bands/"
 @load(load_fold*"conf_bands_Nodes_$(N)_T_$(T)_Sample_$(Nsample)_Ns_" * dgpType * "_$(Nsteps1)_$(Nsteps2)_MPLE_target_$(targetAllTv)_N_conf_$(Nconf_bands)_robust.jld" ,
     stats_T, changeStats_T,estParSS_T,sampledMat_T ,parMatDgp_T,
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
        tmp, ~ = DynNets.score_driven_filter(model,gasParVec,indTvPar )
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
     plot(1:T,parMatDgp_T[parInd,:],"k",linewidth=4)
     plot(1:T,gasFiltPar[parInd,:],"r")
     plt.fill_between(1:T, conf_band[1,:,parInd], y2 =conf_band[2,:,parInd],color =(0.9, 0.2 , 0.2, 0.1)  )#, color='b', alpha=.1)
     #N_bands = length(gasFiltPar_conf_all[n])
     # for i=1:N_bands
     #     plot(1:T,gasFiltPar_conf_all[n][i][parInd,:],"r",alpha = 0.02)
     # end
    parInd = 2
    subplot(1,2,2);
    plot(1:T,parMatDgp_T[parInd,:],"k",linewidth=4)
    plot(1:T,gasFiltPar[parInd,:],"r")
     plt.fill_between(1:T, conf_band[1,:,parInd], y2 =conf_band[2,:,parInd],color =(0.9, 0.2 , 0.2, 0.1)  )#, color='b', alpha=.1)
    namePar1 = "Number of Links"
    namePar2 = "GWESP"
    subplot(1,2,1);    title(namePar1); legend(legTex)
    subplot(1,2,2);  title(namePar2 ); legend(legTex)


#endregion
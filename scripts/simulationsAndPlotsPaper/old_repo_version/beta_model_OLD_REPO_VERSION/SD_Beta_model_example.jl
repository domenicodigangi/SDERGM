

# Script that allows the user to sample from a beta model with fitnesses time varying,
# according to one of the dgps described in the paper, and then filter with the
# SD beta model

#add packages to LOAD path
push!(LOAD_PATH,"./DynNets/src/")

using  Utilities, DynNets, StaticNets, StatsBase


#%% Sample DGP

N = 10 # number of nodes
NTV = round.(Int,N/2) # Number of nodes with time varying fitnesses
T = 250 #
dynType = "SIN"
minDeg, maxDeg = degb = [3,8]
unifDeg = round(Int,(maxDeg-  minDeg)/2)

#targetStatic = true=> use a single static estimate for targeting in SD estimate
#targetStatic = false=> use a the mean of single snapshots  estimates for targeting in SD estimate
targetStatic = true
# Define time varying fitneesse
dynFitDgp, indsTVnodes =  StaticNets.dgpDynamic(StaticNets.fooErgmDirBin1,dynType,N,T;
                                  NTV = NTV, degIOUncMeans =unifDeg*ones(2N),
                                 degb = degb  )

# Sample time series of networks from the dynamical fitnesses
dynDegsSam = zeros(size(dynFitDgp))
for t=1:T
    matSamp = StaticNets.sampl(StaticNets.ErgmDirBin1(dynFitDgp[:,t]),1;  parGroupsIO = dynFitDgp[:,t])
    dynDegsSam[:,t] = [sumSq(matSamp,2);sumSq(matSamp,1)]
end
estFitSS_T =  StaticNets.estimate( StaticNets.SnapSeqNetDirBin1(dynDegsSam); identPost = true,identIter= false,targetErr = 0.001 )
for t=1:T
    dynFitDgp[:,t] =  StaticNets.identify(StaticNets.ErgmDirBin1(dynDegsSam[:,t]),dynFitDgp[:,t];idType = "firstZero" )
    estFitSS_T[:,t] =  StaticNets.identify(StaticNets.ErgmDirBin1(dynDegsSam[:,t]),estFitSS_T[:,t];idType = "firstZero" )
end

degsIO_T = dynDegsSam
indsGroups = [Int.(1:2N),Int.(1:2N)] #[Int.(1:2N),ones(2N)] #
modGasDirBin1_3Npars = DynNets.SdErgmDirBin1(degsIO_T, [zeros(2N), zeros(2N) ,zeros(2N)],
                                            indsGroups, "FISHER-DIAG")
rmseSS =sqrt.(meanSq((dynFitDgp - estFitSS_T).^2,2))
if targetStatic
    estFitStatic = DynNets.estimateSnapSeq(modGasDirBin1_3Npars, degsIO_T = mean(degsIO_T,dims=2) )
    estFitSS_T = repeat(estFitStatic, 1, 5)
end

#%% Estimate SD beta model ---Can take a while (~30 min for N=10)
@time parGas_3Npars, ~  =  DynNets.estimateTarg(modGasDirBin1_3Npars;SSest =estFitSS_T)
gasFiltFit,~ = DynNets.score_driven_filter_or_dgp(modGasDirBin1_3Npars,[parGas_3Npars[1];parGas_3Npars[2];parGas_3Npars[3]])
rmseGas =sqrt.(meanSq((dynFitDgp - gasFiltFit).^2,2))

# plot the resulting filtered path
indsPlot = [ 5]
for Nplot =indsPlot
    plot(1:T,dynFitDgp[Nplot,:])
    plot(1:T,filFitGas[Nplot,:,:],".r")
    plot(1:T,estFitSS[Nplot,:,:],".b")
end

#%%Add Confidence bands accounting for parameter uncertainty
using ForwardDiff

# TO BE ADJUSTED TO THE CASE WITHOUT PARAMETERS RESTRICTIONS1!!
# NOW IT HANDLES ONLY THE RESTRICTED CASE, AND MIGHT CONTAIN BUGS

Ttrain = T# round(Int, T/2) #70 106 #
#Estimate gas model on  windows of length Ttrain
modGasDirBin1_N_p_2 = SdErgmDirBin1(degsIO_T[:,1:Ttrain],"FISHER-DIAG", "ONE_PAR_ALL")
gasParEstOnTrain,~ = estimateTarg(modGasDirBin1_N_p_2;SSest = estFitSS_T[:,1:Ttrain] )

#f(x) = DynNets.score_driven_filter_or_dgp(modGasDirBin1_N_p_2, x;groupsInds = modGasDirBin1_N_p_2.groupsInds)[2]
GBA = 1
allpar0 = array_2_vec_all_par(modGasDirBin1_N_p_2, gasParEstOnTrain)
BA_re_hat = allpar0[end-1:end]
function UnrestrAB(vecRePar::Array{<:Real,1}) #restrict a vector of only A and B
      diag_B_Re = vecRePar[1:GBA]
      diag_A_Re = vecRePar[GBA+1:end]
      diag_B_Un = log.(diag_B_Re ./ (1 .- diag_B_Re ))
      diag_A_Un = log.(diag_A_Re)
      vecUnPar =  [diag_B_Un; diag_A_Un]
      return vecUnPar
end
function loglike(Model::DynNets.SdErgmDirBin1, degs_t::Array{<:Real,1},
                        f_t::Array{<:Real,1})

    thetas_mat_t_exp, exp_mat_t = StaticNets.expMatrix2(StaticNets.fooErgmBin1,f_t)
    exp_deg_t = sum(exp_mat_t,dims = 2)

    loglike_t = sum(f_t.*degs_t) -  sum(UpperTriangular(log.(1 .+ thetas_mat_t_exp))) #  sum(log.(1 + thetas_mat_t_exp))
    return loglike_t
 end
function logLike_t(Model, obsT, vReGasPar)
 groupsInds = Model.groupsInds
      N2,T = size(obsT);N = round(Int,N2/2)
      NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model, groupsInds)
      # Organize parameters of the GAS update equation
      WGroupsIO = vReGasPar[1:NGW]
      #    StaticNets.expMatrix2(StaticNets.fooErgmDirBin1,WGroupsIO )
      W_allIO = WGroupsIO[groupsInds[1]]
      BgasGroups  = vReGasPar[NGW+1:NGW+GBA]
      AgasGroups  = vReGasPar[NGW+GBA+1:NGW+2GBA]
      #distribute nodes among groups with equal parameters
      AgasIO = AgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
      BgasIO = BgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
      WgasIO   = W_allIO[indTvNodesIO]
      UMallNodesIO = W_allIO
      UMallNodesIO  =  WgasIO ./ (1 .- BgasIO)
      ftotIO_t =  identify(Model,UMallNodesIO)
      I_tm1 = UniformScaling(N)
      loglike_t = zero(Real)
      for t=1:T-1
            degsIO_t = obsT[:,t] # vector of in and out degrees
            ftotIO_tp1,loglike_t = predict_score_driven_par(Model,N,degsIO_t,ftotIO_t,I_tm1,
                                                  indTvNodesIO,WgasIO,BgasIO,AgasIO)
            ftotIO_t = ftotIO_tp1
      end
      return  loglike_t::T where T <:Real
end
function logLike_T(Model, obsT, vReGasPar)
      groupsInds = Model.groupsInds
      N2,T = size(obsT);N = round(Int,N2/2)
      NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model, groupsInds)
      # Organize parameters of the GAS update equation
      WGroupsIO = vReGasPar[1:NGW]
      #    StaticNets.expMatrix2(StaticNets.fooErgmDirBin1,WGroupsIO )
      W_allIO = WGroupsIO[groupsInds[1]]
      BgasGroups  = vReGasPar[NGW+1:NGW+GBA]
      AgasGroups  = vReGasPar[NGW+GBA+1:NGW+2GBA]
      #distribute nodes among groups with equal parameters
      AgasIO = AgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
      BgasIO = BgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
      WgasIO   = W_allIO[indTvNodesIO]
      # start values equal the unconditional mean, but  constant ones remain equal to the unconditional mean, hence initialize as:
      UMallNodesIO = W_allIO
      UMallNodesIO  =  WgasIO ./ (1 .- BgasIO)
      ftotIO_t =  identify(Model,UMallNodesIO)
      I_tm1 = UniformScaling(N)
      loglike_T = zero(Real)
      for t=1:T-1
            degsIO_t = obsT[:,t] # vector of in and out degrees
            ftotIO_tp1,loglike_t = predict_score_driven_par(Model,N,degsIO_t,ftotIO_t,I_tm1,
                                                  indTvNodesIO,WgasIO,BgasIO,AgasIO)
            ftotIO_t = ftotIO_tp1
            loglike_T += loglike_t
      end
      return  loglike_T::T where T <:Real
end


obsT = modGasDirBin1_N_p_2.obsT
T_train = size(obsT)[2]
f_T(x) = logLike_T(modGasDirBin1_N_p_2::SdErgmDirBin1, obsT, x)
# @time ForwardDiff.hessian(f_t, allpar0)
#@time FiniteDiff.finite_difference_hessian(f_T, allpar0)

covA0 =  ForwardDiff.hessian(f_T, allpar0)./T_train
# estimate outer product of scores
vec_of_f_t =[x -> logLike_t(modGasDirBin1_N_p_2::SdErgmDirBin1, obsT[:, 1:t], x) for t in 1:Ttrain]
# run once to precompile the functions
[f(allpar0) for f in vec_of_f_t]
vec_g_t = [ForwardDiff.gradient(f, allpar0) for f in vec_of_f_t[2:end]]
global covB0 = zeros(size(covA0))
 for g_t in vec_g_t
      global covB0 += g_t * g_t'
 end
 global covB0 = covB0./ T_train
#%%
diagCorrect = 1e-3
covMatHat = pinv(covA0) * covB0 * pinv(covA0)
      covMatHat = Symmetric(covMatHat)

covMatHat_corr = (covMatHat .+ Diagonal(diagCorrect.*ones(size(covMatHat)[1])))

parDistr = MvNormal(allpar0, diag(covMatHat_corr))

samp = rand(parDistr, 2000)
indsToKeep = (0 .< samp[end-1, :] .< 1) .& (samp[end, :] .>0)
sum(indsToKeep)
#hist(samp[end,:])
N_samp_par = 1000
allpar_sample = samp[:, indsToKeep][:, 1:N_samp_par]

#%%
ftotIO_0 = meanSq(estimateSnapSeq(modGasDirBin1_N_p_2; degsIO_T = obsT[:, 1:5]), 2)
ftot_T_est, _ = score_driven_filter_or_dgp(modGasDirBin1_N_p_2, allpar0,ftotIO_0 = ftotIO_0)
ftot_T, _ = score_driven_filter_or_dgp(modGasDirBin1_N_p_2, allpar_sample[:,3],  ftotIO_0 = ftotIO_0)
plot(ftot_T_est')
#plot(ftot_T')

Nterms = 2*N

global gasFiltPar_conf = zeros(N_samp_par, Nterms, T_train)
      for i= 1:N_samp_par
       # @show(i)
       global gasFiltPar_conf[i,:,:], _ = score_driven_filter_or_dgp(modGasDirBin1_N_p_2, allpar_sample[:,i],  ftotIO_0 = ftotIO_0)
      end

indsToRem = [any(isnan.(gasFiltPar_conf[i, : , :])) for i in 1:N_samp_par ]
sum(indsToRem)
gasFiltPar_conf = gasFiltPar_conf[ .!indsToRem, : , : ]

using Statistics
# Compute confidence bands
quant_vals = [0.95, 0.05]
 conf_band = zeros(length(quant_vals),T_train, Nterms)
 for t=1:T_train
     for k=1:Nterms
      filt_t = gasFiltPar_conf[:,k,t]
      conf_band[:,t,k] = Statistics.quantile(filt_t[.!isnan.(filt_t)],quant_vals)
    end
 end

figure()
      legTex = ["DGP";"SDERGM"; "95%"]
      node = 11
      plot(1:T_train,ftot_T_est[node, :],"r")
      plt.fill_between(1:T_train, conf_band[1, :, node], y2 =conf_band[2,:,node],color =(0.9, 0.2 , 0.2, 0.1)  )#, color='b', alpha=.1)

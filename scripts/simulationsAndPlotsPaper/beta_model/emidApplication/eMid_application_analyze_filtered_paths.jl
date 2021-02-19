"""
Analyze filtered fitnesses for emid data. Code created as part of the revision of
   SDERGM paper
"""

## Load dataj
using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase#,DynNets
using PyCall; pygui(:qt); using PyPlot


#estimate and save for half the datase (after LTRO) or whole data?
halfPeriod = false
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod ? periodEndStr =  "2012_03_12_to_2015_02_27.jld" : periodEndStr =
                                                "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,
      inactiveBanks,
      YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)

N_test = 50

# single snapshot estimates
matY_T = YeMidWeekly_T[1:N_test,1:N_test,3:end]
 matA_T = matY_T.>0
    N = length(matA_T[1,:,1])
    T = size(matA_T)[3]
    Ttrain = 100#round(Int, T/2) #70 106 #
    threshVar = 0.00#5
    degsIO_T = [sumSq(matA_T,2);sumSq(matA_T,1)]
 allFitSS =  StaticNets.estimate( StaticNets.SnapSeqNetDirBin1(degsIO_T);
            identPost = false,identIter= true )


#plot(allFitSS[:,1])

#Estimate gas model on  windows of length Ttrain
model = GasNetModelDirBin1(degsIO_T[:,1:Ttrain],"FISHER-DIAG", "ONE_PAR_ALL")
gasParEstOnTrain,~ = estimateTarg(model;SSest = allFitSS[:,1:Ttrain] )

#f(x) = DynNets.score_driven_filter_or_dgp(model, x;groupsInds = model.groupsInds)[2]
GBA = 1
allpar0 = array2VecGasPar(model, gasParEstOnTrain)
BA_re_hat = allpar0[end-1:end]
function UnrestrAB(vecRePar::Array{<:Real,1}) #restrict a vector of only A and B
      diag_B_Re = vecRePar[1:GBA]
      diag_A_Re = vecRePar[GBA+1:end]
      diag_B_Un = log.(diag_B_Re ./ (1 .- diag_B_Re ))
      diag_A_Un = log.(diag_A_Re)
      vecUnPar =  [diag_B_Un; diag_A_Un]
      return vecUnPar
end

BA_un_hat =  UnrestrAB(BA_re_hat)
w_hat = allpar0[1:2*N]
#%%
function loglike(Model::DynNets.GasNetModelDirBin1, degs_t::Array{<:Real,1},
                        f_t::Array{<:Real,1})

    thetas_mat_t_exp, exp_mat_t = StaticNets.expMatrix2(StaticNets.fooNetModelBin1,f_t)
    exp_deg_t = sum(exp_mat_t,dims = 2)

    loglike_t = sum(f_t.*degs_t) -  sum(UpperTriangular(log.(1 .+ thetas_mat_t_exp))) #  sum(log.(1 + thetas_mat_t_exp))
    return loglike_t
 end



using ForwardDiff
obsT = model.obsT
T_train = size(obsT)[2]
f_T(x) = logLike_T(model::GasNetModelDirBin1, obsT, x)
# @time ForwardDiff.hessian(f_t, allpar0)
#@time FiniteDiff.finite_difference_hessian(f_T, allpar0)

covA0 =  ForwardDiff.hessian(f_T, allpar0)./T_train
# estimate outer product of scores
vec_of_f_t =[x -> logLike_t(model::GasNetModelDirBin1, obsT[:, 1:t], x) for t in 1:Ttrain]
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
ftotIO_0 = meanSq(estimateSnapSeq(model; degsIO_T = obsT[:, 1:5]), 2)
ftot_T_est, _ = score_driven_filter_or_dgp(model, allpar0,ftotIO_0 = ftotIO_0)
ftot_T, _ = score_driven_filter_or_dgp(model, allpar_sample[:,3],  ftotIO_0 = ftotIO_0)
plot(ftot_T_est')
#plot(ftot_T')

Nterms = 2*N

global gasFiltPar_conf = zeros(N_samp_par, Nterms, T_train)
      for i= 1:N_samp_par
       # @show(i)
       global gasFiltPar_conf[i,:,:], _ = score_driven_filter_or_dgp(model, allpar_sample[:,i],  ftotIO_0 = ftotIO_0)
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

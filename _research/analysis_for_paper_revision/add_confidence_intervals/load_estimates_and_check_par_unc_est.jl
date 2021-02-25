"""
Load estimates from known DGPs and check that the parameter uncertainty is estimated sufficiently well
"""

#region import and models
using Pkg
Pkg.activate(".") 
Pkg.instantiate() 
using DrWatson
using DataFrames
using ScoreDrivenERGM
using PyPlot
pygui(true)
using Statistics
using Logging
using SharedArrays
#endregion

import ScoreDrivenERGM:DynNets


@time df = collect_results( datadir("sims", "samDgpFiltSD_est_conf"))
a = 1

names(df)


# select results for one set of parameters
begin
T = 300
N = 300
modelTag = string(DynNets.GasNetModelDirBin0Rec0_pmle())
dgpType = "SD"
nSample = 100
res = filter([:modelTag, :T, :N, :dgpSettings, :nSample] => (m,t,n, d, s) -> all((m==modelTag, t==T, n==N, d.type == dgpType, s == nSample)), df)[1,:]

# visual inspection of filter, dgp and conf bands
if true
    for nPlot = 1:20
        if !any(res.errInds[:,nPlot])
            DynNets.plot_filtered_and_conf_bands(res.model, N, res.allfVecT_filt[:,:,nPlot], res.allConfBandsFiltPar[:,:,:,:,nPlot] ; parDgpT=res.allParDgpT[:,:,nPlot], confBands2 = res.allConfBandsPar[:,:,:,:,nPlot] )
        end
    end
end
end

#region Distribution of static SD parameters
begin


indPar = 6
fig, ax = subplots(2,2)
ax[1,1].hist(res.allvEstSdResPar[2, :], 30)
ylim = ax[1,1].get_ylim()
ax[1,1].vlines(res.dgpSettings.opt.B, ylim[1], ylim[2], color = "r" )
ax[2,1].hist(res.allvEstSdResPar[5, :], 30)
ylim = ax[2,1].get_ylim()
ax[2,1].vlines(res.dgpSettings.opt.B, ylim[1], ylim[2], color = "r" )
ax[1,2].hist(res.allvEstSdResPar[3, :], 30)
ylim = ax[1,2].get_ylim()
try 
    ax[1,2].vlines(res.dgpSettings.opt.A, ylim[1], ylim[2], color = "r" )
catch
end
ax[2,2].hist(res.allvEstSdResPar[6, :], 30)
ylim = ax[2,2].get_ylim()
ax[2,2].vlines(res.dgpSettings.opt.A, ylim[1], ylim[2], color = "r" )
suptitle("Estimators' distributions for Static SD parameters  nSample = $nSample , T = $T, N = $N")
tight_layout()

end
#endregion

#region distribution of white estimators
begin
# unrestrict static SD estimates
allvEstSdUnPar =  mapslices(x->DynNets.unrestrict_all_par(res.model, res.model.indTvPar, x),  res.allvEstSdResPar, dims=1)
covParBoot = cov(permutedims(allvEstSdUnPar))
allCovWhite = SharedArray(zeros(6,6,nSample))
allCovErr = SharedArray{Bool}(nSample)
Threads.@threads for n=1:nSample
   allCovWhite[:,:,n], allCovErr[n] = DynNets.white_estimate_cov_mat_static_sd_par(res.model, N, res.allObsT[n], res.model.indTvPar, res.allftot_0[:, n], res.allvEstSdResPar[:,n])
end 

parNames = ["w_theta", "B_theta", "A_theta", "w_eta", "B_eta", "A_eta"]

fig, ax = subplots(3,2)
for i = 1:6
    j=i
    ax[i].hist(allCovWhite[i,j,:], range=quantile(allCovWhite[i,j,:], [0.01, 0.99]), 20)
    ylim = ax[i].get_ylim()
    ax[i].vlines(covParBoot[i,j], ylim[1], ylim[2], color = "r" )
    ax[i].set_title(parNames[i])
end
suptitle("White Estimators of Diagonal Elements of Covariance Matrix \n Par Bootstrap Estimate in Red. nSample = $nSample , T = $T, N = $N")
tight_layout()

end
#endregion

# Domande
# 1. C'é un bias di campione finito. Lo stimatore, MLE o White, della covarianza degli stimatori dei parametri é ok? 
# 2. Visto che posso costruire la distribuzione dei parametri con parametric bootstrap, perché non usare quella?
# 3. Che relazione c'é tra la distribuzione da parametric bootstrap e la mia idea di usare non parametric bootstrap campionando le osservazioni (che portano in dote la storia fino a quel punto, ma considerando la likelihood della sola osservazione)? 


#endregion


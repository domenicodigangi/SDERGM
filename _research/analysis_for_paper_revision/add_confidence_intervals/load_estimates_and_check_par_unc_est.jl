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
import ScoreDrivenERGM:DynNets, Utilities
using PyPlot
pygui(true)
using Statistics
using Logging
using SharedArrays
#endregion

import ScoreDrivenERGM:DynNets


@elapsed df = collect_results( datadir("sims", "samDgpFiltSD_est"))
@elapsed df = collect_results( datadir("sims", "samDgpFiltSD_conf"))
df["modelTag"] = string.(df["model"]) 

# select results for one set of parameters
begin
T = 3000
N = 100
nSample = 50
indQuant = 1
modelTag = string(DynNets.GasNetModelDirBin0Rec0_pmle(scoreScalingType = "FISH_D"))
modelTag = string(DynNets.GasNetModelDirBin0Rec0_mle(scoreScalingType = "FISH_D"))
parUncMethod = "WHITE-MLE" 

dgpSetting = DynNets.list_example_dgp_settings(DynNets.GasNetModelDirBin0Rec0_mle()).dgpSetARlowlow


         

# # visual inspection of filter, dgp and conf bands
res = filter([:modelTag, :T, :N, :dgp, :S, :m] => (mTag,t,n, d, s, m) -> all((mTag==modelTag, t==T, n==N, d == dgpSetting, s == nSample, m == parUncMethod)), df)

nrow(res) != 1 ? error("$N, $T,  $(size(res))") : res = res[1,:]
    for nPlot = 1:1
        if !any(res.errInds[:,nPlot])
            DynNets.plot_filtered_and_conf_bands(res.model, N, res.allfVecT_filt[:,1:end,nPlot], res.allConfBandsFiltPar[:,1:end,:,:,nPlot] ; parDgpTIn=res.allParDgpT[:,1:T,nPlot], confBands2In = res.allConfBandsPar[:,1:end,:,:,nPlot] , offset = 1, indBand = indQuant)
        end
    end

using Statistics
using SharedArrays
# unrestrict static SD estimates

if false

allCovWhite = res.allmvSDUnParEstCovWhite
covParBoot = cov(vEstSdUnParBootDist')
covParBoot = cov(Utilities.drop_bad_un_estimates(vEstSdUnParBootDist)')

parNames = [ "B_theta", "A_theta",  "B_eta", "A_eta"]
fig, ax = subplots(2,2)
for i = 3:6
    j=i
    ax[i-2].hist(log10.(allCovWhite[i,j,:]), range=quantile(log10.(allCovWhite[i,j,:]), [0.05, 0.95]), 20)
    ylim = ax[i-2].get_ylim()
    ax[i-2].vlines(log10(covParBoot[i,j]), ylim[1], ylim[2], color = "r" )
    ax[i-2].set_title("Var( $(parNames[i-2]))")
end
suptitle("Empirical Distribution of White Estimators \n $(res.modelTag) \n  nSample = $nSample , T = $T, N = $N, A=$(res.dgp.opt.A[1])")
tight_layout()

end






res = filter([:modelTag, :T, :N, :dgpSettings, :nSample] => (mTag,t,n, d, s) -> all((mTag==modelTag, t==T, n==N, d == dgpSetting, s == nSample)), df)

nrow(res) != 1 ? error("$N, $T,  $(size(res))") : res = res[1,:]


vEstSdResParBootDist =  res.allvEstSdResPar
vEstSdUnParBootDist = mapslices(x -> DynNets.unrestrict_all_par(res.model, res.model.indTvPar, x), vEstSdResParBootDist, dims=1)

parNames = ["w_theta", "B_theta", "A_theta", "w_eta", "B_eta", "A_eta"]
if false
sample1 = vEstSdResParBootDist
sample1 = vEstSdUnParBootDist
sample1 = Utilities.drop_bad_un_estimates(vEstSdUnParBootDist)
fig, ax = subplots(3,2)
for i = 1:6
    ax[i].hist(sample1[i,:], range=quantile(sample1[i,:], [0.01, 0.99]), 20, alpha = 0.4)
    ax[i].set_title(parNames[i])
    ylim = ax[2,1].get_ylim()
    # ax[i].vlines(res.vEstSdResPar[i], ylim[1], ylim[2], color = "r" )

end
end
if false
fig, ax = subplots(2,2)
ax[1,1].hist(res.allvEstSdResPar[2, :], 30)
ylim = ax[1,1].get_ylim()
ax[1,1].vlines(res.dgpSettings.opt.B, ylim[1], ylim[2], color = "r" )
ax[2,1].hist(res.allvEstSdResPar[5, :], 30)
ylim = ax[2,1].get_ylim()
ax[2,1].vlines(res.dgpSettings.opt.B, ylim[1], ylim[2], color = "r" )
ax[1,2].hist(res.allvEstSdResPar[3, :], 30)
ylim = ax[1,2].get_ylim()
ax[1,2].vlines(res.dgpSettings.opt.A, ylim[1], ylim[2], color = "r" )
ax[2,2].hist(res.allvEstSdResPar[6, :], 30)
ax[2,2].vlines(res.dgpSettings.opt.A, ylim[1], ylim[2], color = "r" )
ylim = ax[2,2].get_ylim()
# ax[2,2].vlines(res.dgpSettings.opt.A, ylim[1], ylim[2], color = "r" )
suptitle("Estimators' distributions for Static SD parameters  nSample = $nSample , T = $T, N = $N")
tight_layout()

end
end



# Domande
# 1. C'é un bias di campione finito. Lo stimatore, MLE o White, della covarianza degli stimatori dei parametri é ok? 
# 2. Visto che posso costruire la distribuzione dei parametri con parametric bootstrap, perché non usare quella?
# 3. Che relazione c'é tra la distribuzione da parametric bootstrap e la mia idea di usare non parametric bootstrap campionando le osservazioni (che portano in dote la storia fino a quel punto, ma considerando la likelihood della sola osservazione)? 


#endregion



"""
Test script for SD dirBin0Rec0 model_mle: one parameter for total number of links and one for reciprocity
"""

# test sampling

include("...\\..\\..\\..\\..\\add_load_paths.jl")
using StaticNets
using DynNets
using PyPlot

sn = StaticNets
## define parameters
N=30
θ_0 = 1
η_0 = -1.5
T =100

## Sample SD gdp
model_mle = fooGasNetModelDirBin0Rec0_mle
indTvPar = trues(2)
B0=0.95
A0 = 0.01
aResGasPar_0 = [[θ_0*(1-B0), B0, A0], [η_0*(1-B0), B0, A0] ]
vResGasPar_0 = DynNets.array2VecGasPar(model_mle, aResGasPar_0, indTvPar)
vResGasPar = vResGasPar_0
ftot_0= zeros(Real,2)

fig, ax = subplots(2,1)
for i=1:100
fVecT_dgp , A_T_dgp, sVecT_dgp = score_driven_filter( model_mle,  vResGasPar_0, indTvPar; dgpNT = (N,T))
ax[1].plot(fVecT_dgp[1,:])
ax[1].set_ylim([θ_0 - 5, θ_0 + 5])
ax[2].plot(fVecT_dgp[2,:])
ax[2].set_ylim([η_0 - 5, η_0 + 5])
end


## filter SD
fVecT_dgp , A_T_dgp, sVecT_dgp = score_driven_filter( model_mle,  vResGasPar_0, indTvPar; dgpNT = (N,T))
stats_T_dgp = [statsFromMat(model_mle, A_T_dgp[:,:,t]) for t in 1:T ]
change_stats_T_dgp = change_stats(model_mle, A_T_dgp)

fVecT_filt , target_fun_val_T, sVecT_filt = score_driven_filter( model_mle,  vResGasPar_0, indTvPar; obsT = stats_T_dgp)
fVecT_filt_p , target_fun_val_T_p, sVecT_filt_p = score_driven_filter( model_mple,  vResGasPar_0, indTvPar; obsT = change_stats_T_dgp)

fig, ax = subplots(2,1)
ax[1].plot(fVecT_filt[1,:], "b")
ax[1].plot(fVecT_filt_p[1,:], "r")
ax[1].plot(fVecT_dgp[1,:], "k")
ax[1].set_ylim([θ_0 - 5, θ_0 + 5])
ax[2].plot(fVecT_filt[2,:], "b")
ax[2].plot(fVecT_filt_p[2,:], "r")
ax[2].plot(fVecT_dgp[2,:], "k")
ax[2].set_ylim([η_0 - 5, η_0 + 5])


## estimate SD
model_mple = fooGasNetModelDirBin0Rec0_pmle

indTargPar = indTvPar
fVecT_dgp , A_T_dgp, sVecT_dgp = score_driven_filter( model_mle,  vResGasPar_0, indTvPar; dgpNT = (N,T))
stats_T_dgp = [statsFromMat(model_mle, A_T_dgp[:,:,t]) for t in 1:T ]
change_stats_T_dgp = change_stats(model_mle, A_T_dgp)

estimate(model_mle; indTvPar=indTvPar, indTargPar=indTargPar, obsT = stats_T_dgp)

estPar_mple, conv_flag,UM_mple , ftot_0_mple = estimate(model_mple; indTvPar=indTvPar, indTargPar=indTargPar, obsT = change_stats_T_dgp)
vResEstPar_mple = DynNets.array2VecGasPar(model_mple, estPar_mple, indTvPar)
fVecT_filt_p , target_fun_val_T_p, sVecT_filt_p = score_driven_filter( model_mple,  vResEstPar_mple, indTvPar; obsT = change_stats_T_dgp)








#

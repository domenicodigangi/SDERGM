
"""
Test script for SD dirBin0Rec0 model: one parameter for total number of links and one for reciprocity
"""

# test sampling
using Revise
using StaticNets, DynNets

## define parameters
N=30
θ_0 = -3.1
η_0 = 2.7


## Sample SD gdp
model = DynNets.fooGasNetModelDirBin0Rec0
indTvPar = trues(2)
B0=0.95
A0 = 0.01
aResGasPar_0 = [[θ_0*(1-B0), B0, A0], [η_0*(1-B0), B0, A0] ]
vResGasPar_0 = DynNets.array2VecGasPar(model, aResGasPar_0, indTvPar)

StaticNets.logLikelihood(StaticNets.fooNetModelDirBin0Rec0, 12, 5, [0, 0])
fVecT_dgp , A_T_dgp = DynNets.gasFilter( model, vResGasPar_0, indTvPar; dgpNT = (10,100))


## filter SD

## estimate SD


















#

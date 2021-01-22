module StaticNets

using Distributions, StatsBase,Optim, LineSearches, StatsFuns,Roots,MLBase, Statistics, LinearAlgebra, Random

using ..AReg
using ..Utilities
using ..Scalings
using ..ErgmRcall



## STATIC NETWORK MODEL
abstract type NetModel end
abstract type NetModelW <: NetModel end
abstract type NetModelBinW <: NetModel end
abstract type NetModelWcount <: NetModelW end
abstract type NetModelBin <: NetModel end
abstract type NetModelDirBin <: NetModelBin end

#constants
targetErrValStaticNets = 1e-2
export targetErrValStaticNets

targetErrValStaticNetsW = 1e-5
export targetErrValStaticNetsW

bigConstVal = 10^6
export bigConstVal

maxLargeVal =  1e40# 1e-10 *sqrt(prevfloat(Inf))
export maxLargeVal

minSmallVal = 1e2*eps()
export minSmallVal

include("./StaticNets_models/StaticNets_Bin1.jl")
include("./StaticNets_models/StaticNets_DirBin1.jl")
include("./StaticNets_models/StaticNets_DirBin0Rec0.jl")
include("./StaticNets_models/StaticNets_DirBinErgmPseudo.jl")


end

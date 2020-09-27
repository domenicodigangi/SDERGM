__precompile__(false)
module StaNets

using Distributions, StatsBase,Optim, LineSearches, StatsFuns,Roots,MLBase, Statistics, LinearAlgebra, Random
# using PyCall;# pygui(:qt);
# using PyPlot

export maxLargeVal

## STATIC NETWORK MODEL
abstract type NetModel end
abstract type NetModelW <: NetModel end
abstract type NetModelBinW <: NetModel end
abstract type NetModelWcount <: NetModelW end
abstract type NetModelBin <: NetModel end

#constants
targetErrValStaNets = 1e-2
targetErrValStaNetsW = 1e-5
bigConstVal = 10^6
maxLargeVal =  1e40# 1e-10 *sqrt(prevfloat(Inf))
minSmallVal = 1e2*eps()

include("./AReg.jl")

include("./HelperFunDom.jl")

include("./StaNets_Bin1.jl")
include("./StaNets_DirBin1.jl")
include("./StaNets_DirBin1Rec0.jl")


end

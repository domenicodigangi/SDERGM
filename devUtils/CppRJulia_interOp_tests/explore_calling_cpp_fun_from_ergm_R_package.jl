
"""
Can we avoid stepping trough R to leverage functionality form ergm R package? Most of it is written in C++.

Not doable because Cxx.jl does not work with juia 1.5


"""

using Cxx
import ScoreDrivenExponentialRandomGraphs:ErgmRcall

using ScoreDrivenExponentialRandomGraphs.Scalings
using ScoreDrivenExponentialRandomGraphs.Utilities
using ScoreDrivenExponentialRandomGraphs:StaticNets

using PyPlot
using RCall
using Statistics

ErgmRcall.clean_start_RCall()
ergmTermsString = "edges +  mutual"
R"""options(warn=-1) """

model = StaticNets.fooNetModelDirBin0Rec0


##-------------------- Test and COmpare MLE and MPLE estimates
#Compare for a single value of the parameters
N=100
Ldgp = n_pox_dir_links(N)/4
Rdgp = n_pox_pairs(N)/4-100

θ_0, η_0 = Tuple(StaticNets.estimate(model, Ldgp, Rdgp, N))

nSample = 100
diadProb = StaticNets.diadProbFromPars(model, [θ_0, η_0])
A_vec = [StaticNets.samplSingMatCan(model, diadProb, N) for i=1:nSample]
statsVec = reduce(hcat,[StaticNets.statsFromMat(model, A) for A in A_vec])


jnum = 10

cxx"""
        void printme(int x) {
            std::cout << x << std::endl;
        }
    """
@cxx printme(jnum)



function c_list(model, )
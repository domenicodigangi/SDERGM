
"""
Test script for dirBin0Rec0 model: one parameter for total number of links and one for reciprocity. 
Investigate the existence of regions of the θ η plane where the latters cannot be estimated for finite N  
"""

using ScoreDrivenExponentialRandomGraphs
using ScoreDrivenExponentialRandomGraphs.Utilities

import ScoreDrivenExponentialRandomGraphs.StaticNets: fooNetModelDirBin0Rec0, ergm_par_from_mean_vals,diadProbFromPars , samplSingMatCan, statsFromMat, get_mple, estimate,NetModelDirBin0Rec0, exp_val_stats, sample_ergm

using PyPlot
pygui(true)
using RCall
using Statistics


ErgmRcall.clean_start_RCall()
ergmTermsString = "edges +  mutual"
model = fooNetModelDirBin0Rec0
    

function sample_dgp_static_estimate_mle(model::NetModelDirBin0Rec0, N, parDgpSeq, nSample)

    nPars = size(parDgpSeq)[2]
    parMle = zeros(nPars, 2, nSample)
    obs = zeros(nPars, 3, nSample)

    for i =1:nPars
    
        θ_0, η_0 = parDgpSeq[:,i]

        A_vec = sample_ergm(model, N, θ_0, η_0, nSample )        

        obs[i,:,:] = reduce(hcat,[statsFromMat(model, A) for A in A_vec])
        parMle[i,:,:] = reduce(hcat,[estimate(model, A) for A in A_vec])


    end

    return parMle, obs
end




begin 
figure()
Nvals = [20, 30, 50, 70, 90]
for N in Nvals 
    θ_0 = -5
    parDgpSeq = reduce(hcat, [ [θ_0, η ] for η = 2:0.05:7])
    parMle, obs = sample_dgp_static_estimate_mle(model, N, parDgpSeq, 150 )
    fracZeroRec = mean(obs[:,2,:].==0;dims=2)
    xplot = parDgpSeq[2,:]
    plot(xplot, fracZeroRec)

end 
    title("Fraction of samples with zero reciprocal pairs θ_dgp = $θ_0 ")
    xlabel( "η_dgp")
    legend(Nvals)
    grid()
end
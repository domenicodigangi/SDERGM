
"""
Test script for dirBin0Rec0 model: one parameter for total number of links and one for reciprocity. 
Investigate the existence of regions of the θ η plane where the latters cannot be estimated for finite N  
"""

using ScoreDrivenERGM
using ScoreDrivenERGM.Utilities

import ScoreDrivenERGM.StaticNets: ErgmDirBin0Rec0(), ergm_par_from_mean_vals,diadProbFromPars , samplSingMatCan, stats_from_mat, get_one_mple, estimate,ErgmDirBin0Rec0, exp_val_stats, sample_ergm, estimate

using PyPlot
pygui(true)
using RCall
using Statistics


ErgmRcall.clean_start_RCall()
ergmTermsString = "edges +  mutual"
model = ErgmDirBin0Rec0()
    

function sample_mats_sequence_static_estimate_mle(model::ErgmDirBin0Rec0, N, parDgpSeq, nSample)

    T = size(parDgpSeq)[2]
    parMle = zeros(T, 2, nSample)
    obs = zeros(T, 3, nSample)

    for i =1:T
    
        A_vec = sample_ergm(model, N, parDgpSeq[:,i], nSample )        

        obs[i,:,:] = reduce(hcat,[StaticNets.stats_from_mat(model, A) for A in A_vec])
        parMle[i,:,:] = reduce(hcat,[StaticNets.estimate(model, A) for A in A_vec])
    end

    return parMle, obs
end



N=30
begin 
gix, ax = subplots(2,1)
Nvals = [30, 50, 70, 90]
θ_0 = -5
for N in Nvals 
    parDgpSeq = reduce(hcat, [ [θ_0, η ] for η = 2:0.05:7])
    parMle, obs = sample_mats_sequence_static_estimate_mle(model, N, parDgpSeq, 150 )
    fracZeroRec = mean(obs[:,2,:].==0;dims=2)
    fracUb = mean(obs[:,2,:].==obs[:,1,:]/2;dims=2)
    xplot = parDgpSeq[2,:]
    ax[1].plot(xplot, fracZeroRec)
    ax[2].plot(xplot, fracUb)

end 
    ax[1].set_title("Fraction of samples with zero reciprocal pairs θ_dgp = $θ_0 ")
    ax[1].set_xlabel( "η_dgp")
    ax[1].legend(Nvals)
    ax[1].grid()
    ax[2].set_title("Fraction of samples with max reciprocal pairs θ_dgp = $θ_0 ")
    ax[2].set_xlabel( "η_dgp")
    ax[2].legend(Nvals)
    ax[2].grid()
    tight_layout()
end
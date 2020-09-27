
"""
Test script for dirBin0Rec0 model: one parameter for total number of links and one for reciprocity
"""

# test sampling

N=30
    θ_0 = -3.1
    η_0 = 2.7
    model = fooNetModelDirBin0Rec0
    diadProb = diadProbFromPars(model, [θ_0, η_0])
    A_vec = [samplSingMatCan(model, diadProb, N) for i=1:100]

hist([sum(A) for A in A_vec])
hist([sum(A'.*A) for A in A_vec])

A = A_vec[1]


# Test Estimates
par_est = reduce(hcat,[estimat(model, A) for A in A_vec])
hist(par_est[2,:], density=true )
 vlines(η_0, 0, 1)
hist(par_est[1,:], density=true )
 vlines(θ_0, 0, 1)

# Are likelihood and loglikelihod equal??

















#

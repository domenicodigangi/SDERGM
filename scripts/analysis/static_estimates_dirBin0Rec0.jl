
"""
Test script for dirBin0Rec0 model: one parameter for total number of links and one for reciprocity
"""

include("../../add_load_paths.jl")

using StaticNets
using ergmRcall
# test sampling

nSample = 500
N=30
    θ_0 = -3.1
    η_0 = 2
    model = fooNetModelDirBin0Rec0
    diadProb = diadProbFromPars(model, [θ_0, η_0])
    A_vec = [samplSingMatCan(model, diadProb, N) for i=1:nSample]

using PyPlot
hist([sum(A) for A in A_vec])
hist([sum(A'.*A) for A in A_vec])


# Test Estimates
par_est = reduce(hcat,[estimate(model, A) for A in A_vec])
hist(par_est[2,:], density=true )
 vlines(η_0, 0, 1, "r")
hist(par_est[1,:], density=true )
 vlines(θ_0, 0, 1, "r")

## Are loglikelihood and pseudo-loglikelihood equal??

formula_ergm = net ~ edges +  gwesp(decay = 0.25,fixed = TRUE)


# Write a funciton that, given a matrix returns the change statistics wrt to a given formula
install_req_R_packages()
clean_start_RCall()
T = size(parMatDgp_T)[2]
# For each t, sample the ERGM with parameters corresponding to the DGP at time t
# and run a single snapeshot estimate
@rput T; @rput parMatDgp_T;@rput N;@rput Nsample

reval("formula_ergm = net ~ "*formula_ergm_str)

    #create an empty network, the formula defining ergm, sample the ensemble and
    # store the sufficient statistics and change statistics in R
    R"
    net <- network.initialize(N)
    sampledMat_T_R =    array(0, dim=c(N,N,T,Nsample))
    changeStats_T_R = list()
    stats_T_R = list()
    for(t in 1:T){
        changeStats_t_R = list()
        stats_t_R = list()
        print(t)
        for(n in 1:Nsample){
                net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parMatDgp_T[,t],control = control.simulate.formula(MCMC.burnin = 100000))
                sampledMat_T_R[,,t,n] <- as.matrix.network( net)
                print(c(t,n,parMatDgp_T[,t]))
                chStat_t <- ergmMPLE(formula_ergm)
                changeStats_t_R[[n]] <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
                stats_t_R[[n]] <- summary(formula_ergm)
                }
                changeStats_T_R[[t]] <-changeStats_t_R
                stats_T_R[[t]] <- stats_t_R
            }"



# import sampled networks in julia
sampledMat_T = BitArray( @rget(sampledMat_T_R))
changeStats_T = @rget changeStats_T_R;# tmp = Array{Array{Float64,2}}(T); for t=1:T tmp[t] =  changeStats_T[t];end;changeStats_T = tmp
stats_T = @rget(stats_T_R); #tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = stats_T[t];end;stats_T = tmp


















#

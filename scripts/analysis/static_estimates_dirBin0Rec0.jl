
"""
Test script for dirBin0Rec0 model: one parameter for total number of links and one for reciprocity
"""

include("../../add_load_paths.jl")

using StaticNets
using ergmRcall
# test sampling

nSample = 1500
N=30
    θ_0 = -3.1
    η_0 = 2
    par_dgp = [θ_0, η_0] 
    model = fooNetModelDirBin0Rec0
    diadProb = diadProbFromPars(model, [θ_0, η_0])
    A_vec = [samplSingMatCan(model, diadProb, N) for i=1:nSample]

sum([sum(diag(A)) for A in A_vec])

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


# Write a funciton that, given a matrix returns the change statistics wrt to a given formula
ergmRcall.clean_start_RCall()
A = A_vec[1]
par = par_est[:,1]

using RCall


ergmTermsString = "edges +  mutual"
function get_change_stats(A::Matrix{T} where T<:Integer, ergmTermsString::String)
    @rput A
    reval("formula_ergm = net ~ "* ergmTermsString)
    # store the sufficient statistics and change statistics in R
    R"""
        net <- network(A)
        chStat_t <- ergmMPLE(formula_ergm)
        changeStats_t_R <- cbind(chStat_t$response, chStat_t$predictor,     chStat_t$weights)
            """

    changeStats = @rget changeStats_t_R;# tmp = Array{Array{Float64,2}}(T); for 
    return changeStats
end


changeStat = get_change_stats(A_vec[1],ergmTermsString) 


pseudo_loglikelihood(par_est[:,1], changeStat)


L,R = statsFromMat(model, A)

logLikelihood(model,L,R,N, par)

n=500
x_vals = LinRange(-0.2,0.2,n) 
par_vec_1 =  [par_dgp .+ [0, x] for x in x_vals]
par_vec_2 =  [par_dgp .+ [x, 0] for x in x_vals]



plot(x_vals,[logLikelihood(model, L,R,N, par) for par in par_vec_1])
plot(x_vals,[logLikelihood(model, L,R,N, par) for par in par_vec_2])







#

module AReg
using PolynomialRoots
using StatsBase

export ARp, ARp_est, fit , simulate,simulateAR1, unc_mean
##
"Autoregressive Structure (only normal errors are considered for the moment)"
mutable struct ARp
    vecpar::Vector{Real}
    sigma::Real
    start_vals::Vector{Real}
end

mutable struct ARp_est
    ar::ARp
    obs::Vector{Real}
    estimates::Vector{Real}
end


function unc_mean(arp::ARp)
    "compute the unconditional mean of an ARp type"
    p = length(arp.vecpar)
    p==0 ? (return 0) : (return arp.vecpar[1]/(1-sum(arp.vecpar[2:end])))
end
function is_ar_stationary(arp::ARp)
    coeff = arp.vecpar
    #coefficients of the characteristic polynomial
    char_pol_coeff = [coeff[1]; - coeff[2:end]]
    #compute the roots of the char poly
    tmp =abs2.(roots(Array{Float64}(char_pol_coeff)))
    #if all the roots lie outside the unit circle then the process is stationary
    all(tmp.>1) ? (return true,unc_mean(arp)) : (return false)
end



##
"fit an AR(p) with predetermined p, the starting values are not estimated
To Do : need to add a method for the observations only.. it should return the
ARp_est with optimal p
"
function fitARp(obs::Vector{<:Real},arp::ARp;demean=true)
    p = length(arp.vecpar)-1 # number of autoregressive terms (constant excluded)
    ## Ols estimation of the AR(p) coefficients
    Y = obs[p+1:end]
    T = length(Y)
    X = ones(T,p+1)
    #println(p)
    for i=1:p
        X[:,p+2-i] = obs[p+1-i:end-i]
    end

    est_par = (X'*X)\(X'Y)
    #est_par = (X'*X)\(X'Y)
    return ARp_est(arp,obs,est_par)
end
##

"Sample an AR process of lenght T and  parameters specified in vecpar starting
from phi_0 to phi_p where p is the order of the process"
function simulate(ar::ARp,T::Int=1000)
    sigma = ar.sigma
    vecpar = ar.vecpar
    start_vals = ar.start_vals
    #Y_T = zeros(T)
    p = length(vecpar)-1
    epsilon_T = sigma.* randn(T)
    Y_T = epsilon_T

    if isempty(start_vals) # if no start values are given pick them at ranUtilities
        if p>0
            start_vals = sigma.*randn(p)
        else
            start_vals = []
        end
    end
    lagged_vals = start_vals
    for t = 1:T
        Y_T[t] +=  vecpar'*[1;lagged_vals ]
        if p>0
            next_lagged_vals = circshift(lagged_vals,-1)
            next_lagged_vals[end] = Y_T[t]
        else
            next_lagged_vals = Y_T[t]
        end
        lagged_vals = next_lagged_vals
    end
    return Y_T
end
function simulateAR1(vecPar::Vector{<:Real},sigma::Float64;T::Int=1000)
    #start on the unconditional mean
    μ = vecPar[1]./(1-vecPar[2])
    #Y_T = zeros(T)
    epsilon_T = sigma.* randn(T)
    Y_T = epsilon_T
    lagged_val = μ
    for t = 1:T
        Y_T[t] +=  vecPar[1] + vecPar[2] * lagged_val
        lagged_val = Y_T[t]
    end
    return Y_T
end


end

using RCall
ErgmRcall.clean_start_RCall()
R"""options(warn=-1) """

#---------------------------------Binary DIRECTED Networks with one Reciprocity par

struct  NetModelDirBin0Rec0 <: NetModelBin #Bin stands for (Binary) Adjacency matrix
    "ERGM for directed networks with totan links and tortal reciprocated links as statistics "
    obs:: Array{<:Real,1} # [total links, total reciprocated links]
    N::Int
    vPar::Array{<:Real,1} # One parameter per statistic
end
fooNetModelDirBin0Rec0 =  NetModelDirBin0Rec0(ones(2), 10, zeros(2))
export fooNetModelDirBin0Rec0
NetModelDirBin0Rec0(A:: Matrix) =
    NetModelDirBin0Rec0(obs, size(A)[1], zeros(2) )

function identify(Model::NetModelDirBin0Rec0,par::Array{<:Real,1};idType ="" )
    return par
  end

function diadProbFromPars(Model::NetModelDirBin0Rec0, par)
    """
    givent the 2model's parameters compute the probabilities for each state of the diad
    """
    θ = par[1]
    η = par[2]
    eθ = exp(θ)
    e2θplusη = exp(2*θ + η)
    den = 1 + 2*eθ + e2θplusη
    p00 = 1/den
    p01 = eθ /den
    p10 = p01 
    p11 = e2θplusη/den
    diadProbsVec = [p00, p01, p10, p11]
    return diadProbsVec
end


"""
given the vector of diad states probabilities sample one random matrix from the corresponding pdf
"""
function samplSingMatCan(Model::NetModelDirBin0Rec0, diadProbsVec::Array{<:Real,1}, N)
    out = zeros(Int8,N,N)
    #display((maximum(expMat),minimum(expMat)))
    for c = 1:N
        for r=1:c-1
            diadState = rand(Categorical(diadProbsVec))
            if diadState==2
                out[r,c] = 1
            elseif diadState==3
                out[c,r] = 1
            elseif diadState==4
                out[r,c] = 1
                out[c,r] = 1
            end
        end
    end
    return out
 end


function sample_ergm(model::NetModelDirBin0Rec0, N, θ_0, η_0, nSample)

    diadProb = diadProbFromPars(model, [θ_0, η_0])

    A_vec = [samplSingMatCan(model, diadProb, N) for i=1:nSample]

    return A_vec
end



function statsFromMat(Model::NetModelDirBin0Rec0, A ::Matrix{<:Real})
    L = sum(A)
    R = sum(A'.*A)/2
    N=size(A)[1]

return [L, R, N]
end

function exp_val_stats(model::NetModelDirBin0Rec0, θ, η, N)
    x = exp(θ)
    x2y = exp(2*θ + η)
    z = 1 + 2*x + x2y
    Nlinks = n_pox_dir_links(N)
    α = (x + x2y)/z
    β = 0.5 * (x2y)/z
    L = Nlinks * α
    R = Nlinks * β
    return L, R
end

"""
return the ergm parameters that fix on average the input values for the ergm statistics:
    - L is  ∑_i>j A_ij + A_ji
    - R is  of ∑_i>j A_ij * A_ji /2
"""
function ergm_par_from_mean_vals(model::NetModelDirBin0Rec0, L, R, N)
    Nlinks = n_pox_dir_links(N)
    # average values per pair
    α = L/Nlinks
    β = R/Nlinks

    x = (2β-α)/(2α-2β-1)
    θ = log(x)

    y = 2(2*β^2  + β- 2*α*β )/((2β-α)^2)
    η =  log(y)

    return θ, η 
end



function estimate(model::NetModelDirBin0Rec0, L, R, N)
    # α = L / (N*(N-1))
    # β = R / (N*(N-1)) 
    # tmp = 1+2*β-2*α
    # x = (α - 2*β)/tmp
    # # y = (α + (2α-1)*x )/(x^2)
    # y = (α + (2*α-1)*x )/((1-α)*x^2)
    # # y = (2*β*(1-2*x))/(x^2*(1-2*β))
    # # y = 1 + (2β-4β*α-8β^2*α+4β*α-α^2)
    # θ_est = log(x)
    # η_est = log(y)
    # # η_est = log( 2*β*tmp*(1-α) /((α-2*β)^2)   )
   
    θ_est, η_est = ergm_par_from_mean_vals(model, L, R, N)
    vPar = [θ_est, η_est]
    return vPar
end



function estimate(Model::NetModelDirBin0Rec0, A)
    L, R, N = statsFromMat(Model, A) 
    return estimate(Model, L, R, N)
end





function logLikelihood(Model::NetModelDirBin0Rec0, L, R, N, par)
    θ, η = par
    z = 1 + 2*exp(θ) + exp(2*θ+η)
    return L * θ + R*η - (N*(N-1)/2)*log(z)
end

function logLikelihood(Model::NetModelDirBin0Rec0, A::Matrix, par)
    L, R, N = statsFromMat(Model, A)
    return logLikelihood(Model, L, R, N, par)
end

function change_stats(Model::NetModelDirBin0Rec0, A::Matrix)
    ergmTermsString = "edges +  mutual"
    return ErgmRcall.get_change_stats(A,ergmTermsString)
end




# pseudolikelihood function from change statistics following strauss and ikeda
function pseudo_loglikelihood_strauss_ikeda(Model::NetModelDirBin0Rec0, par, changeStat, response, weights)
    logit_P = sum(par.*changeStat', dims=1)      
    P = inv_logit.(logit_P)    
    logPVec = log.([response[i] == zero(response[i]) ? 1 - P[i] : P[i] for i=1:length(response) ])
    logPTot = sum(logPVec.*weights)
    return  logPTot
end

function pseudo_loglikelihood_from_sdergm(Model::NetModelDirBin0Rec0, par::Array{<:Real,1}, changeStat::Array{<:Real,2})
    indPres = (changeStat[:,1]).>0 # change stats of matrix elements that are 
    tmpMatPar = par' .* changeStat[:,2:end-1]
    p_ij =  exp.( .- sum(tmpMatPar,dims = 2))
    mult =changeStat[:,end]
    logpseudolike_t =  .-  sum(mult .* log.(1 .+ 1 ./ p_ij ))
    π = (1 ./ (1 .+ p_ij  ))
    for p =1:length(par)
        δ_p = changeStat[:,p+1] # p+1 because to skip the first column
        tmp1 = sum( mult[indPres] .* δ_p[indPres] )
        logpseudolike_t += par[p]*tmp1
    end
    return logpseudolike_t
end


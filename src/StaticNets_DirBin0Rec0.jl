
#---------------------------------Binary DIRECTED Networks with one Reciprocity par

struct  NetModelDirBin0Rec0 <: NetModelBin #Bin stands for (Binary) Adjacency matrix
    "ERGM for directed networks with totan links and tortal reciprocated links as statistics "
    obs:: Array{<:Real,1} # [total links, total reciprocated links]
    N::Int
    vPar::Array{<:Real,1} # One parameter per statistic
end
fooNetModelDirBin0Rec0 =  NetModelDirBin0Rec0(ones(2), 10, zeros(2))
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
    e2θplus2η = exp(2*(θ + η))
    den = 1 + 2*eθ + e2θplus2η
    p00 = 1/den
    p01 = eθ /den
    p11 = e2θplus2η/den
    diadProbsVec = [p00, p01, p01, p11]
    return diadProbsVec
end

function samplSingMatCan(Model::NetModelDirBin0Rec0, diadProbsVec::Array{<:Real,1}, N)
    """
    given the vector of diad states probabilities sample one random matrix from the corresponding pdf
    """
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

function statsFromMat(Model::NetModelDirBin0Rec0, A ::Matrix{<:Real})
    L = sum(A)
    R = sum(A'.*A)/2

return [L, R]
end

function estimate(Model::NetModelDirBin0Rec0, L, R, N )
    L_bar = L / (N*(N-1))
    R_bar = 2*R / (N*(N-1)) # when I coputed it I defined R^bar as twice the number of average mutual pairs. This is the only reason for the two in front
    θ_est = log((R_bar - L_bar)/(2*L_bar - R_bar -1) )
    η_est = log( (- R_bar)*(2*L_bar - R_bar -1)/((1- R_bar)*(R_bar-L_bar)^2)  )/2
    vPar = [θ_est, η_est]
    return vPar
end

function estimate(Model::NetModelDirBin0Rec0, A ::Matrix{<:Real})
    N = size(A)[1]
    L, R = statsFromMat(Model, A)
    return estimate(Model, L, R, N )
end

function logLikelihood(Model::NetModelDirBin0Rec0, L, R, N, par)
    θ, η = par
    return L * θ + R*η -  (N*(N-1)/2)*log(1 + 2*exp(θ) + exp(2*θ + η))
end

function pseudoLogLikelihood(Model::NetModelDirBin0Rec0, L, R, par)
    θ, η = par
    return L * θ + R*η - (N*(N-1))*log(1 + 2*exp(θ) + exp(2*(θ + η)))
end

function pseudo_loglikelihood(par::Array{<:Real,1}, changeStat::Array{<:Real,2})
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
export pseudo_loglikelihood
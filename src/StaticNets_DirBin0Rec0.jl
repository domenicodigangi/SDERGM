
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
    e2θplusη = exp(2*θ + η)
    den = 1 + 2*eθ + e2θplusη
    p00 = 1/den
    p01 = eθ /den
    p10 = p01 
    p11 = e2θplusη/den
    diadProbsVec = [p00, p01, p10, p11]
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

function estimate(Model::NetModelDirBin0Rec0, A)
    L, R = statsFromMat(Model, A) 
    N = size(A)[1]
    L_bar = L / (N*(N-1))
    R_bar = R / (N*(N-1)) 

    tmp = 1+2*R_bar-2*L_bar
    x = (L_bar - 2*R_bar)/tmp
    # y = (L_bar + (2L_bar-1)*x )/(x^2)
    y = (L_bar + (2*L_bar-1)*x )/((1-L_bar)*x^2)
    # y = (2*R_bar*(1-2*x))/(x^2*(1-2*R_bar))
    # y = 1 + (2R_bar-4R_bar*L_bar-8R_bar^2*L_bar+4R_bar*L_bar-L_bar^2)
    θ_est = log(x)
    η_est = log(y)
    # η_est = log( 2*R_bar*tmp*(1-L_bar) /((L_bar-2*R_bar)^2)   )
    vPar = [θ_est, η_est]
    return vPar
end






function pseudo_loglikelihood_from_sdergm(par::Array{<:Real,1}, changeStat::Array{<:Real,2})
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
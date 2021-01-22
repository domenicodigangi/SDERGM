using RCall
ErgmRcall.clean_start_RCall()
R"""options(warn=-1) """


"""
A generic ERGM for directed networks to be estimated using the pseudo likelihood 
"""
struct  NetModelDirBinErgmPml <: NetModelDirBin 
    ergmTermsString::String # needs to be compatible with R ergm package notation and names
    nErgmPar::Int
end
export NetModelDirBinErgmPml

NetModelDirBinErgmPml(ergmTermsString)  = NetModelDirBinErgmPml(fixes, 1 + count(i->(i=='+'), "fixes"))

name(x::NetModelDirBinErgmPml) = "NetModelDirBinErgmPml($(x.ergmTermsString))"

function change_stats(model::NetModelDirBinErgmPml, A::Matrix)
    ErgmRcall.get_change_stats(A,model.ergmTermsString)
end


function pseudo_loglikelihood_strauss_ikeda(model::NetModelDirBinErgmPml, par, changeStat, response, weights)
    logit_P = sum(par.*changeStat', dims=1)      
    P = inv_logit.(logit_P)    
    logPVec = log.([response[i] == zero(response[i]) ? 1 - P[i] : P[i] for i=1:length(response) ])
    logPTot = sum(logPVec.*weights)
    return  logPTot
end


estimate(model::NetModelDirBinErgmPml, A::Matrix) = ErgmRcall.get_mple(model.ergmTermsString, A)

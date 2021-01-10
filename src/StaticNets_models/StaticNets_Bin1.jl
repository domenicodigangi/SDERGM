
##------------------ Binary UNDIRECTED Networks

struct  NetModelBin1 <: NetModelBin #Bin stands for (Binary) Adjacency matrix
    "ERGM with degrees as statistics and possibility to have groups of nodes
    associated with the same parameter. "
    obs:: Array{<:Real,1} #[degrees]
    Par::Array{<:Real,1} # One parameter per Node  UNRestricted
    groupsInds::Array{<:Real,1} #The index of the groups each node belongs to
 end
fooNetModelBin1 =  NetModelBin1(ones(3),zeros(3),ones(3))
NetModelBin1(deg:: Array{<:Real,1}) =  NetModelBin1(deg,zeros(length(deg)),Vector(1:length(deg)) )
function DegSeq2graphDegSeq(Model::NetModelBin1,Seq::Array{<:Real,1})
    "From a degree sequence find another one that is certanly graphicable.
     I choose a naive way to ensure graphicability.
      The matrix associated with the starting degree sequence is very hubs oriented
      the largest node catches as many links as it needs to fill all its strubs "

    Seq = sort(Seq,rev=true)#decreasing sort
    N = length(Seq)
    stubs = ceil.(Int,copy(Seq))
    mat = zeros(Int64,N,N)
    for i in N:-1:1
       j = 1
       while (stubs[i]>0)&(j<=N)
           if (mat[i,j] == 0)& (i!=j)& (stubs[j]>0)
               mat[i,j] = 1
               mat[j,i] = 1
               stubs[i] -= 1
               stubs[j] -= 1
               #println(j)
           else
               #println((i, j,mat[i,j],stubs[j]))
           end
           j += 1
       end
    end
    #display(stubs)
    #display(mat)
    ((sum(mat[Array{Bool}(eye(mat))])) != 0)  ?  error() : ()
    graphSeq= sumSq(mat,2)
    return graphSeq
 end
function linSpacedPar(Model::NetModelBin1,Nnodes::Int; Ngroups = Nnodes,deltaN::Int=3,graphConstr = true)
    "Generate an equally spaced degree sequence. Ensure that it is graphicable
    and then estimate the parameters "

          #generate syntetic degree distributions, equally spaced, and estimate the model
    if graphConstr
        tmpdegs  = DegSeq2graphDegSeq(fooNetModelBin1, Array{Real,1}(range(1+deltaN,stop=Nnodes - round(0.1*Nnodes)-deltaN, length=Nnodes)))
    else
        #stima modello a Ngroups fitness costanti con degree dei nodi linearly spaced
        tmpdegs  =   Array{Real,1}(range(deltaN,stop=Nnodes-deltaN, length=Nnodes) )
        tmpdegs = sort(tmpdegs,rev =true)
    end
    #if Nnodes==Ngroups
    groupsInds = distributeAinVecN(Vector(1:Ngroups),Nnodes)
    out,it,estMod = estimate(NetModelBin1(tmpdegs, zeros(length(tmpdegs)), groupsInds))
    #else
    return out,tmpdegs
 end
function samplSingMatCan(Model::NetModelBin1,expMat::Array{<:Real,2})
    """
    Being p a matrix of entries in [0,1], samples one bernoulli r.v. for eache matrix element
    """
    N1,N2 = size(expMat)
    N1!=N2  ?  error() : N=N1
    out = zeros(Int8,N,N)
    for c = 1:N
        for r=1:c-1
            out[r,c] = rand(Bernoulli(expMat[r,c]))
        end
    end
    out = Symmetric(out)
    return out
 end

function expMatrix2(Model::NetModelBin1,par::Array{<:Real,1})
    "Given the vector of model parameters, return the product matrix(often useful
    in likelihood computation) and the expected matrix"
    parMat = putZeroDiag( exp.(Symmetric(par  .+ par') ))
    expMat = parMat ./ (1 .+ parMat)
    return parMat, expMat
 end
expMatrix(Model::NetModelBin1,par::Array{<:Real,1}) = expMatrix2(Model,par)[2] # return only the expected matrix given the parameters
function expValStatsFromMat(Model::NetModelBin1,expMat::Array{<:Real,2}  )
    "Expected value of the statistic for a model, given the expected matrix"
    expDeg = sumSq(expMat,2)
    expValStats = expDeg
    return expValStats
 end
function firstOrderCond(Model::NetModelBin1;deg::Array{<:Real,1} = Model.obs, par::Array{<:Real,1}=Model.Par, groupsInds::Array{<:Real,1} = Model.groupsInds )
    "Given the model, the degree sequence, the parameters and groups assignements,
    return the First order conditions. Gradient of the loglikelihood, or system
    of differences between expected and observed values of the statistics (either
    degrees of single nodes or degree of a group). "
    N = length(deg)
    groupsNames = unique(groupsInds)
    NG = length(groupsNames)

    if N==NG
        expMat = expMatrix(Model,par)
        expDeg = sumSq(expMat,2)
        foc = deg - expDeg
    elseif N<NG
        error()
    else
        tmpDic = countmap(groupsInds)
        # cumulative degree of each group
        degGroups = zeros(NG) ; for i =1:N degGroups[groupsInds[i]] += deg[i] end
        nodesPar = par[groupsInds]
        expMat = expMatrix(Model,nodesPar)
        #display(expMat)
        expDeg = sumSq(expMat,2)
        expDegGroups = zeros(NG) ; for i =1:N expDegGroups[groupsInds[i]] += expDeg[i] end
        foc = degGroups - expDegGroups
    end
    return foc
 end
function identify(Model::NetModelBin1,Par::Array{<:Real,1})
    "Given a vector of parameters, return the transformed vector that verifies
    an identification condition. Do not do anything if the model is identified."
    return Par
  end
function estimate(Model::NetModelBin1; deg::Array{<:Real,1} = Model.obs,
            groupsInds = Model.groupsInds, targetErr::Real=targetErrValStaticNets,bigConst = bigConstVal)
    "Given model type, observations of the statistics and groups assignments,
    estimate the parameters."
    # The parameters inside this function are Restricted to be in [0,1] hence
    # they need to be converted (log(par)) in order to be used outside
    N = length(deg)
    groupsNames = unique(groupsInds)
    NG = length(groupsNames)

    tmpDic = countmap(groupsInds)
    #number of nodes per group
    NgroupsMemb = sortrows([[i for i in keys(tmpDic)] [  tmpDic[i] for i in keys(tmpDic) ] ])
    # agrregate degree of each group
    degGroups = zeros(NG) ; for i =1:N degGroups[groupsInds[i]] += deg[i] end

    #sort the groups for degrees
    indsDegGroupsOrder = Int.(sortrows([degGroups Vector(1:NG)], rev = true)[:,2])

    Nzer = sum(degGroups.==0)
    L = sum(deg)
    LperLink  = L/(N*N-N)
    #println((LperLink,N,L))
    #uniform parameter as starting value
    unifstartval =  sqrt((LperLink/((1-LperLink))));
    #analytically solvable case
    if NG==1
        outParVec = unifstartval*ones(NG)
        outMod = NetModelBin1(Model.obs,outParVec,groupsInds)
        i = 0
    else
        par = unifstartval*ones(NG)
        eps = 1e-8
        i=0
        maxIt = 150
        relErr = 1
        tol_x = 10*Base.eps()
        function tmpFindZero(objFun::Function,brak::Array{<:Real,1},start::Real)
        #    wrapper find_zero options
            #println((brak[1],brak[2]))
            #println((objFun(brak[1]),objFun(brak[2])))
            out = find_zero(objFun,(brak[1],brak[2]),Bisection())
            #println((out,objFun(out)))
            return out
        end
        while (i<maxIt)&(relErr>targetErr)
            parMax = maximum(par)
            parMin = minimum(par[par.!=0])
            for nInd=1:NG
                #solve first the gradient of nodes with largest degrees
                ng = indsDegGroupsOrder[nInd]
                #put to zero each component of the gradient1
                function objfun_ng(par_ng::Real)
                     parVecOthGroups = par_ng.*par[1:end .!=ng]
                     expGroupDeg = NgroupsMemb[ng,2]*( sum( NgroupsMemb[1:end .!=ng ,2].*( parVecOthGroups ./ (1 + parVecOthGroups)) )    + (NgroupsMemb[ng,2]-1)* (par_ng^2)/(1+(par_ng^2)) )
                     return  degGroups[ng] - expGroupDeg
                end
                #find zero of given component

                if degGroups[ng]== (N-1-Nzer)*NgroupsMemb[ng,2]

                    par[ng] = (1/parMin)*(1-0.1*targetErr)/( 0.1*targetErr)
                elseif degGroups[ng]==0
                    par[ng] = 0
                else
                     tmpBnd = (degGroups[ng]/(NgroupsMemb[ng,2]*(N-1-Nzer)) )/(1-degGroups[ng]/(NgroupsMemb[ng,2]*(N-1-Nzer)) )

                    brakets = [ 0,2* tmpBnd/parMin]
                    par[ng] = tmpFindZero(objfun_ng,brakets,par[ng])
                end
                #println(par)
                #parUN = log.(par)
                #display(firstOrderCond(Model;deg = deg,par = parUN, groupsInds = groupsInds )')

            end
            # The parameters inside this function are Restricted to be in [0,1] hence
            # they need to be converted (log(par)) in order to be used outside
            parUN = log.(par)
            parUN[par.==0] = -bigConstVal
            err = firstOrderCond(Model;deg = deg,par = parUN, groupsInds = groupsInds )
            #    tmpMat = par.*par';expMat = putZeroDiag!(tmpMat ./ (1+tmpMat)) ;expDeg = sumSq(expMat,2)
            relErrVec = abs.(err ./ degGroups)
            relErr =  maximum(relErrVec[degGroups.!=0])
            # display(err')
             #println(relErr)
            i+=1
         end
         if i<maxIt
             # do nothing, optimization succeded
         else
             println(deg)
             if NG==N
                 error("with one node per group the optimization should succeed")
             else
                 error("Add a maximum likelihood refinement for the groups case if the map does not succeeed, as it might very well happen")
             end
         end
        parUN = log.(par)
        parUN[par.==0] = -bigConst
        outParVec = parUN
        outMod = NetModelBin1(Model.obs,outParVec,groupsInds)
    end
    return outParVec, i , outMod
 end
function sampl(Model::NetModelBin1,Nsample::Int;  par::Array{<:Real,1}=Model.Par ,groupsInds = Model.groupsInds)
    "Sample a given model."
    nodesPar = par[groupsInds]
    N = length(nodesPar)
    expMat = expMatrix(Model,nodesPar)
    SampleMats =zeros(Int8,N,N,Nsample)
    for s=1:Nsample
        SampleMats[:,:,s] = samplSingMatCan(Model,expMat)
    end
    return SampleMats
 end

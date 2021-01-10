##
#Binary DIRECTED Networks
##

struct  NetModelDirBin1 <: NetModelBin #Bin stands for (Binary) Adjacency matrix
    "ERGM for directed networks with in and out degrees as statistics and
    possibility to have groups of nodes associated with a single pair of in out
    parameters. "
    obs:: Array{<:Real,1} #[InDegs;OutDegs] degrees of each node (even if groups are present)
    Par::Array{<:Real,1} # One parameter per group of nodes  UNRESTRICTED VERSION  [InPar;OutPar]
    groupsInds::Array{<:Real,1} #The index of the groups each node belongs to
    # two parameters per node (In and out)
 end
fooNetModelDirBin1 =  NetModelDirBin1(ones(6),zeros(6),ones(3))
NetModelDirBin1(deg:: Array{<:Real,1}) =  NetModelDirBin1(deg,zeros(length(deg[:])) , Int.(1:round(Int,length(deg[:])/2)))

function identify(Model::NetModelDirBin1,parIO::Array{<:Real,1};idType ="equalIOsums" )
    "Given a vector of parameters, return the transformed vector that verifies
    an identification condition. Do not do anything if the model is identified."
    #set the first of the in parameters equal to one (Restricted version)
    N,parI,parO = splitVec(parIO)

    if idType == "equalIOsums"
        Δ = sum(parI[isfinite.(parI)]) - sum(parO[isfinite.(parO)])
        shift = Δ/(2N)
    elseif idType == "firstZero"
        shift = parI[1]
    end
    parIO = [ parI .- shift ;  parO .+ shift ]
    return parIO
  end


function DegSeq2graphDegSeq(Model::NetModelDirBin1,SeqIO::Array{<:Real,1})
    #From a degree sequence find another one that is certanly graphicable
    #The matrix associated with the starting degree sequence is very hubs oriented

    # the largest node catches as many links as it needs to fill all its strubs
    N, SeqI,SeqO = splitVec(SeqIO)
    sum(SeqI) == sum(SeqO) ? () : error()
    inds = Vector(1:N)
    StubsI_inds = Int.(sortrows([SeqI inds],rev=true))#decreasing sort
    StubsO_inds = Int.(sortrows([SeqO inds],rev=true))

    mat = zeros(Int64,N,N)
    for i in 1:N
       j = 1
       while (StubsI_inds[i,1]>0)&(j<N)
           ind_i = StubsI_inds[i,2]
           ind_j = StubsO_inds[j,2]

           if (mat[ind_i,ind_j] == 0)& (ind_i!=ind_j) &(StubsO_inds[j,1]>0)

               mat[ind_i,ind_j] = 1
               StubsI_inds[i,1]  -= 1
               StubsO_inds[j,1]  -= 1
           end
           j += 1
       end
    end
    ((sum(mat[Array{Bool}(eye(mat))])) != 0)  ? error() : ()
    graphSeq= [sumSq(mat,2); sumSq(mat,1)]
    return graphSeq
 end
function zeroDegParFun(Model::NetModelDirBin1,N::Int;degsIO=zeros(Int,10),
            parIO=zeros(Float64,10), method="AVGSPACING" ,smallProb = 1e-1)
    #define a number big enough to play the role of Inf for
    # purposes of sampling N(N-1) bernoulli rvs with prob 1/(1+exp(  bigNumb))
    # prob of sampling degree > 0 (N-1)/(1+exp(bigNumb)) < 1e-6
    if uppercase(method) == "INF"
        # zero degrees have -Inf fitnesses
        zeroDegPar = - Inf #log((N-1)/smallProb - 1)
    elseif uppercase(method) == "SMALL"
        # such that the probability of observing a link is very small
        zeroDegPar = -10
    elseif uppercase(method) == "AVGSPACING"
        if (sum(degsIO)==0)|(sum(parIO) ==0)
            error()
        end
        ~,degI ,degO = splitVec(degsIO)
        ~,parI ,parO = splitVec(parIO)
        imin = findmin(degI[degI.!=0])[2]
        omin = findmin(degO[degO.!=0])[2]
        # find a proxy of the one degree difference between fitnesses
        diffDegs = degI[degI.!=0] .- degI[degI.!=0][imin]
        avgParStepI = mean(((parI[degI.!=0] .- parI[degI.!=0][imin])[diffDegs.!=0]) ./ diffDegs[diffDegs.!=0])
        !isfinite(avgParStepI) ? println(((parI .- parI[imin])[diffDegs.!=0]) ./ diffDegs[diffDegs.!=0])error() : ()
        diffDegs = degO[degO.!=0] .- degO[degO.!=0][omin]
        avgParStepO = mean((parO[degO.!=0] .- parO[degO.!=0][omin])[diffDegs.!=0] ./ diffDegs[diffDegs.!=0])
        !isfinite(avgParStepO) ? println(((parI .- parI[imin])[diffDegs.!=0]) ./ diffDegs[diffDegs.!=0])error() : ()
        #subtract a multiple of this average step to fitnesses of smallest node
        zeroDegParI = parI[degI.!=0][imin] .- avgParStepI .* degI[degI.!=0][imin]
        zeroDegParO = parO[degO.!=0][omin] .- avgParStepO .* degO[degO.!=0][omin]
    end
    return zeroDegParI,zeroDegParO
 end
function linSpacedPar(Model::NetModelDirBin1,Nnodes::Int;Ngroups = Nnodes,deltaN::Int=3,graphConstr = true,degMax = Nnodes-2)
    "for a given static model return a widely spaced sequence of groups parameters"
    N=Nnodes
    #generate syntetic degree distributions, equally spaced, and estimate the model
    tmpdegsIO = round.([ Vector(range(deltaN,stop=degMax,length=N));Vector(range(deltaN,stop=degMax,length=N))])
    if graphConstr
        tmpdegsIO  =   DegSeq2graphDegSeq(fooNetModelDirBin1,tmpdegsIO)
    end

    if Ngroups>N#if more groups than nodes are required throw an error
        error()
    else
        groupsInds = distributeAinVecN(Vector(1:Ngroups),Nnodes)
    end
    out,it,estMod = estimate(NetModelDirBin1(tmpdegsIO))
    return out,tmpdegsIO
 end
function samplSingMatCan(Model::NetModelDirBin1,expMat::Array{<:Real,2})
    """
    given the expected matrix sample one ranUtilities matrix from the corresponding pdf
    """
    N1,N2 = size(expMat)
    N1!=N2 ? error() : N=N1
    out = zeros(Int8,N,N)
    #display((maximum(expMat),minimum(expMat)))
    for c = 1:N
        for r=1:N
            if r!=c
                out[r,c] = rand(Bernoulli(expMat[r,c]))
            end
        end
    end
    return out
 end
function expMatrix2(Model::NetModelDirBin1, parNodesIO::Array{<:Real,1}  )
    "Given the vector of model parameters (groups), return the product matrix(often useful
    in likelihood computation) and the expected matrix"
    N,parI,parO = splitVec(parNodesIO)
    parMat =  Utilities.putZeroDiag_no_mut(exp.(parI  .+ parO') )
    expMat = parMat ./ (1 .+ parMat)
    return parMat,expMat
 end
expMatrix(Model::NetModelDirBin1, parNodesIO::Array{<:Real,1}  ) = expMatrix2(Model, parNodesIO )[2]# return only the expected matrix given the parameters
function expValStatsFromMat(Model::NetModelDirBin1,expMat::Array{<:Real,2}  )
    "Expected value of the statistic for a model, given the expected matrix"
    expDeg = [sumSq(expMat,2);sumSq(expMat,1)]
    expValStats = expDeg
    return expValStats
 end
expValStats(Model::NetModelDirBin1, parNodesIO::Array{<:Real,1} ) = expValStatsFromMat(Model, expMatrix(Model,parNodesIO )) # if only the paramters are given as input
function firstOrderCond(Model::NetModelDirBin1;degIO::Array{<:Real,1} = Model.obs,
                    degGroups::Array{<:Real,1} = zeros(4), parGroupsIO::Array{<:Real,1}=Model.Par, groupsInds::Array{<:Real,1} = Model.groupsInds )
    "Given the model, the degree sequence, the parameters and groups assignements,
    return the First order conditions. Gradient of the loglikelihood, or system
    of differences between expected and observed values of the statistics (either
    degrees(In and Out) of single nodes or of a group). "
    N,degI,degO = splitVec(degIO)
    groupsNames = unique(groupsInds)
    NG = length(groupsNames)
    NGcheck,parI,parO = splitVec(parGroupsIO)
    NG==NGcheck ?  () : error()
    if sum(degGroups)==0
        degGroups = zeros(2NG) ; for i =1:N degGroups[groupsInds[i]] += degI[i];degGroups[NG + groupsInds[i]] += degO[i]; end
    end
    parNodesIO = [parI[groupsInds];parO[groupsInds]]
    N,expDegI,expDegO = splitVec(expValStats(Model,parNodesIO))
    expGroupsDeg = zeros(2NG) ; for i =1:N expGroupsDeg[groupsInds[i]] += expDegI[i];expGroupsDeg[NG + groupsInds[i]] += expDegO[i]; end
    #println(groupsInds)
    #println(degGroups)
    #println(expGroupsDeg)
    foc = degGroups .- expGroupsDeg
    return foc
 end


function estimateSlow(Model::NetModelDirBin1; degIO::Array{<:Real,1} = Model.obs,
                    groupsInds = Int.(1:round(length(degIO)/2)), targetErr::Real=targetErrValStaticNets,
                    startVals =zeros(2))
    "Given model type, observations of the statistics and groups assignments,
    estimate the parameters."
    N,degI,degO = splitVec(degIO)
    groupsNames = unique(groupsInds)
    NG = length(groupsNames)
    tmpDic = countmap(groupsInds)
    #number of nodes per group
    NgroupsMemb = sortrows([[i for i in keys(tmpDic)] [  tmpDic[i] for i in keys(tmpDic) ] ])
    # agrregate degree of each group
    degGroups = zeros(2NG) ; for i =1:N degGroups[groupsInds[i]] += degI[i]; degGroups[NG + groupsInds[i]] += degO[i]; end
    ~,degGroupsI , degGroupsO = splitVec(degGroups)

    #sort the groups for degrees
    indsDegGroupsOrderI = Int.(sortrows([degGroupsI Vector(1:NG)], rev = true)[:,2])
    indsDegGroupsOrderO = Int.(sortrows([degGroupsO Vector(1:NG)], rev = true)[:,2])

    NzerI = sum(degGroupsI.==0)
    NzerO = sum(degGroupsO.==0)
    (sum(degGroupsI) - sum(degGroupsO) )< 10^Base.eps() ?  L =sum(degGroupsO) : error()
    LperLink  = L/(N*N-N)
    # The parameters inside this function are Restricted to be in [0,1] hence
    # they need to be converted (log(par)) in order to be used outside

    #uniform parameter as starting value
    unifstartval =  sqrt((LperLink/((1-LperLink))));
    if NG==1
        outParUN = [ 0 ; log(unifstartval)]
        outMod = NetModelDirBin1(Model.obs,outParUN,groupsInds)
        i = 0
    else
        parI = unifstartval*ones(NG)
        parO = unifstartval*ones(NG)
        #bigNumber = 1e3 ./ targetErr
        eps =  targetErr*1e-12 #10*Base.eps()

        i=0
        maxIt = 25
        relErr = 1
        tol_x = 10*Base.eps()
        sum(startVals)==0 ? () : (parI=startVals[1:N];parO=startVals[1+N:end];)
        function tmpFindZero(objFun::Function,brak::Array{<:Real,1},start::Real)
            #    wrapper find_zero options
            #println((brak[1],brak[2]))
            #println((objFun(brak[1]),objFun(brak[2])))
            #println()
            out = find_zero(objFun,(brak[1],brak[2]),Bisection())
            #println((out,objFun(out)))
            return out
        end
        while (i<maxIt)&(relErr>targetErr)
            parOMax = maximum(parO)
            parOMin = minimum(parO[parO.!=0])
            for nInd=1:NG
                ng = indsDegGroupsOrderI[nInd]
            #    println((1,n))
                #put to zero each component of the gradient1
                function objfunI_ng(parI_ng::Real)
                    parVecOthGroups = parI_ng.*parO[1:end .!=ng]
                    parO_ng = parO[ng]
                    expGroupDegI = NgroupsMemb[ng,2]*( sum( NgroupsMemb[1:end .!=ng ,2].*( parVecOthGroups ./ (1 + parVecOthGroups)) )    + (NgroupsMemb[ng,2]-1)* (parI_ng*parO_ng)/(1+(parI_ng*parO_ng)) )
                    return  degGroupsI[ng] - expGroupDegI

                end
                #find zero of given component

                if degGroupsI[ng] == NgroupsMemb[ng,2]*(N-1-NzerO)
                    #if a group has the maximum possible in degree set par to a
                    #very high value
                    # println((1,ng))
                     parI[ng] = (1-eps)/(eps* parOMin)
                    # objfunI_ngMinEps(parI_ng::Real) = (objfunI_ng(parI_ng::Real) - eps)
                    # tmpBnd = bigNumber #(degGroupsO[ng]/(NgroupsMemb[ng,2]*(N-1-NzerI)) )/(1-degGroupsO[ng]/(NgroupsMemb[ng,2]*(N-1-NzerI)) )
                    # brakets = [ 0,2* tmpBnd/parOMin]
                    # parI[ng] = tmpFindZero(objfunI_ngMinEps,brakets,parI[ng])
                elseif degGroupsI[ng]==0
                    #error("Remove zero col and rows")
                    parI[ng] = 0
                else
                    tmpBnd = (degGroupsI[ng]/(NgroupsMemb[ng,2]*(N-1-NzerO)) )/(1-degGroupsI[ng]/(NgroupsMemb[ng,2]*(N-1-NzerO)) )
                    #println((tmpBnd,degGroupsI[ng],NgroupsMemb[ng,2],NzerO))
                    brakets = [ 0,2* tmpBnd/parOMin]
                    parI[ng] = tmpFindZero(objfunI_ng,brakets,parI[ng])
                end
            end
            #parUn = log.([parI;parO])
            #err = firstOrderCond(Model;deg = deg,par = parUn , groupsInds = groupsInds)
            #println(err)
            parIMax = maximum(parI)
            parIMin = minimum(parI[parI.!=0])
            for nInd=1:NG
                ng = indsDegGroupsOrderO[nInd]
                #println((2,n))
                #put to zero each component of the gradient1
                function objfunO_ng(parO_ng::Real)
                    parVecOthGroups = parO_ng.*parI[1:end .!=ng]
                    parI_ng = parI[ng]
                    expGroupDegO = NgroupsMemb[ng,2]*( sum( NgroupsMemb[1:end .!=ng ,2].*( parVecOthGroups ./ (1 + parVecOthGroups)) )    + (NgroupsMemb[ng,2]-1)* (parI_ng*parO_ng)/(1+(parI_ng*parO_ng)) )
                    return  degGroupsO[ng] - expGroupDegO
                end
                #find zero of given component
                if degGroupsO[ng] == NgroupsMemb[ng,2]*(N-1-NzerI)
                    #if a group has the maximum possible out degree set par to a
                    #very high value
                    # println((2,ng))
                    parO[ng] = (1-eps)/(eps* parIMin)
                    # objfunO_ngMinEps(parO_ng::Real) = (objfunO_ng(parO_ng::Real) - eps)
                    # tmpBnd = bigNumber #(degGroupsO[ng]/(NgroupsMemb[ng,2]*(N-1-NzerI)) )/(1-degGroupsO[ng]/(NgroupsMemb[ng,2]*(N-1-NzerI)) )
                    # brakets = [ 0,2* tmpBnd/parIMin]
                    # parO[ng] = tmpFindZero(objfunO_ngMinEps,brakets,parO[ng])
                elseif degGroupsO[ng]==0
                    #error("Remove zero col and rows")
                    parO[ng] = 0
                else
                    tmpBnd = (degGroupsO[ng]/(NgroupsMemb[ng,2]*(N-1-NzerI)) )/(1-degGroupsO[ng]/(NgroupsMemb[ng,2]*(N-1-NzerI)) )
                    brakets = [ 0,2* tmpBnd/parIMin]
                    parO[ng] = tmpFindZero(objfunO_ng,brakets,parO[ng])
                end

            end
            parReIO = [parI;parO]
            parUnIO = -bigConstVal.*ones(parReIO)
            tmpInds = parReIO.>0
            #show(tmpInds)
            parUnIO[tmpInds] = log.(parReIO[tmpInds])
            err = firstOrderCond(Model;degIO = degIO,parGroupsIO = parUnIO , groupsInds = groupsInds)
        #    tmpMat = par.*par';expMat = putZeroDiag!(tmpMat ./ (1+tmpMat)) ;expDeg = sumSq(expMat,2)
            relErrVec = abs.(err ./ degGroups)
            #println(relErrVec)
            relErr =  maximum(relErrVec[degGroups.!=0])
            isnan(relErr) ? error() : ()
            # display(err')
             #println((mean(relErrVec[deg.!=0]),relErr))
            i+=1
        end
        #println([parI parO])
        #i<maxIt ? () : (println(i);error())
        parIO = [parI;parO]
        parUN = log.(parIO)
    end
    return parUN, i
 end
function estimate(model::NetModelDirBin1; degIO::Array{<:Real,1} = model.obs,
                    targetErr::Real=targetErrValStaticNets,identIter = false)
    ## Write and test a new way to estimate DirBOIn1 following Chatterjee, Diaconis and  Sly
    N,degI,degO = splitVec(degIO)
    ldegI = log.(degI);ldegO = log.(degO)
    (sum(degI) - sum(degO) )< 10^Base.eps() ?  L =sum(degO) : error()
    LperLink  = L/(N*N-N)
    unifstartval =  0.5*log((LperLink/((1-LperLink))))
    #inds = (ldegI.==0)
    #function that defines the iteration
    function phiFunc(vParUn::Array{<:Real,1})
        ~,vParReI,vParReO = splitVec(exp.(vParUn))
        matI = putZeroDiag(1 ./ (vParReI  .+ (1 ./ vParReO)' ))
        matO = putZeroDiag(1 ./ ( (1 ./ vParReI)  .+ vParReO'  ))
        outPhiI = (ldegI - log.(sumSq(matI,2)))
        outPhiO = (ldegO - log.(sumSq(matO,1)))
        #println(vParReI[ldegI.==0])
        #println(log.(sumSq(matI,2))[inds])
        #println(matI[inds,:][1,1:20])
        #println(matI[inds,:][2,1:20])
        return [outPhiI;outPhiO]
     end
    function checkFit(vParUn::Vector{<:Real})
        expMat = expMatrix(model,vParUn)
        nnzInds = degIO.!=0
        errsIO = [sumSq(expMat,2) - degI; sumSq(expMat,1) - degO][nnzInds]
        relErrsIO = abs.(errsIO) ./ (degIO[nnzInds])
        return prod(relErrsIO .< targetErr)
     end

    maxIt = 20000
    it=0
    it2 = 0
    vParUn = [ldegI;ldegO]# 0*unifstartval.*ones(2N)
    #starting from an identified set of par. and then identify at each step
    #amounts to doing an identified step. i.e. moving summing an identified
    #vector
    identIter ? vParUn =  identify(model,vParUn) : ()
    while (!checkFit(vParUn)) & (it<maxIt)
        vParUn = phiFunc(vParUn); it+=1
        identIter ? vParUn =  identify(model,vParUn) : ()
        #println(vParUn[degIO.==1][1:5])
    end
    #if did not obtain the required precision try with the slower method
    if it==maxIt
        outSlow = estimateSlow(model; degIO = degIO,targetErr = targetErr,
                            startVals =exp.(vParUn))

        vParUn = identify(model,outSlow[1])
        it2 = outSlow[2]
    end
    #zeroInds = degIO.==0
    #vParUn[zeroInds] = bigNegParFun(model,N)
    return vParUn,it,it2
 end
 #

function sampl(Model::NetModelDirBin1,Nsample::Int;  parGroupsIO::Array{<:Real,1}=Model.Par, groupsInds = Model.groupsInds )
    NGcheck,parI,parO = splitVec(parGroupsIO)
    parNodesIO = [parI[groupsInds];parO[groupsInds]]
    expMat = expMatrix(Model, parNodesIO)
    N = length(expMat[:,1])
    SampleMats =zeros(Int8,N,N,Nsample)
    for s=1:Nsample
        SampleMats[:,:,s] = samplSingMatCan(Model,expMat)
    end
    return SampleMats
 end



# Snapshot sequence Functions-----------------------------------------------------------------

struct SnapSeqNetDirBin1
    obsT:: Array{<:Real,2}  # In degrees and out Degrees for each t . The last
                            # dimension is always time. 2NxT
    Par::Array{<:Real,2} # In and Out Parameters [InW_T;outW_T] for each t
 end
SnapSeqNetDirBin1(obs::Array{<:Real,2}) = SnapSeqNetDirBin1(obs,zero(obs))
fooSnapSeqNetDirBin1 = SnapSeqNetDirBin1(zeros(2,3))

function estimate(Model::SnapSeqNetDirBin1; degsIO_T::Array{<:Real,2}=Model.obsT,
                    targetErr::Real=1e-4,identPost=false,identIter=false,
                    zeroDegFit="small")
   N2,T = size(degsIO_T);N = round(Int,N2/2)
   UnParT = zeros(2N,T)
   prog = 0
   for t = 1:T
       degs_t = degsIO_T[:,t]
      # println(degs_t)
       UnPar_t , Niter1,Niter2  = estimate(fooNetModelDirBin1; degIO = degs_t ,
                                    targetErr =  targetErr,identIter=identIter)
       zeroParI,zeroParO = zeroDegParFun(fooNetModelDirBin1,N;
                                        degsIO=degsIO_T[:,t],parIO = UnPar_t )
      # println((zeroParI,zeroParO))
       !isfinite(zeroParI) ? error() : ()
       !isfinite(zeroParO) ? error() : ()
       UnPar_t[findall(degs_t[1:N].==0)] .= zeroParI
       UnPar_t[N .+ findall(degs_t[1+N:2N].==0)] .= zeroParO
       identPost ? UnParT[:,t] = identify(fooNetModelDirBin1,UnPar_t) :
                                    UnParT[:,t] = UnPar_t
       #round(t/T,1)>prog ?  (prog = round(t/T,1);println(prog)) : ()

        #println((t,Niter1,Niter2))
   end
   sum(.!isfinite.(UnParT))>0 ? error() : ()

 return UnParT
 end

function sampl(Model::SnapSeqNetDirBin1, vParUnIO_T::Array{<:Real,2}; Nsample::Int=1  )
    N2,T = size(vParUnIO_T);N = round(Int,N2/2)
    SampleMats =zeros(Int8,Nsample,N,N,T)
    for s=1:Nsample
        for t = 1:T
            SampleMats[s,:,:,t] = samplSingMatCan(fooNetModelDirBin1,expMatrix(fooNetModelDirBin1,vParUnIO_T[:,t]))
        end
    end
    Nsample == 1  ?  (return squeeze(SampleMats,1)) : (return SampleMats )
  end

 #---------------------------------Binary DIRECTED Networks with predetermined COVARIATES

 #
 # struct  NetModelDirBin1X0 <: NetModelBin #Bin stands for (Binary) Adjacency matrix
 #                                        # X for regressors and 0(lower level pox) for link specific parameters
 #    "ERGM with in and out node specific parameters and dependence on link
 #    specific regressors trough link specific parameters (0)  "
 #    Amats:: BitArray{3}   #Different realizations of the adjqcency matrix Needed!
 #                        #Because the covariates parameters are link specific. Time
 #                        #is the last dimension
 #    regressors0:: Array{<:Real,3}
 #    parDegs:: Array{<:Real,2} # One parameter per node  UNRESTRICTED VERSION  [InPar OutPar]
 #    parRegs:: Array{<:Real,2} #one parameter for each link associated with each regressor
 # end
 # fooNetModelDirBin1X0 =  NetModelDirBin1X0(tmpArrayAdjmatrEqual(6,10),ones(6,6,10),zeros(6,2),ones(6,6))
 # NetModelDirBin1X0(Amats:: BitArray{3},regs0::Array{<:Real,3}) =  (N =length(Amats[:,1,1]); NetModelDirBin1X0(Amats,regs0, zeros(N ,2) , zeros(N,N) ))

function defineConstDegs(degsIO_T;thVarConst = 0.005)
        N2 = length(degsIO_T[:,1])
        isConst =  [StatsBase.trimvar(degsIO_T[i,:];prop = 0.001)  for i =1:N2].< thVarConst

        ~,isConI,isConO = splitVec(isConst)
        @show sum(isConst)
        return sum(isConst),isConI,isConO
    end
#Some useful wrap functions
function estManyAR(obs::Array{<:Real,2};p::Int=1,isCon::BitArray{1} = falses(length(obs[:,1])))
    #obs is a matrix with a time series for each row (time is last dimension)
    D,T = size(obs)
    parEst = zeros(D,p+1)
    arp = AReg.ARp(ones(p+1),.1,[])

    for d=1:D
        if !isCon[d]
            parEst[d,:] = AReg.fitARp(obs[d,:],arp).estimates
        else
            #if the ts is costant, use only the unconditional mean with zero persistence
            parEst[d,1] = mean(obs[d,:])
        end
    end
    return parEst
 end
function oneStepForAR1(obs_t::Float64,AR1Par::Vector{<:Real};steps::Int = 1)
    #given w and β AR1 parameters, compute the one step ahead forecast from obs_t
    if steps==1
        obstpSteps = AR1Par[1] + AR1Par[2] * obs_t
    else
        obstpSteps = AR1Par[1] * sum(AR1Par[2].^Vector(0:steps-1)) + AR1Par[2].^steps * obs_t
    end
    return obstpSteps
 end
#Estimate AR1 on the training sample
function forecastEvalAR1Net(allFit::Array{Float64,2}, obsNet_T::BitArray{3},Ttrain::Int;mat2rem=falses(obsNet_T[:,:,1]) )
    N2,T = size(allFit)
    Ttest = T-Ttrain
    trainFit = allFit[:,1:Ttrain]
    degsIO_T = [sumSq(obsNet_T,2);sumSq(obsNet_T,1)]
    Nconst,isConIn,isConOut = defineConstDegs(degsIO_T[:,1:Ttrain],thVarConst =0.005 )
    estAR1= estManyAR(trainFit;isCon = [isConIn;isConOut])
    ## One step ahead forecasts for each fitness for each time in Test sample
    foreFit = zeros(allFit)
    for n=1:N2 ,t=1:T-1
        #compute forecasting for all links, diagonal will be disregarded after
        foreFit[n,t+1] = oneStepForAR1(allFit[n,t],estAR1[n,:])
    end
    tp,fp,AUC,realVals,foreVals =    nowCastEvalFitNet(foreFit[:,Ttrain+1:end],obsNet_T[:,:,Ttrain+1:end];mat2rem=mat2rem)
    return tp,fp,AUC,realVals,foreVals
 end
function nowCastEvalFitNet(allFit::Array{Float64,2}, obsNet_T::BitArray{3};
            expMat_T :: Array{Float64,3}=zeros(10,10,50),mat2rem=falses(obsNet_T[:,:,1]), shift::Int=0,plotFlag = true )



    sum(expMat_T) == 0  ? inputExpMat =false : inputExpMat=true
    N2,T = size(allFit)
    inputExpMat  ? T = size(expMat_T)[3] : ()
    println((N2,T))
    #println(size(expMat_T))
    obsNet_Test_T= obsNet_T

    # se l'expected matrices sono date come input nel seguto non calcolarle
    #ricorda che se expMat_T vengono da multi steps ahead forecasting,non va usato
    # il temporal shift
    inputExpMat ? (shift!=0 ? error : ()) : ()
    ##use input fitnesses
    testFit = allFit
    #disregard diagonal elements and those that are in mat2rem
    noDiagIndnnc = putZeroDiag((.!mat2rem)    )
    @show Nlinksnnc =sum(noDiagIndnnc)
    if shift>=0
            println(Nlinksnnc*(T-shift))
            println(T)
        foreVals = 100 * ones(Nlinksnnc*(T-shift))
        realVals = falses( Nlinksnnc*(T-shift))
        lastInd=1
        for t=1:T-shift
            adjMat =  obsNet_Test_T[:,:,t+shift]
            if inputExpMat
                expMat =  expMat_T[:,:,t]
                #println(size(expMat))
            else
                expMat = expMatrix(fooNetModelDirBin1,testFit[:,t])
            end
            foreVals[lastInd:lastInd+Nlinksnnc-1] = expMat[noDiagIndnnc]
            realVals[lastInd:lastInd+Nlinksnnc-1] = adjMat[noDiagIndnnc]

            lastInd += Nlinksnnc
        end
    else
        foreVals = 100 * ones(Nlinksnnc*(T+shift))
        realVals = falses( Nlinksnnc*(T+shift))
        lastInd=1
        for t=1:T+shift
            adjMat =  obsNet_Test_T[:,:,t]
            if inputExpMat
                expMat =  expMat_T[:,:,t]
            else
                expMat = expMatrix(fooNetModelDirBin1,testFit[:,t-shift])
            end


            foreVals[lastInd:lastInd+Nlinksnnc-1] = expMat[noDiagIndnnc]
            realVals[lastInd:lastInd+Nlinksnnc-1] = adjMat[noDiagIndnnc]

            lastInd += Nlinksnnc
        end
    end
 # println(realVals)
    tp,fp,AUC = rocCurve(realVals,foreVals;plotFlag = plotFlag)
    return tp,fp,AUC,realVals,foreVals
 end

function forecastEvalAR1Net(      obsMat_allT::BitArray{3},Ttrain::Int;identPost = false,identIter=true,mat2rem = falses(obsMat_allT[:,:,1]) )
        #if only the matrices obs are given then estimate the one fit model for each t
        degsIO_T = [sumSq(obsMat_allT,2);sumSq(obsMat_allT,1)]
        snapModData = SnapSeqNetDirBin1(degsIO_T)
        allFit =  estimate(snapModData; identPost = identPost,identIter=identIter)
        return   forecastEvalAR1Net(allFit, obsMat_allT,Ttrain;mat2rem = mat2rem)

end

function dgpFitVarN(N::Int,degMin::Int,degMax::Int;exponent = 1)

    # exponent determines the convexity of the distributions of degrees between
    # deg min and deg max. increasing exponent --> decreasing density
    # if we let deg max grow with N than
    #generate a degree distribution
    tmp = (range(degMin,stop=degMax,length=N)).^exponent

    degIOUncMeans = repeat(round.(Int,Vector(tmp.*(degMax/tmp[end]) )),2)
    degIOUncMeans[degIOUncMeans.<degMin] = degMin
    meanFit ,~ = StaticNets.estimate(StaticNets.NetModelDirBin1(degIOUncMeans),degIO = degIOUncMeans)
    expDegs = StaticNets.expValStats(fooNetModelDirBin1,meanFit)
    return meanFit, expDegs
end

function dgpDynamic(model::NetModelDirBin1,mode,N::Int,T::Int;NTV = 2,degIOUncMeans::Array{Float64,1} =zeros(10),
                    degb = [10, 40])
    if sum(degIOUncMeans)==0 # if the unconditional means of degrees are not Given then fix them
        degIOUncMeans = range(1,stop=10,length=N)
        error()
    else
        N = round(Int,length(degIOUncMeans)/2)

    end
    #data unconditional mean definisci ampiezza della variazione, che sarà ampiezza seno e varianza AR1


    ## generate the dynamical fitnesses

    #println((N,NTV))
    indsTV = round.(Int,range(2, stop=N, length=NTV))
    indLogTV = falses(N);indLogTV[indsTV] .= true; indLogTV = repeat(indLogTV,2)

    meanFit ,~ = StaticNets.estimate(StaticNets.NetModelDirBin1(degIOUncMeans),degIO = degIOUncMeans)
    #define the oscillations bounds for the fitnesses
    periodOsc = round(Int,T/2)

    # all fits start at the same point with different fases and amplitudes
    #define a variation
    ~,meanfitI,meanfitO = splitVec(meanFit)
    #shuffle ranUtilitiesly the fases of the sines but fix the ranUtilities seed for reproducibility

    #solve for x :deglb =  (N-1) * 1/(1+exp( - x - meanfitO))
    fitbI = - mean(meanfitO) .- log.( (N-1) ./ degb .-1)
    fitbO = - mean(meanfitI) .- log.( (N-1) ./ degb .-1)
    ampI = ones(NTV).*Vector(range(0.4,stop=1, length=NTV).*(diff(fitbI)))#(diff(fitbI))#
    ampO = ones(NTV).* Vector(range(0.4,stop=1,length=NTV).*(diff(fitbO))) # (diff(fitbO)) #
    #center the dynamical parameters at the middle of the bounds
    meanFit[indsTV] =  ampI ./ 2  .+ fitbI[1]
    meanFit[N.+indsTV] = ampO ./ 2  .+ fitbO[1]

    times = Array{Float64,1}(1:T)
    dynFit = repeat(meanFit,1,T)

    if uppercase(mode) ==  uppercase("sin")
        println(NTV)
        phases =  rand( MersenneTwister(1234),2NTV)*((2pi))
        for n=1:NTV
            dynFit[indsTV[n],:]   = meanFit[indsTV[n]]  .+ ampI[n]*sin.( times.*((2pi/periodOsc)) .+ phases[n] )
            dynFit[N+indsTV[n],:] = meanFit[N+indsTV[n]]  .+ ampO[n]*sin.( times.*((2pi/periodOsc)) .+ phases[n+NTV] )
        end
    elseif uppercase(mode) ==  uppercase("const")
        for n=1:NTV
            dynFit[indsTV[n],:]   = meanFit[indsTV[n]] .*ones(T)
            dynFit[N+indsTV[n],:] = meanFit[N+indsTV[n]]  .*ones(T)
        end
    elseif uppercase(mode) ==  uppercase("steps")
        Nsteps = 2
        phases =  rand( MersenneTwister(1234),1:T,2NTV)

        for n=1:NTV
            minPar = meanFit[indsTV[n] ] .- ampI[n]/2
            maxPar = meanFit[indsTV[n]]  .+ ampI[n]/2
            #steps between min and max values
            dynFit[indsTV[n],:]  =  randSteps(minPar,maxPar,Nsteps,T)
            maxPar = meanFit[indsTV[n]+N] + ampO[n]/2
            #steps between min and max values
            dynFit[indsTV[n]+N,:]  = randSteps(minPar,maxPar,Nsteps,T)
            #circyularly shift time index by the ranUtilities phase of each
            if false
                #circyularly shift time index by the ranUtilities phase of each
                dynFit[indsTV[n],:]  = circshift(dynFit[indsTV[n],:]  ,phases[n])
                dynFit[indsTV[n]+N,:]  = circshift(dynFit[indsTV[n]+N,:],phases[n])
            end
            minPar = meanFit[indsTV[n]+ N] - ampO[n]/2
        end


    elseif uppercase(mode) == uppercase("ar1")

        sigma = 0.01
         B = 0.999 #rand([0.7,0.99])

        for n=1:NTV
            scalType = "uniform"
            minPar = meanFit[indsTV[n] ] - ampI[n]/2
            maxPar = meanFit[indsTV[n]] + ampI[n]/2

            dynFit[indsTV[n],:]  = dgpAR(meanFit[indsTV[n]],B,sigma,T,minMax=[minPar,maxPar];scaling = scalType)
            #println((meanFit[indsTV[n]]  ,[minPar,maxPar]) )
            minPar = meanFit[indsTV[n]+ N] - ampO[n]/2
            maxPar = meanFit[indsTV[n]+N] + ampO[n]/2

            dynFit[indsTV[n]+N,:]   = dgpAR(meanFit[indsTV[n]+N],B,sigma,T,minMax=[minPar,maxPar];scaling = scalType)


        end

    end


    #

    for t=1:T dynFit[:,t] = identify(model,dynFit[:,t]) end
    return dynFit, indLogTV
end

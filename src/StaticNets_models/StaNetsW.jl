#__precompile__()
module StaticNetsW

using Utilities,Distributions, StatsBase,Optim, LineSearches, StatsFuns,Roots

export splitVec, distributeAinB, distributeAinB!, distributeAinVecN,maxLargeVal

## STATIC NETWORK MODEL
abstract type NetModel end
abstract type NetModelW <: NetModel end
abstract type NetModelBinW <: NetModel end
abstract type NetModelWcount <: NetModelW end
abstract type NetModelBin <: NetModel end

#constants
targetErrValStaticNets = 0.005
targetErrValStaticNetsW = 1e-5
bigConstVal = 10^6
maxLargeVal =  1e40# 1e-10 *sqrt(prevfloat(Inf))
minSmallVal = 1e2*eps()


#-------------------------Count valued Weighted Directed Networks---------------------------------------


struct  NetModelDirWCount1 <: NetModelW #Bin stands for (Binary) Adjacency matrix
    "ERGM for directed networks with geometrically(NB??) distributed count weights,
    in and out strenghts as statistics and possibility to have groups of nodes
    associated with a single pair of in out parameters. "
    obs:: Array{<:Int,1} #[InStr;OutStr] degrees of each node (even if groups are present)
    Par::Array{Array{<:Real,1},1} # One parameter per group of nodes  UNRESTRICTED VERSION  [[InPar;OutPar], scalePar ]
    groupsInds::Array{<:Real,1} #The index of the groups each node belongs to
    # two parameters per node (In and out)
 end
fooNetModelDirWCount1 =  NetModelDirWCount1(ones(Int,6),[zeros(6),zeros(1)],ones(3))
NetModelDirWCount1(StrIO:: Array{<:Int,1}) =  NetModelDirWCount1(StrIO,zeros(length(StrIO[:])),[ones(round(Int,length(StrIO[:])/2)),ones(1)])

function DegSeq2graphDegSeq(Model::NetModelDirWCount1, StrIO::Array{<:Real,1}; Amat::BitArray{2} = trues(round((Int,length(StrIO)/2)),round(Int,(length(StrIO)/2)) ))
    #From a degree sequence find another one that is certanly graphicable
    #The matrix associated with the starting degree sequence is very hubs oriented

    # One stub for each unit weight
    N, StrI,StrO = splitVec(StrIO)
    sum(StrI) == sum(StrO)?(): (println(sum(StrI) - sum(StrO));error())
    inds = Vector(1:N)
    StubsI_inds = Int.(sortrows([StrI inds],rev=true))#decreasing sort
    StubsO_inds = Int.(sortrows([StrO inds],rev=true))

    mat = zeros(Int64,N,N)
    for i in 1:N
        j=1
        while j < N
            ind_i = StubsI_inds[i,2]
            ind_j = StubsO_inds[j,2]
            while (StubsI_inds[i,1]>0)& (StubsO_inds[j,1]>0)&(ind_i!=ind_j)&Amat[i,j]
                mat[ind_i,ind_j] += 1
                StubsI_inds[i,1]  -= 1
                StubsO_inds[j,1]  -= 1
                #println(i,j, StubsI_inds[i,1],StubsO_inds[j,1]  )
            end
            j += 1
        end
    end
    #display(StubsI_inds[:,1])
    ((sum(mat[Array{Bool}(eye(mat))])) != 0) ?error():()
    graphSeq= [squeeze(sum(mat,2),2); squeeze(sum(mat,1),1)]
    return graphSeq
 end


#-------------------------Countinuosly Weighted   Directed Networks Fixed A---------------------------------------


struct  NetModelDirW1Afixed <: NetModelW #Bin stands for (Binary) Adjacency matrix
    "ERGM for directed networks with Bin A fixed and Gamma distributed
     Continuous  Weights and finite probability of observing zeros.
     Thes statistics are in and out strenghts, and a binary adjacency matrix.  "
    obs:: Array{<:Real,1}  # [InStr;OutStr]
    Amat:: BitArray{2} #Binary Adjacency matrix Amat
    Par::Array{Array{<:Real,1},1} # One parameter per group of nodes  UNRESTRICTED VERSION  [ [InPar;OutPar], ScalePar]
 end
fooNetModelDirW1Afixed =  NetModelDirW1Afixed(ones(6),trues(3,3),[zeros(6),zeros(1)])
NetModelDirW1Afixed(StrIO::Array{<:Real,1},Amat::BitArray{2}) =  NetModelDirW1Afixed(StrIO,Amat,[zeros(Float64,length(StrIO[:])), zeros(1)]  )

function DegSeq2graphDegSeq(Model::NetModelDirW1Afixed,StrIO::Array{<:Real,1};Amat::BitArray{2}=Model.Amat)
    #From a strengths sequence find another one that is certanly graphicable

    # pe deg2graph deg posso importare le funzioni in c++ scritte da veraart oppure fare
    # una discretizzazione (ad esempio di modo che la strength piu' piccola sia divisa in 10N
    # unità) e poi passarla in pasto alla versione Discreta pesata
   N,StrI,StrO = splitVec(StrIO)
    sum(StrI) == sum(StrO)?():error()

    minContinuousVal = minimum(StrIO[StrIO.!=0])
    # choose the minimum unit weight such that the smallest node has  10N stubs
    unitW = minContinuousVal/(10N)
    #display(minContinuousVal)
    StrCountIO = ceil(Int,StrIO./(unitW))

    graphSeq = unitW .* DegSeq2graphDegSeq(fooNetModelDirWCount1,StrCountIO;Amat = Amat)
    return graphSeq
 end

bndPar2uBndPar(Model::NetModelDirW1Afixed, bndParNodesIO::Array{<:Real,1} ;indMin::Int = -1 ) =
    bndPar2uBndPar(Model, bndParNodesIO ;indMin = indMin)

uBndPar2bndPar(Model::NetModelDirW1Afixed, uBndParNodesIO::Array{<:Real,1} ;indMin::Int = -1 ) =
    uBndPar2bndPar(Model, uBndParNodesIO ;indMin = indMin)

samplSingMatCan(Model::NetModelDirW1Afixed,expMat::Array{<:Real,2},Amat::BitArray{2};α=1.0::Real) =
    (expMat[.!Amat] = 0; samplSingMatCan(Model::NetModelDirW1,expMat;α=α) )

function expMatrix2(Model::NetModelDirW1Afixed, parNodesIO::Array{<:Real,1},Amat::BitArray{2} ;logalpha = 0.0 )
    "Given the vector of model parameters (groups), return the product matrix(often useful
    in likelihood computation) and the expected matrix"
    α = exp( logalpha )
    N,parI,parO = splitVec(parNodesIO)
    parMat = exp.(parI) .+ exp.(parO)'
    expMat = putZeroDiag(α./( parMat) )
    expMat[.!Amat] = 0
    infInd = isfinite.(expMat)
    all(infInd) ?() :( display([squeeze(prod(infInd,2),2) squeeze(prod(infInd,1),1)]); error())
    return parMat,expMat
 end
expMatrix(Model::NetModelDirW1Afixed, parNodesIO::Array{<:Real,1} ,Amat::BitArray{2};logalpha = 0.0  ) =
        expMatrix2(Model, parNodesIO,Amat;logalpha = logalpha )[2]# return only the expected matrix given the parameters
expValStatsFromMat(Model::NetModelDirW1Afixed,expMat::Array{<:Real,2}  ) =
                        expValStatsFromMat(fooNetModelDirBin1,expMat)

expValStats(Model::NetModelDirW1Afixed, parNodesIO::Array{<:Real,1},Amat::BitArray{2} ) = expValStatsFromMat(Model, expMatrix(Model,parNodesIO ,Amat)) # if only the paramters are given as input
function firstOrderCond(Model::NetModelDirW1Afixed;strIO::Array{<:Real,1} = Model.obs,parNodesIO::Array{<:Real,1} = Model.Par, Amat::BitArray{2}=Model.Amat )
    "Given the model, the strengths sequence, the parameters and Adjacency matrx ,
    return the First order conditions. Gradient of the loglikelihood, or system
    of differences between expected and observed values of the statistics (either
    degrees(In and Out) of single nodes or of a group). "
    foc = strIO .-expValStats(Model,parNodesIO,Amat)
    return foc
 end
function identify!(Model::NetModelDirW1Afixed,parIO::Array{<:Real,1})
    "Given a vector of parameters, return the transformed vector that verifies
    an identification condition. Do not do anything if the model is identified."
    #set the first of the in parameters equal to one (Restricted version)
    # N,parI,parO = splitVec(parIO)
    # parIO = log.([ exp.(parI) - exp.(parI[1]) ;  exp.(parO) + exp.(parI[1]) ])
    return parIO
  end

function estimateSlow(Model::NetModelDirW1Afixed; strIO::Array{<:Real,1} = Model.obs,
                Amat::BitArray{2}=falses(10,10) , targetErr  =targetErrValStaticNetsW,bigConst   = bigConstVal,identIter = false,startVals::Array{Float64,1} =zeros(10) )
    "Given model type, observations of the statistics and groups assignments,
    estimate the parameters."
    N,strI,strO = splitVec(strIO)
    (sum(strI) - sum(strO) )< 10^Base.eps()? S =sum(strO) :error()
    SperLink  = S/(N*N-N)
    # The parameters inside this function are Restricted to be in [0,1] hence
    # they need to be converted (log(par)) in order to be used outside

    # if no Amat is given assume the one implied by the strengths sequence
    sum(Amat ) == 0? Amat = (strI.>0).*(strO'.>0):()
    #uniform parameter as starting value
    unifstartval = 0.5/SperLink;
    bigConstVal = 1e4
    if sum(startVals) == 0
        local parI = unifstartval*ones(Float64,N)
        local parO = unifstartval*ones(Float64,N)
    else
        N,parI,parO = splitVec(startVals)
    end

    #bigNumber = 1e3./targetErr
    eps =  targetErr*1e-12 #10*Base.eps()
    i=0
    maxIt = 500
    relErr = 1
    tol_x = 10*Base.eps()

    function tmpFindZero(objFun::Function,brak::Array{<:Real,1},start::Real)
        #    wrapper find_zero options
        println((brak[1],brak[2]))
        println((objFun(brak[1]),objFun(brak[2])))
        #println()
        out = find_zero(objFun,(brak[1],brak[2]),Bisection())
        println((out,objFun(out)))
        return out
    end
    degsI = sumSq(Amat,2)
    degsO = sumSq(Amat,1)
    zerIndI = (strI .== 0 )
    (zerIndI == (degsI.==0) )? () : error()
    zerIndO = (strO .== 0 )
    @show (zerIndO)
    @show  (degsO.==0)
    (zerIndO == (degsO.==0) )? () : error()

    while (i<maxIt)&(relErr>targetErr)
        for n=1:N
            indGzO= (parO.>0).& (.!zerIndO); indGzO[n] = false
            indLzO= (parO.<0).& (.!zerIndO); indLzO[n] = false
            lBndParI = maximum( [(sum(indGzO)>0? -minimum(parO[indGzO]): [] ); (sum(indLzO)>0 ? -minimum(parO[indLzO]) : []) ] )

            println((1,n))
            #put to zero each component of the gradient1
            function objfunI_n(parI_n::Real)
                parVecOthGroups = parI_n .+ parO[Amat[:,n]]
                expStrI =  sum(  1./( parVecOthGroups) )   ##parO_n = parO[n]   + (NgroupsMemb[ng,2]-1)* ( 1/(parI_ng + parO_ng)) )
                return  strI[n] - expStrI

            end
            #find zero of given component

            if strI[n] <= 1e-3
                #error("Remove zero col and rows")
                parI[n] = bigConstVal
                #println(0)
            else
                tmpBnd = 1/(strI[n]/(N*(N-1)) )
                #println((lBndParI))
                brakets =[ lBndParI + 1/S ;  bigConstVal]
                #ERRORE PROBABILMENTE IN SCAMIO TRA MATRICE A E A'
                #DEVO ANCHE CORREGGERE IL LOWER BOUND PER TENERE CONTO DGLI
                #ELEMENTI NULLI DELLA AMAT
                parI[n] = tmpFindZero(objfunI_n,brakets,parI[n])
            end
        end
        parUnIO = ([parI;parO])
        identIter? parUnIO = identify(fooNetModelDirW1Afixed,parUnIO):()
        #println(parI)
        err = firstOrderCond(Model; strIO = strIO,parGroupsIO = parUnIO , groupsInds = groupsInds)
        #display(err)


        for n=1:N

            indGzI=(parI.>0).& (.!zerIndI); indGzI[n] = false
            indLzI=(parI.<0).& (.!zerIndI); indLzI[n] = false
            lBndParO = maximum( [(sum(indGzI)>0? -minimum(parI[indGzI]): [] ); (sum(indLzI)>0 ? -minimum(parI[indLzI]) : []) ] )

            #println((2,nInd))
            #put to zero each component of the gradient1
            function objfunO_n(parO_n::Real)
                parVecOthGroups = parO_n .+ parI[1:end ][Amat[n,1:end]]
                expStrO = sum(1./(parVecOthGroups))#NgroupsMemb[ng,2]* sum( NgroupsMemb[1:end .!=ng ,2].*(1./( parVecOthGroups)) )    + (NgroupsMemb[ng,2]-1)* (1/(parI_ng + parO_ng)) )
                return  strO[n] - expStrO
            end
            #find zero of given component
            if strO[n]<= 1e-3
                #error("Remove zero col and rows")
                parO[n] = bigConstVal
                #println(bigConstVal)
                #println(0)
            else
                tmpBnd = 1/(strGroupsO[n]/(NgroupsMemb[n,2]*(N-1)) )
                #println((lBndParO))
                brakets =[lBndParO + 1/S ;  bigConstVal ]
                parO[n] = tmpFindZero(objfunO_n,brakets,parO[n])
            end

        end
        parUnIO = ([parI;parO])
        #println(parO)
        err = firstOrderCond(Model; strIO = strIO,parGroupsIO = parUnIO , groupsInds = groupsInds)
        #display(err)
    #    tmpMat = par.*par';expMat = putZeroDiag!(tmpMat./(1+tmpMat)) ;expDeg = squeeze(sum(expMat,2),2)
        relErrVec = abs.(err./strGroups)[strGroups .>= 1e-3]
        #println(relErrVec)
        #take into account also groups with zero strength. multiply by target error
        # implies that the expected strength for groups with zero strength
        # has to be lower than 1
        absErrRescVec = abs.(err[strGroups.==0]) .* targetErr
        relErr = maximum([relErrVec;absErrRescVec])
        isnan(relErr)?error():()
        # display(err')
         #println((mean(relErrVec[deg.!=0]),relErr))
        i+=1
        #i==1?error():()
    end
    #println([parI parO])
    i<maxIt?():(println(i);error())
    # the parameter of the first index has to be unique because I want it to be the only infinite
    strIO[1]==0 ? parI[1] = bigConstVal*50:()
    parIO = [parI;parO]
    #println(unique(parIO))
    α = exp(Model.Par[2])
    parIOout = parIO #.* α

    outPar = identify!(Model,parIO)
    outMod = NetModelDirW1(strIO,[outPar,Model.Par[2]],groupsInds)

    return outPar, i , outMod
 end
function sampl(Model::NetModelDirW1Afixed,Nsample::Int;  parNodesIO::Array{Array{<:Real,1}}=Model.Par,Amat::BitArray{2}=Model.Amat  )
    parNodesIO = [parI[groupsInds];parO[groupsInds]]
    expMat = expMatrix(Model, parNodesIO,Amat)
    N = length(expMat[:,1])
    SampleMats =zeros(Float64,N,N,Nsample)
    logα = parArray[2][1]
    α = exp(logα)
    for s=1:Nsample
        SampleMats[:,:,s] = samplSingMatCan(Model,expMat,Amat,α)
    end
    return SampleMats
 end
function estimateIPFMat(strIO::Array{<:Real,1},A::Array{<:Real,2} ; targetErr  =100 * targetErrValStaticNetsW)
    "Given In and Out strength sequences and a binary matrix A, return the matrix obtained with IPF "
    N,sI,sO = splitVec(strIO)
    nnzInds = strIO.!=0
    ~,nnzIndsI,nnzIndsO = splitVec(nnzInds)
    W = Float64.(A)
    it = 0;maxAbsRelErr = 10
    while maxAbsRelErr > targetErr
        sIprev = sum(W,2)

        W[nnzIndsI,:] = W[nnzIndsI,:].*(sI[nnzIndsI]./sIprev[nnzIndsI])
        sOprev = sum(W,1)
        W[:,nnzIndsO] = (sO[nnzIndsO]./sOprev[nnzIndsO])'.*W[:,nnzIndsO]
        strIO_IPF = expValStatsFromMat(Model,W)
        err = strIO_IPF - strIO
        relErr = err[nnzInds]./strIO[nnzInds]
        maxAbsRelErr = maximum(abs.(relErr))
        it +=1
    end
    return W,it
 end

 #-------------------------Countinuosly Weighted Directed Networks---------------------------------------


 struct  NetModelDirW1 <: NetModelW #Bin stands for (Binary) Adjacency matrix
     "ERGM for directed networks with Gamma distributed Continuous Weights,
     in and out strenghts as statistics and possibility to have groups of nodes
     associated with a single pair of in out parameters. "
     obs:: Array{<:Real,1} #[InStr;OutStr] degrees of each node (even if groups are present)
     Par::Array{Array{<:Real,1},1} # One parameter per group of nodes  UNRESTRICTED VERSION  [ [InPar;OutPar], ScalePar]
     groupsInds::Array{Int,1} #The index of the groups each node belongs to
     # two parameters per node (In and out)
  end
 fooNetModelDirW1 =  NetModelDirW1(ones(6),[zeros(6),zeros(1)],ones(3))
 NetModelDirW1(StrIO:: Array{<:Real,1}) =  NetModelDirW1(StrIO,[zeros(Float64,length(StrIO[:])), zeros(1)] ,Int.(1:round(Int,length(StrIO[:])/2))  )

 function DegSeq2graphDegSeq(Model::NetModelDirW1,StrIO::Array{<:Real,1})
     #From a degree sequence find another one that is certanly graphicable

     # pe deg2graph deg posso importare le funzioni in c++ scritte da veraart oppure fare
     # una discretizzazione (ad esempio di modo che la strength piu' piccola sia divisa in 10N
     # unità) e poi passarla in pasto alla versione Discreta pesata
     N,StrI,StrO = splitVec(StrIO)
     sum(StrI) == sum(StrO)?():error()

     minContinuousVal = minimum(StrIO[StrIO.!=0])
     # choose the minimum unit weight such that the smallest node has  10N stubs
     unitW = minContinuousVal/(10N)
     #display(minContinuousVal)
     StrCountIO = ceil(Int,StrIO./(unitW))
     graphSeq = unitW .* DegSeq2graphDegSeq(fooNetModelDirWCount1,StrCountIO)
     return graphSeq
  end

 function bndPar2uBndPar(Model::NetModelDirW1, bndParNodesIO::Array{<:Real,1} ;indMin::Int = -1 )
     N,bndParI, bndParO = splitVec(bndParNodesIO)
     bndParI[1] == 0 ? () :error()
     uBndParO = log.(bndParO)
     indMin == -1 ? indMin = indmin(bndParO):()
     #println(indMin)
     # Define the unbounded version of the In par s.t. it can really takeall allowed values
     # without violating λ_i +λ_j > 0
     uBndParI = log.(bndParO[indMin] .+ bndParI )
     uBndParI[1] = -Inf
     uBndParIO = [uBndParI;uBndParO]
     return  uBndParIO
  end
 function uBndPar2bndPar(Model::NetModelDirW1, uBndParNodesIO::Array{<:Real,1} ;indMin::Int = -1 )
     N,uBndParI, uBndParO = splitVec(uBndParNodesIO)
     uBndParI[1] == -Inf ? () :error()
     expLambdaO = exp.(uBndParO)
     indMin == -1 ? indMin = indmin(expLambdaO):()
     # Define the unbounded version of the In par s.t. it can really takeall allowed values
     # without violating λ_i +λ_j > 0
     bndParI = - expLambdaO[indMin] .+ exp.(uBndParI)
     #identification
     bndParI[1] = 0
     bndParO =   expLambdaO + minSmallVal

     #bndParO[indzeroO] = minSmallVal
     bndParIO = [bndParI;bndParO]
     indInf = bndParIO.>= maxLargeVal

     bndParIO[indInf] = maxLargeVal
     return  bndParIO
  end
 function linSpacedPar(Model::NetModelDirW1,Nnodes::Int;Ngroups = Nnodes,deltaN::Int=3,graphConstr = true)
     error()
     "for a given static model return a widely spaced sequence of groups parameters"
     N=Nnodes
     #generate syntetic degree distributions, equally spaced, and estimate the model
     tmpdegsIO = round.([ Vector(linspace(deltaN,N-1-deltaN,N));Vector(linspace(deltaN,N-1-deltaN,N))])
     if graphConstr
         tmpdegsIO  =   DegSeq2graphDegSeq(fooNetModelDirBin1,tmpdegsIO)
     end

     if Ngroups>N#if more groups than nodes are required throw an error
         error()
     else
         groupsInds = distributeAinVecN(Vector(1:Ngroups),Nnodes)
     end
     out,it,estMod = estimate(NetModelDirW1(tmpdegsIO);groupsInds = groupsInds)
     return out,tmpdegsIO
  end
 function samplSingMatCan(Model::NetModelDirW1,expMat::Array{<:Real,2};α=1.0::Real)
     """
     given the expected matrix sample one ranUtilities matrix from the corresponding pdf
     """
     N1,N2 = size(expMat)
     N1!=N2?error():N=N1
     out = zeros(Float64,N,N)
     for c = 1:N
         for r=1:N
             if r!=c
                 #The Gamma in Distributions takes as input α and θ = expVal/α
                 #println(expMat[r,c]/α)
                 expMat[r,c]/α > 0 ? () :  (println(expMat[r,c]);  error() )
                 out[r,c] = rand(Gamma(α,expMat[r,c]/α))
             end
         end
     end
     return out
  end
 function expMatrix2(Model::NetModelDirW1, parNodesIO::Array{<:Real,1} ;logalpha = 0.0 )
     "Given the vector of model parameters (groups), return the product matrix(often useful
     in likelihood computation) and the expected matrix"
     α = exp( logalpha )
     N,parI,parO = splitVec(parNodesIO)
     parMat = parI .+ parO'
     expMat = putZeroDiag(α./( parMat) )
     infInd = isfinite.(expMat)
     all(infInd) ?() :( display([squeeze(prod(infInd,2),2) squeeze(prod(infInd,1),1)]); error())
     return parMat,expMat
  end
 expMatrix(Model::NetModelDirW1, parNodesIO::Array{<:Real,1} ;logalpha = 0.0  ) =
         expMatrix2(Model, parNodesIO;logalpha = logalpha )[2]# return only the expected matrix given the parameters
 expValStatsFromMat(Model::NetModelDirW1,expMat::Array{<:Real,2}  ) =
                         expValStatsFromMat(fooNetModelDirBin1,expMat)

 expValStats(Model::NetModelDirW1, parNodesIO::Array{<:Real,1} ) = expValStatsFromMat(Model, expMatrix(Model,parNodesIO )) # if only the paramters are given as input
 function firstOrderCond(Model::NetModelDirW1;strIO::Array{<:Real,1} = Model.obs,
                     strGroups::Array{<:Real,1} = zeros(4), parGroupsIO::Array{<:Real,1}=Model.Par, groupsInds::Array{<:Real,1} = Model.groupsInds )
     "Given the model, the degree sequence, the parameters and groups assignements,
     return the First order conditions. Gradient of the loglikelihood, or system
     of differences between expected and observed values of the statistics (either
     degrees(In and Out) of single nodes or of a group). "
     N,strI,strO = splitVec(strIO)
     groupsNames = unique(groupsInds)
     NG = length(groupsNames)
     NGcheck,parI,parO = splitVec(parGroupsIO)
     NG==NGcheck? ():error()
     if sum(strGroups)==0
         strGroups = zeros(2NG) ; for i =1:N strGroups[groupsInds[i]] += strI[i];strGroups[NG + groupsInds[i]] += strO[i]; end
     end
     parNodesIO = [parI[groupsInds];parO[groupsInds]]
     N,expStrI,expStrO = splitVec(expValStats(Model,parNodesIO))
     expGroupsStr = zeros(2NG) ; for i =1:N expGroupsStr[groupsInds[i]] += expStrI[i];expGroupsStr[NG + groupsInds[i]] += expStrO[i]; end
     # println(groupsInds)
     # println(size(strGroups))
     # println(size(expGroupsStr))
     foc = strGroups .- expGroupsStr
     return foc
  end
 function identify!(Model::NetModelDirW1,parGroupsIO::Array{<:Real,1})
     "Given a vector of parameters, return the transformed vector that verifies
     an identification condition. Do not do anything if the model is identified."
     #set the first of the in parameters equal to one (Restricted version)
     N,parI,parO = splitVec(parGroupsIO)
     parIO = [ parI - parI[1] ;  parO + parI[1] ]
     return parIO
   end

 function estimate(Model::NetModelDirW1; strIO::Array{<:Real,1} = Model.obs,
                     groupsInds = Model.groupsInds, targetErr  =targetErrValStaticNetsW,bigConst   = bigConstVal)
     "Given model type, observations of the statistics and groups assignments,
     estimate the parameters."
     N,strI,strO = splitVec(strIO)
     groupsNames = unique(groupsInds)
     NG = length(groupsNames)
     tmpDic = countmap(groupsInds)
     #number of nodes per group
     NgroupsMemb = sortrows([[i for i in keys(tmpDic)] [  tmpDic[i] for i in keys(tmpDic) ] ])
     # agrregate degree of each group
     strGroups = zeros(Float64,2NG) ; for i =1:N strGroups[groupsInds[i]] += strI[i]; strGroups[NG + groupsInds[i]] += strO[i]; end
     ~,strGroupsI , strGroupsO = splitVec(strGroups)
     zerIndI = strGroupsI.==0
     zerIndO = strGroupsO.==0
     #sort the groups for degrees
     indsStrGroupsOrderI = Int.(sortrows([strGroupsI Vector(1:NG)], rev = true)[:,2])
     indsStrGroupsOrderO = Int.(sortrows([strGroupsO Vector(1:NG)], rev = true)[:,2])

     (sum(strGroupsI) - sum(strGroupsO) )< 10^Base.eps()? S =sum(strGroupsO) :error()
     SperLink  = S/(N*N-N)
     # The parameters inside this function are Restricted to be in [0,1] hence
     # they need to be converted (log(par)) in order to be used outside

     #uniform parameter as starting value
     unifstartval = 0.5/SperLink;
     bigConstVal = 1e4
     local parI = unifstartval*ones(Float64,NG)
     local parO = unifstartval*ones(Float64,NG)
     if NG==1
         error()
         parI[1] = 0.0
         parO[1] = log(unifstartval)
         outParUN = [ parI; parO ]
         outMod = NetModelDirW1(Model.obs,outParUN,groupsInds)
         i = 0
     else
         #bigNumber = 1e3./targetErr
         eps =  targetErr*1e-12 #10*Base.eps()

         i=0
         maxIt = 500
         relErr = 1
         tol_x = 10*Base.eps()

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

             for nInd=1:NG
                 ng = indsStrGroupsOrderI[nInd]

                 indGzO= (parO.>0).& (.!zerIndO); indGzO[ng] = false
                 indLzO= (parO.<0).& (.!zerIndO); indLzO[ng] = false
                 lBndParI = maximum( [(sum(indGzO)>0? -minimum(parO[indGzO]): [] ); (sum(indLzO)>0 ? -minimum(parO[indLzO]) : []) ] )

                 #println((1,nInd))
                 #put to zero each component of the gradient1
                 function objfunI_ng(parI_ng::Real)
                     parVecOthGroups = parI_ng .+ parO[1:end .!=ng]
                     parO_ng = parO[ng]
                     expGroupStrI = NgroupsMemb[ng,2]*( sum( NgroupsMemb[1:end .!=ng ,2].*( 1./( parVecOthGroups)) )    + (NgroupsMemb[ng,2]-1)* ( 1/(parI_ng + parO_ng)) )
                     return  strGroupsI[ng] - expGroupStrI

                 end
                 #find zero of given component

                 if strGroupsI[ng] <= 1e-3
                     #error("Remove zero col and rows")
                     parI[ng] = bigConstVal
                     #println(0)
                 else
                     tmpBnd = 1/(strGroupsI[ng]/(NgroupsMemb[ng,2]*(N-1)) )
                     #println((lBndParI))
                     brakets =[ lBndParI + 1/S ;  bigConstVal]
                     parI[ng] = tmpFindZero(objfunI_ng,brakets,parI[ng])
                 end
             end
             parUnIO = ([parI;parO])
             #println(parI)
             err = firstOrderCond(Model; strIO = strIO,parGroupsIO = parUnIO , groupsInds = groupsInds)
             #display(err)


             for nInd=1:NG
                 ng = indsStrGroupsOrderO[nInd]

                 indGzI=(parI.>0).& (.!zerIndI); indGzI[ng] = false
                 indLzI=(parI.<0).& (.!zerIndI); indLzI[ng] = false
                 lBndParO = maximum( [(sum(indGzI)>0? -minimum(parI[indGzI]): [] ); (sum(indLzI)>0 ? -minimum(parI[indLzI]) : []) ] )

                 #println((2,nInd))
                 #put to zero each component of the gradient1
                 function objfunO_ng(parO_ng::Real)
                     parVecOthGroups = parO_ng .+ parI[1:end .!=ng]
                     parI_ng = parI[ng]
                     expGroupStrO = NgroupsMemb[ng,2]*( sum( NgroupsMemb[1:end .!=ng ,2].*(1./( parVecOthGroups)) )    + (NgroupsMemb[ng,2]-1)* (1/(parI_ng + parO_ng)) )
                     return  strGroupsO[ng] - expGroupStrO
                 end
                 #find zero of given component
                 if strGroupsO[ng]<= 1e-3
                     #error("Remove zero col and rows")
                     parO[ng] = bigConstVal
                     #println(bigConstVal)
                     #println(0)
                 else
                     tmpBnd = 1/(strGroupsO[ng]/(NgroupsMemb[ng,2]*(N-1)) )
                     #println((lBndParO))
                     brakets =[lBndParO + 1/S ;  bigConstVal ]
                     parO[ng] = tmpFindZero(objfunO_ng,brakets,parO[ng])
                 end

             end
             parUnIO = ([parI;parO])
             #println(parO)
             err = firstOrderCond(Model; strIO = strIO,parGroupsIO = parUnIO , groupsInds = groupsInds)
             #display(err)
         #    tmpMat = par.*par';expMat = putZeroDiag!(tmpMat./(1+tmpMat)) ;expDeg = squeeze(sum(expMat,2),2)
             relErrVec = abs.(err./strGroups)[strGroups .>= 1e-3]
             #println(relErrVec)
             #take into account also groups with zero strength. multiply by target error
             # implies that the expected strength for groups with zero strength
             # has to be lower than 1
             absErrRescVec = abs.(err[strGroups.==0]) .* targetErr
             relErr = maximum([relErrVec;absErrRescVec])
             isnan(relErr)?error():()
             # display(err')
              #println((mean(relErrVec[deg.!=0]),relErr))
             i+=1
             #i==1?error():()
         end
         #println([parI parO])
         i<maxIt?():(println(i);error())
         # the parameter of the first index has to be unique because I want it to be the only infinite
         strIO[1]==0 ? parI[1] = bigConstVal*50:()
         parIO = [parI;parO]
         #println(unique(parIO))
         α = exp(Model.Par[2])
         parIOout = parIO #.* α

         outPar = identify!(Model,parIO)
         outMod = NetModelDirW1(strIO,[outPar,Model.Par[2]],groupsInds)
     end
     return outPar, i , outMod
  end
 function sampl(Model::NetModelDirW1,Nsample::Int;  parArray::Array{Array{<:Real,1}}=Model.Par, groupsInds = Model.groupsInds )
     parGroupsIO = parArray[1]
     NGcheck,parI,parO = splitVec(parGroupsIO)
     parNodesIO = [parI[groupsInds];parO[groupsInds]]
     expMat = expMatrix(Model, parNodesIO)
     N = length(expMat[:,1])
     SampleMats =zeros(Float64,N,N,Nsample)
     logα = parArray[2][1]
     α = exp(logα)
     for s=1:Nsample
         SampleMats[:,:,s] = samplSingMatCan(Model,expMat,α)
     end
     return SampleMats
  end

#--------------------------- Snapshot sequence functions

function zeroStrParFun(Model::NetModelDirW1Afixed,N::Int;method="SMALL" ,smallProb = 1e-1)
    #define a number big enough to play the role of Inf for
    # purposes of sampling N(N-1) bernoulli rvs with prob 1/(1+exp(  bigNumb))
    # prob of sampling degree > 0 (N-1)/(1+exp(bigNumb)) < 1e-6
    if uppercase(method) == "INF"
        # zero degrees have -Inf fitnesses
        zeroDegPar = - Inf #log((N-1)/smallProb - 1)
    elseif uppercase(method) == "10AS21"
        # such that the difference  fit(2) - fit(1) == fit(1)-fit(0) , where
        # fit(n) means the fitness associated with a node of deg n.
        # when more than one choice is possible (the fitness depends on both degrees)
        # than choose fit(2) and fit(1) that maximize the difference

    elseif uppercase(method) == "SMALL"
        # such that the probability of observing a link is very small
        zeroDegPar = 100
    end

    return zeroDegPar
 end

function estimate(model::NetModelDirW1Afixed; strIO::Array{<:Real,1} = Model.obs,
                Amat::BitArray{2}=falses(10,10) , targetErr  =0.001,identIter = false)
    ## Write and test a new way to estimate DirW1 extending the iteration of
    # Chatterjee, Diaconis and  Sly to weighted networks
    N,strI,strO = splitVec(strIO)
    sum(Amat) == 0 ?  Amat = (strI.>0).*(strO'.>0) :()
    lstrI = log.(strI);lstrO = log.(strO)
    (sum(strI) - sum(strO) )< 10^Base.eps()? S =sum(strO) :error()
    SperLink  = S/(N*N-N)
    unifstartval =  -0.5*log(SperLink)
    nnzInds = strIO.!=0
    #inds = (ldegI.==0)
    #function that defines the iteration
    function phiFunc(vParUn::Array{<:Real,1})
        ~,vParUnI,vParUnO = splitVec((vParUn))
        matI = exp.( -vParUnI .+ vParUnO' )
        matO =  1./matI
        mat2sumI =  putZeroDiag(1./(1+matI ))

        mat2sumI[.!Amat ] =0
        mat2sumO =  putZeroDiag(1./(1+matO ))
        mat2sumO[.!Amat ] =0
    #    println(putZeroDiag(1./(1+matI )))
        outPhiI = ( - lstrI + log.( sumSq( mat2sumI ,2) ) )
        outPhiO = ( - lstrO + log.( sumSq( mat2sumO,1) ) )
        #println(vParReI[ldegI.==0])
        #println(log.(sumSq(matI,2))[inds])
        #println(matI[inds,:][1,1:20])
        #println(matI[inds,:][2,1:20])
        return [outPhiI;outPhiO]
     end
    function checkFit(vParUn::Vector{<:Real})
        #println(vParUn)
        expMat = expMatrix(model,vParUn,Amat)

        expStrIO =  [sumSq(expMat,2) ; sumSq(expMat,1) ]
        errsIO = (expStrIO - strIO)[nnzInds]
        relErrsIO = abs.(errsIO)./(strIO[nnzInds])
        # mInd = indmax(relErrsIO)
        # println((mInd,relErrsIO[mInd],strIO[nnzInds][mInd],expStrIO[nnzInds][mInd],vParUn[nnzInds][mInd] ))
        # println(sum(expMat[mInd,:]) )
        # tmp = Vector(1:2N)
        # println(tmp[nnzInds][mInd])

        return prod(relErrsIO .< targetErr)
     end

    maxIt = 200
    it=0
    it2 = 0
    vParUn = [-lstrI;-lstrO]# 0*unifstartval.*ones(2N)
    #starting from an identified set of par. and then identify at each step
    #amounts to doing an identified step. i.e. moving summing an identified
    #vector
    identIter?vParUn =  identify!(model,vParUn):()
    while (!checkFit(vParUn)) & (it<maxIt)
        println(it)
        vParUn = phiFunc(vParUn); it+=1
        identIter?vParUn =  identify!(model,vParUn):()
    end
    #if did not obtain the required precision try with the slower method
    vParUn[.!nnzInds] = zeroStrParFun(model,N)
    if it==maxIt

        startSlowVals = exp.(vParUn)
        outSlow = estimateSlow(model; strIO = strIO,targetErr = targetErr,
                            startVals =startSlowVals)

        vParUn = identify!(model,outSlow[1])
        it2 = outSlow[2]
    end
    #zeroInds = degIO.==0
    #vParUn[zeroInds] = bigNegParFun(model,N)
    return vParUn,it,it2
 end
 #

struct SnapSeqNetDirW1
    obsT:: Array{<:Real,2}  # In and Out strengths for each t . The last
                            # dimension is always time. 2NxT
    Par::Array{<:Real,2} # In and Out Parameters [InW_T;outW_T] for each t
 end
SnapSeqNetDirW1(obs::Array{<:Real,2}) = SnapSeqNetDirW1(obs,zeros(obs))
fooSnapSeqNetDirW1 = SnapSeqNetDirW1(zeros(2,3))
function estimate(Model::SnapSeqNetDirW1; strIO_T::Array{<:Real,2}=Model.obsT,
                        targetErr::Real=1e-3,identPost=false,identIter=false,zeroDegFit="small")
   N2,T = size(strIO_T);N = round(Int,N2/2)
   UnParT = zeros(2N,T)
   prog = 0
   for t = 1:T
       str_t = strIO_T[:,t]
       UnPar_t , Niter1,Niter2  = estimate(fooNetModelDirW1Afixed; strIO = str_t ,Amat = falses(10,10), targetErr =  targetErr,identIter=identIter)
       identPost?UnParT[:,t] = identify!(fooNetModelDirW1Afixed,UnPar_t): UnParT[:,t] = UnPar_t
        round(t/T,1)>prog? (prog = round(t/T,1);println(prog)):()
        #println((t,Niter1,Niter2))
   end
   UnParT[strIO_T.==0] = zeroDegParFun(fooNetModelDirW1Afixed,N;method=zeroDegFit )
 return UnParT
 end

 function sampl(Model::SnapSeqNetDirW1, vParUnIO_T::Array{<:Real,2}; Nsample::Int=1  )
    N2,T = size(vParUnIO_T);N = round(Int,N2/2)
    SampleMats =zeros(Float64,Nsample,N,N,T)
    for s=1:Nsample
        for t = 1:T
            SampleMats[s,:,:,t] = samplSingMatCan(fooNetModelDirW1Afixed,expMatrix(fooNetModelDirW1Afixed,vParUnIO_T[:,t]))
        end
    end
    Nsample == 1 ? (return squeeze(SampleMats,1)):(return SampleMats )
  end


 #---------------------------------Binary DIRECTED Networks with predetermined COVARIATES



end































##

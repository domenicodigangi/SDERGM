
#__precompile__()
"Various functions regarding Dynamic Binary Networks and fitness models"
module DynNetsW
#export sample_bern_mat, sample_fitTV, sampleDgpNetsTv

using Distributions, StatsBase,Optim, LineSearches, StatsFuns, Roots
using StaticNets,Utilities
using ForwardDiff

abstract type GasNetModel end
abstract type GasNetModelW <: GasNetModel end
abstract type GasNetModelWcount <: GasNetModelW end
abstract type GasNetModelBin <: GasNetModel end

#constants

#-----------------WEIGHTED DIRECTED NETWORKS------------


# Directed NEtworks
struct  GasNetModelDirW1 <: GasNetModelW
     """ A Logistic Gas model for Directed Weighted networks and probability depending only
             on time varying parameters.
             ''\p^{(t)}_{ij}  = logistic( \theta^{(t)}_i + \theta^{(t)}_j)''

             # Model Definition and parameters conversion functions
         """
     # node - specific W and group Specific (A B)
     obsT:: Array{<:Real,2} # In strengths degrees and out strenghts for each t . The first dimension is always time. Tx2N
     Par::Array{Array{<:Real,1},1} #[[InW;outW], Agroups , Bgroups]  Group specific parameters
     #Nodes in a group have the same GAS parameters, but different fitnesses

     # each node belongs to a W group, each W group is assigned 2 AB groups, in order
     #to allow different parameters for in and our. The AB groups In and Out can be the same

     groupsInds::Array{Array{Int,1},1} #[[Nodes_2_groupsW, Wgroups_2_groupsABI , Wgroup_2_groupsABO ]  (_2_ stands for "assignement to" )
     # In and out W parameters are assumed to be always different, while In and Out AB
     # can be equal. Moreover, I assume that nodes with the same Wgroup have the same ABgroups
     scoreScalingType::String # String that specifies the rescaling of the score. For a list of possible
     # choices see function scalingMatGas
  end

#inizializza senza specificare assenza di scaling
GasNetModelDirW1( obsT , Par  , groupsInds) =
        GasNetModelDirW1( obsT,Par, groupsInds,"FISHER-DIAG")
#Initialize by observations only
GasNetModelDirW1(obsT:: Array{<:Real,2}) =   (N = round(length(obsT[1,:])/2);
                                                GasNetModelDirW1(obsT,
                                                    [ zeros(2N) ,0.9  * ones(1), 0.01 * ones(1)],
                                                    [Int.(1:N),Int.(ones(N)),Int.(ones(N))],"FISHER-DIAG") )
fooGasNetModelDirW1 = GasNetModelDirW1(ones(100,20))

# Relations between Static and Dynamical models: conventions on storage for
# parameters and observations
StaModType(Model::GasNetModelDirW1) = StaticNets.fooNetModelDirW1# to be substituted with a conversion mechanism
# linSpacedPar(Model::GasNetModelDirBin1,Nnodes::Int;NgroupsW = Nnodes, deltaN::Int=3,graphConstr = true) =
#         StaticNets.linSpacedPar(StaModType(Model),Nnodes;Ngroups = NgroupsW,deltaN=deltaN,graphConstr =graphConstr);
#
# options and conversions of parameters for optimization

function setOptionsOptim(Model::GasNetModelDirW1)
    "Set the options for the optimization required in the estimation of the model.
    For the optimization use the Optim package."
    tol = eps()*10
    opt = Optim.Options(  g_tol = tol,
                     x_tol = tol,
                     f_tol = tol,
                     iterations = 5000,
                     show_trace = true,#false,#
                     show_every=1)

    algo = Newton(; alphaguess = LineSearches.InitialHagerZhang(),
                     linesearch = LineSearches.BackTracking())

    # algo = NewtonTrustRegion(; initial_delta = 0.5,
    #                 delta_hat = 1.0,
    #                 eta = 0.1,
    #                 rho_lower = 0.25,
    #                 rho_upper = 0.75)# Newton(;alphaguess = Optim.LineSearches.InitialStatic(),
      return opt, algo
 end


function NumberOfGroups(Model::GasNetModelDirW1,groupsInds::Array{Array{Int,1},1})
     NGW = length(unique(groupsInds[1]))
     GBA = length( unique( [groupsInds[2][.!(groupsInds[2].==0)] ; groupsInds[3][.!(groupsInds[3].==0)] ]) )
      return NGW, GBA
 end
function NumberOfGroupsAndABindNodes(Model::GasNetModelDirW1,groupsInds::Array{Array{Int,1},1})
     NGW = length(unique(groupsInds[1]))
     GBA = length( unique( [groupsInds[2][.!(groupsInds[2].==0)] ; groupsInds[3][.!(groupsInds[3].==0)] ]) )
     ABgroupsIndNodesIO = [groupsInds[2][groupsInds[1]] ; groupsInds[3][groupsInds[1]] ]
     indTvNodesIO =   (.!(ABgroupsIndNodesIO .==0))#if at least one group ind is zero then both in and out tv par are constant
     return NGW, GBA, ABgroupsIndNodesIO, indTvNodesIO
 end
function vec2ArrayGasPar(Model::GasNetModelDirW1,VecGasPar::Array{<:Real,1})
     groupsInds = Model.groupsInds
     NGW, GBA = NumberOfGroups(Model,groupsInds)
     ArrayGasPar =  [VecGasPar[1:2NGW],ones(GBA).*VecGasPar[2NGW+1:2NGW+GBA],ones(GBA).*VecGasPar[2NGW+GBA+1:2NGW + 2GBA]]
     return ArrayGasPar
     end
function restrictGasPar(Model::GasNetModelDirW1,vecUnGasPar::Array{<:Real,1})
     "From the Unrestricted values of the parameters return the restricted ones
     takes as inputs a vector of parameters that can take any value in R and returns
     the appropriately contrained values, e.g. the coefficient of autoregressive
     component in the GAS update rule has to be positive and <1."
     groupsInds = Model.groupsInds
     NGW, GBA= NumberOfGroups(Model,groupsInds)
     W_un = vecUnGasPar[1:2NGW]
     diag_B_Un = vecUnGasPar[2NGW+1:2NGW+GBA]
     diag_A_Un = vecUnGasPar[2NGW+ GBA+1:2NGW+2GBA]
     W_Re = W_un
     diag_A_Re = exp.(diag_A_Un)
     diag_B_Re = 1./(1+exp.(.-diag_B_Un))
     vecRePar =  [W_Re; diag_B_Re; diag_A_Re]
     return vecRePar
     end
function unRestrictGasPar( Model::GasNetModelDirW1,vecReGasPar::Array{<:Real,1})
     "From the restricted values of the parameters return the unrestricted ones
     takes as inputs a vector of parameters that can take any value in R and returns
     the appropriately contrained values, i.e. the coefficient of autoregressive
     component in the GAS update rule has to be positive and <1."
     groupsInds = Model.groupsInds
     NGW, GBA= NumberOfGroups(Model,groupsInds)
     W_Re = vecReGasPar[1:2NGW]
     diag_B_Re = vecReGasPar[2NGW+1:2NGW+GBA]
     diag_A_Re = vecReGasPar[2NGW+GBA+1:end]

     W_Un = W_Re
     diag_B_Un = log.(diag_B_Re./(1-diag_B_Re ))
     diag_A_Un = log.(diag_A_Re)

     vecUnPar =  [W_Re;diag_B_Un; diag_A_Un]
     return vecUnPar
     end


function uBndPar2bndPar_T(Model::GasNetModelDirW1, uBndParNodesIO_T::Array{<:Real,2}  )
    T,N = size(uBndParNodesIO_T)
    outBndPar_T = zeros(uBndParNodesIO_T)
    for t = 1:T
        uBndParNodesIO_t = uBndParNodesIO_T[t,:]
        outBndPar_T[t,:] = StaticNets.uBndPar2bndPar(StaModType(Model), uBndParNodesIO_t )
    end
    return outBndPar_T
end
#Gas Filter Functions

function scalingMatGas(Model::GasNetModelDirW1,expMat::Array{<:Real,2},I_tm1::Array{<:Real,2})
    "Return the matrix required for the scaling of the score, given the expected
     matrix and the Scaling matrix at previous time. "
    if uppercase(Model.scoreScalingType) == "FISHER-EWMA"
         error()
        # λ = 0.8
        #
        # I = expMat.*(1-expMat)
        # diagI = sum(I,2)
        # [I[i,i] = diagI[i] for i=1:length(diagI) ]
        # I_t =  λ*I + (1-λ) *I_tm1
        # scaling = sqrtm(I) ##
    else
        error()
    end
    return scaling
 end

function scaleScoreGas(Model::GasNetModelDirW1,  ftotIO_t::Array{<:Real,1},fBndTot_t::Array{<:Real,1},gradIO_t::Array{<:Real,1}, expfIO_t::Array{<:Real,1},exp_Mat_t::Array{<:Real,2},
                        indMin::Int64,  WgasIO::Array{<:Real,1},BgasIO::Array{<:Real,1},AgasIO::Array{<:Real,1},strIO_t::Array{<:Real,1},
                        exp_strIO_t::Array{<:Real,1},t::Int; I_tm1::Array{<:Real,2} = zeros(exp_Mat_t))
    "Return the rescaled score according to the scalingtype of the model"
    if uppercase(Model.scoreScalingType) == ""
        s_t = gradIO_t
    elseif uppercase(Model.scoreScalingType) == "FISHER-DIAG"
        #Rescale the score in order to set its variance to one
        tmpMat = (exp_Mat_t.^2)
        N = round(Int,length(gradIO_t)/2)
        diagI_t = ([sumSq(tmpMat,2) ; sumSq(tmpMat,1)])
        tmpDiagEl =   sum(diagI_t[1:N]) - diagI_t[N+indMin] #
        print_diagI_t = copy(diagI_t)
        if false # unhortodox scaling that uses the observed values
            tmpDiagEl =  ( (sum(diagI_t[1:N]) - diagI_t[N+indMin]) + sum(strIO_t[1:N]) )/2#
            diagI_t = (diagI_t + strIO_t.^2)./2
        end
        diagI_t[N + indMin] = tmpDiagEl
        diagI_t = diagI_t.*( expfIO_t.^2)
        rescVec = sqrt.(diagI_t)
        s_t =  gradIO_t./rescVec
        s_t[gradIO_t.==0] = 0

        if false
            println(t)
        printInd = 45
    display((ftotIO_t[1:N][printInd],fBndTot_t[1:N][printInd],strIO_t[1:N][printInd],exp_strIO_t[1:N][printInd],s_t[1:N][printInd]))
    display((ftotIO_t[N+1:end][printInd],fBndTot_t[N+1:end][printInd],strIO_t[N+1:end][printInd],exp_strIO_t[N+1:end][printInd],s_t[1+N:2N][printInd]))
    display((ftotIO_t[N+1:end][indMin],fBndTot_t[N+1:end][indMin],strIO_t[N+1:end][indMin],exp_strIO_t[N+1:end][indMin],s_t[1+N:end][indMin]))
        end
        if false# maximum(abs.(s_t[isfinite.(s_t)])) > 1e4


            indOprint = indmax(abs.(s_t[1+N:2N]))
            @show indOprint
             println("str")
             println(strIO_t[indOprint])
            # display([strIO_t[1:N] strIO_t[N+1:2N] ] )
             println("exp Str")
             println((indmax(exp_strIO_t[1:N]),exp_strIO_t[1:N][indmax(exp_strIO_t[1:N])] ,indmax(exp_strIO_t[N+1:2N]),exp_strIO_t[N+1:2N][indmax(exp_strIO_t[N+1:2N])]  ) )
             #display([exp_strIO_t[1:N] exp_strIO_t[N+1:2N] ] )
            println("fBndTot_t")
            println(fBndTot_t[N+1:end][indOprint])
            #display([fBndTot_t[1:N] fBndTot_t[N+1:2N] ] )
            println("exp_mat_t")
            println(exp_Mat_t[:,indOprint])
            #display(exp_Mat_t)
            println("grad")
             #display([gradIO_t[1:N] gradIO_t[N+1:2N] ] )
           println([maximum(abs.(gradIO_t[2:N])) indmax(abs.(gradIO_t[2:N])); maximum(abs.(gradIO_t[1+N:2N])) indmax(abs.(gradIO_t[1+N:2N]))])
             println(gradIO_t)
             println("diagI_t")
              #display([print_diagI_t[1:N] print_diagI_t[N+1:2N] ] )
              println(" expfIO_t.^2")
               #display([ (expfIO_t.^2)[1:N]  (expfIO_t.^2)[N+1:2N] ] )
              println("tmpDiagEl")
              #display(tmpDiagEl)
            println("resc vect")
            println(rescVec)
            #display([rescVec[1:N] rescVec[N+1:2N]])
           println("score")
           println(s_t[2:end])
           println([maximum(abs.(s_t[2:N])) indmax(abs.(s_t[2:N])); maximum(abs.(s_t[1+N:2N])) indmax(abs.(s_t[1+N:2N]))])
            # println(s_t[1:N])
            #println(s_t[N+1:2N]  )
            #display([s_t[1:N] s_t[N+1:2N] ] )
            println("A")
            # display([AgasIO[1:N] AgasIO[N+1:2N] ] )
             println("B")
             # display([BgasIO[1:N] BgasIO[N+1:2N] ] )
              println("W")
              # display([WgasIO[1:N] WgasIO[N+1:2N] ] )
            maximum(abs.(s_t[isfinite.(s_t)])) > 1e8?  error():()
        end

    elseif uppercase(Model.scoreScalingType) == "FISHER-EWMA"
         error()
        # λ = 0.8
        #
        # I = scalingMatGas(Model,grad_t,exp_mat_t,I_tm1)
        # diagI = sum(I,2)
        # [I[i,i] = diagI[i] for i=1:length(diagI) ]
        # I_t =  λ*I + (1-λ) *I_tm1
        # scalingMat = sqrtm(I) ##
    end

    return s_t, I_tm1
 end
function predict_score_driven_par( Model::GasNetModelDirW1,strIO_t::Array{<:Real,1},
                         ftotIO_t::Array{<:Real,1},I_tm1::Array{<:Real,2},indTvNodesIO::BitArray{1},
                         WgasIO::Array{<:Real,1},BgasIO::Array{<:Real,1},AgasIO::Array{<:Real,1},t::Int ;
                         fBndTot_t::Array{<:Real,1} =  StaticNets.uBndPar2bndPar(StaModType(Model),ftotIO_t) )
     #= likelihood and gradients depend on all the parameters (ftot_t), but
     only the time vaying ones (f_t) are to be updated=#
     fIO_t = ftotIO_t[indTvNodesIO] #Time varying fitnesses
     N,ftotI_t,ftotO_t = splitVec(ftotIO_t)
     indTvNodesI = indTvNodesIO[1:N]

     #find the minimum out parameter
     indMin = indmin(ftotO_t)
     # indMin2 = indmax(strIO_t[N+1:end])
     #
     indTvNodesIO[N+1:end][indMin] ? (scoreAdd = true; indMin = indmin(ftotO_t[indTvNodesIO[N+1:end]]) ) : scoreAdd = false
     #display((BgasIO[1],AgasIO[1]))
     sumParMat, exp_mat_t = StaticNets.expMatrix2(StaModType(Model),fBndTot_t )

     exp_strIO_t = [sumSq(exp_mat_t,2);sumSq(exp_mat_t,1)]
     #The following part is a faster version of logLikelihood_t
     loglike_t =   ((sum(log.(sumParMat)) - sum(fBndTot_t.*strIO_t) ))

     # compute the score of observations at time t wrt parameters at tm1
     expfIO_t = exp.(fIO_t)

     strDiffIO_t =   (exp_strIO_t[indTvNodesIO] -strIO_t[indTvNodesIO])

     gradIO_t = strDiffIO_t.*expfIO_t # only for TV paramaters
    # scoreAdd ?  gradIO_t[N + indMin] +=  -expfIO_t[N + indMin] * sum(strDiffIO_t[1:N]) :()

     sIO_t,I_tm1 =  scaleScoreGas(Model,  ftotIO_t,fBndTot_t,gradIO_t, expfIO_t,exp_mat_t,indMin, WgasIO,BgasIO,AgasIO,strIO_t,exp_strIO_t,t)
     sIO_t[1] = 0 #identification of the update equation
     # sembra che cosi facendo la strength In del primo nodo non contriubisca all'update
     # ma non \'e' cosi. Infatti l'info sul primo in degree è ridondante se considerata
     # insieme a tutti gli altri in degrees e tutti gli out degrees


     if false# indMin_tp1 == indMin
         #update the out parameters
         fO_tp1 = WgasIO[N+1:2N] + BgasIO[N+1:2N].* fIO_t[N+1:2N] + AgasIO[N+1:2N].*sIO_t[N+1:2N]
         #if the minimum out par has changed than transform all the in (static and dynamic)
         indMin_tp1 = indmin(fO_tp1)
         @show (indMin,indMin_tp1)
        # ftotIO_t[1:N] = log.(exp.(ftotI_t) + exp(fO_tp1[indMin_tp1]) - exp(ftotO_t[indMin]) )
        fI_t = ftotIO_t[1:N][indTvNodesI]
        # update the in parameters
        fI_tp1 = WgasIO[1:N] + BgasIO[1:N].* fI_t + AgasIO[1:N].*sIO_t[1:N]
        fIO_tp1 = [fI_tp1 fO_tp1]
    elseif false # unhortodox modification that limits exploding scores
        score_thr_pos = 1
        score_thr_neg = -1
        shockod = AgasIO.*sIO_t
        shockod[(shockod).>score_thr_pos] = score_thr_pos.*( 1+ log.(abs.(shockod[(shockod).>score_thr_pos])))
        shockod[(shockod).<score_thr_neg] = score_thr_neg.*( 1+ log.(abs.(shockod[(shockod).<score_thr_neg])))
        #s_t[abs.(s_t).>score_thr] = sign.(s_t[abs.(s_t).>score_thr]).*score_thr.*( 1+ log.(abs.(s_t[abs.(s_t).>score_thr])))
        fIO_tp1 = WgasIO + BgasIO.* fIO_t + shockod
    elseif true # unhortodox modification that limits exploding scores
        score_thr_pos = 1000
        score_thr_neg = -35
        shockod = sIO_t
        shockod[(shockod).>score_thr_pos] = score_thr_pos.*( 1+ log10.(abs.(shockod[(shockod).>score_thr_pos])))
        shockod[(shockod).<score_thr_neg] = score_thr_neg.*( 1+ log10.(abs.(shockod[(shockod).<score_thr_neg])))
        #s_t[abs.(s_t).>score_thr] = sign.(s_t[abs.(s_t).>score_thr]).*score_thr.*( 1+ log.(abs.(s_t[abs.(s_t).>score_thr])))
        fIO_tp1 = WgasIO + BgasIO.* fIO_t + AgasIO.*shockod
    else
        fIO_tp1 = WgasIO + BgasIO.* fIO_t + AgasIO.*sIO_t

    end

     ftotIO_t[indTvNodesIO] = fIO_tp1 #of all parameters udate dynamic ones with GAS
     ftotIO_tp1 = ftotIO_t
     return ftotIO_tp1,loglike_t,I_tm1
  end
function score_driven_filter( Model::GasNetModelDirW1,
                                vResGasPar::Array{<:Real,1} ;
                                obsT::Array{<:Real,2}=Model.obsT,
                                groupsInds:: Array{Array{Int,1},1}=Model.groupsInds,
                                dgpTN = (0,0))
    """GAS Filter the Dynamic Fitnesses from the Observed degrees, given the GAS parameters
     given T observations for the degrees in TxN vector degsT
     """
     #se parametri per dgp sono dati allora campiona il dgp, altrimenti filtra
    sum(dgpTN) == 0 ? dgp = false : dgp = true
    if dgp
        T,N = dgpTN
     Y_T = zeros(Float64,T,N,N)
    else
     T,N2 = size(obsT);N = round(Int,N2/2)
    end

    NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model,groupsInds)
    # ABgroups Inds for each par (In and Out). IO means that we have a 2N vector
    # of indices: one index for each In par and below one index for each Out par


    # Organize parameters of the GAS update equation
    WGroupsI = vResGasPar[1:NGW]
    WGroupsO = vResGasPar[NGW+1:2NGW]
    W_allIO = [WGroupsI[groupsInds[1]] ; WGroupsO[groupsInds[1]] ]


    BgasGroups  = vResGasPar[2NGW+1:2NGW+GBA]
    AgasGroups  = vResGasPar[2NGW+GBA+1:2NGW+2GBA]
    #distribute nodes among groups with equal parameters
    AgasIO = AgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
    BgasIO = BgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
    WgasIO   = W_allIO[indTvNodesIO]# identify!(Model,W_allIO[indTvNodesIO])

    # start values equal the unconditional mean, but  constant ones remain equal to the unconditional mean, hence initialize as:
    UMallNodesIO = W_allIO
    UMallNodesIO[indTvNodesIO] =  WgasIO./(1-BgasIO)
    if false
        println((vResGasPar[2NGW+1:2NGW+2GBA]))
        println("W")
        display([WgasIO[1:N] WgasIO[N+1:2N] ] )
        println("Unc Means")
         display([UMallNodesIO[1:N] UMallNodesIO[N+1:2N] ] )
     end
    fIOVecT = ones(Real,T,2N)

    ftotIO_t = UMallNodesIO# identify!(Model,UMallNodesIO)
    I_tm1 = eye(Float64,2N)
    loglike = 0

    for t=1:T
        #println(t)
        #display(ftotIO_t)
        fBndTot_t =  StaticNets.uBndPar2bndPar(StaModType(Model),ftotIO_t)

        if dgp
            expMat_t = StaticNets.expMatrix(StaModType(Model),fBndTot_t )# exp.(ftot_t)
            Y_T[t,:,:] = StaticNets.samplSingMatCan(StaModType(Model),expMat_t)
            #display(expMat_t)
            strIO_t =  [ sumSq(Y_T[t,:,:],2); sumSq(Y_T[t,:,:],1)]
        else
            strIO_t = obsT[t,:]  # vector of in and out degrees
        end
        # println("str_t")
        #  display([strIO_t[1:N] strIO_t[N+1:2N] ] )
        #140<=t<=150? display(fBndTot_t'):()
        ftotIO_t,loglike_t,I_tm1 = predict_score_driven_par(Model,strIO_t,ftotIO_t,I_tm1,indTvNodesIO,WgasIO,BgasIO,AgasIO,t ; fBndTot_t = fBndTot_t)
        fIOVecT[t,:] = ftotIO_t #store the filtered parameters from previous iteration
        loglike += loglike_t

    end
    if dgp
        return Y_T , fIOVecT
    else
        return fIOVecT, loglike
    end
    end
score_driven_filter(Model::GasNetModelDirW1) = score_driven_filter(Model,array2VecGasPar(Model,Model.Par))
sampl(Mod::GasNetModelDirW1,T::Int)=( N = length(Mod.groupsInds[1]) ;
                                        tmpTuple = score_driven_filter(Mod,[Mod.Par[1];Mod.Par[2];Mod.Par[3]];  dgpTN = (T,N) );
                                        degsIO_T = [sumSq(tmpTuple[1],3) sumSq(tmpTuple[1],2)];
                                        (GasNetModelDirBin1(degsIO_T,Mod.Par,Mod.groupsInds,Mod.scoreScalingType),tmpTuple[1],tmpTuple[2]) )

# # Estimation
#
function estSingSnap(Model::GasNetModelDirW1, str_t::Array{<:Real,1}; groupsInds = Model.groupsInds, targetErr::Real=targetErrValDynNets)
    hatBndPar,~,~ = StaticNets.estimate(StaModType(Model); strIO = str_t ,groupsInds = groupsInds[1], targetErr =  targetErr)
    hatUbndPar = StaticNets.bndPar2uBndPar(StaModType(Model),hatBndPar)
    return hatUbndPar
 end
function estimateSnapSeq(Model::GasNetModelDirW1; strT::Array{<:Real,2}=Model.obsT,
                     groupsInds = Model.groupsInds,targetErr::Real=targetErrValDynNets,print_t = false)
    T,N2 = size(Model.obsT);N = round(Int,N2/2)
    NGW = length(unique(groupsInds[1]))
    UnParT = zeros(T,2NGW)
    prog = 0
    for t = 1:T
        print_t?println(t):()
        str_t = strT[t,:]
        UnParT[t,:] = estSingSnap(Model,str_t; groupsInds = groupsInds ,targetErr =  targetErr)
        round(t/T,2)>prog ? (prog=round(t/T,2);println(prog) ):()
    end
    return UnParT
 end

function estimate(Model::GasNetModelDirW1;start_values = [zeros(Real,10),zeros(Real,3),zeros(Real,3) ],
                    targeting = false, meanEstSS::Array{Float64,1} = zeros(1) )
    T,N2 = size(Model.obsT);N = round(Int,N2/2)
    groupsInds = Model.groupsInds
    NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model,groupsInds)
    optims_opt, algo = setOptionsOptim(Model)

    #set the starting points for the optimizations if not given
    if sum([sum(start_values[i]) for i =1:3 ])==0

         B0_ReGroups =  0.8 * ones(GBA)
         A0_ReGroups = 0.01 * ones(GBA)

         BgasNodes_0   = B0_ReGroups[ABgroupsIndNodesIO[indTvNodesIO]]
         AgasNodes_0   = A0_ReGroups[ABgroupsIndNodesIO[indTvNodesIO]]
         # choose the starting points for W in order to fix the Unconditional Means
         # equal to the static estimates for the average degrees(over time)
         MeanStrT =  meanSq(Model.obsT,1)
         if targeting
             if sum(meanEstSS) == 0
                 estSS =  DynNets.estimateSnapSeq(Model)
                 meanEstSS = meanSq(estSS,1)
             end
             UMstaticGroups = meanEstSS
         else
             error("This model doe not work well on real data without a targeting")
             UMstaticGroups = StaticNets.bndPar2uBndPar(StaModType(Model), estSingSnap(Model ,MeanStrT; groupsInds = groupsInds ))
         end
        # display([UMstaticGroups[1:N] UMstaticGroups[N+1:2N] ])
         conv_flag1 = true

         GBA>0?( W0_Groups = UMstaticGroups.*(1-B0_ReGroups[[groupsInds[2];groupsInds[3]]])):error()#W0_Groups = UMstaticGroups
         #W0_Groups = W_prob

    else
         W0_Groups = start_values[1]
         B0_ReGroups = start_values[2]
         A0_ReGroups = start_values[3]
         conv_flag1 = true
    end
    if GBA == 0  return [UMstaticGroups, zeros(1), zeros(1)] ,  conv_flag1 end;
    # Estimate  jointly all  Gas parameters
    vGasParAll0_Un = unRestrictGasPar(Model,[W0_Groups ;B0_ReGroups;A0_ReGroups ] )[2:end] # vGasParGroups0_Un#
    # objective function for the optimization
    function objfunGas(vparUn::Array{<:Real,1})# a function of the groups parameters
        #vParRes = restrictGasPar(Model,[-Inf;vparUn])
        GBA,Bun,Aun = splitVec(vparUn) #ottimizza solo A e B
        Bre = 1./(1+exp.(-Bun))
        Are = exp.(Aun)
        W_new = UMstaticGroups.*(1-Bre)
        vParRes = [W_new;Bre;Are]
        #display((vParRes[2],vParRes[end-1],vParRes[end]))
            foo,loglikelValue = score_driven_filter( Model,vParRes)

        return  -loglikelValue
    end

    #Run the optimization without AD if a matrix sqrt is required in the filter
    if uppercase(Model.scoreScalingType) == "FISHER-EWMA"
        ADobjfunGas = objfunGas
    else
        ADobjfunGas = TwiceDifferentiable(objfunGas, vGasParAll0_Un[end-1:end]; autodiff = :forward);
    end
    display(vGasParAll0_Un)

    display(objfunGas(vGasParAll0_Un[end-1:end]))

    #display(ForwardDiff.gradient(objfunGas, vGasParAll0_Un[end-1:end]))
    #error()
    optim_out2  = optimize(objfunGas,vGasParAll0_Un[end-1:end],algo,optims_opt)
    BAestUn = Optim.minimizer(optim_out2)
    BAestRe = [1./(1+exp.(-BAestUn[1:GBA]));exp.(BAestUn[1+GBA:end])]
    arrayGasParHat_Re  = [UMstaticGroups.*(1 - BAestRe[1:GBA]) , BAestRe[1:GBA],BAestRe[GBA+1:end]]
    conv_flag = conv_flag1*Optim.converged(optim_out2)
    println(optim_out2)
    return  arrayGasParHat_Re, conv_flag,Optim.minimum(optim_out2)
    end


end

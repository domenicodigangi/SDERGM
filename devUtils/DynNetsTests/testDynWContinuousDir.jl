N = 10
NG = N#5
NBA = 1
# Set the values of the unconditional means
maxW = 200000 
minW = 10000
tmpStrIO = round.([ Vector(linspace(minW,maxW,N));Vector(linspace(minW,maxW ,N))])
N, StrI,StrO = splitVec(tmpStrIO);sum(StrI) == sum(StrO)?():error()
dgpStrIO  =   StaticNets.DegSeq2graphDegSeq(StaticNets.fooNetModelDirWCount1,tmpStrIO)
groupsInds = [Utilities.distributeAinVecN(Int.(1:NG),N) , Utilities.distributeAinVecN(Int.(1:1),NG) , Utilities.distributeAinVecN(Int.(1:1),NG)]
StaMod = StaticNets.NetModelDirW1(Int.(dgpStrIO),[zeros(Float64,2N),zeros(Float64,1)],groupsInds[1])
estPar,estIt,estMod = StaticNets.estimate(StaMod)
uBndPar = StaticNets.bndPar2uBndPar(StaMod,estPar)
# Test Dgp
A = ones(1)*0.1
B = ones(1)*0.9
W = uBndPar.*(1-B)
T =300
dgpParArr = [W,B,A];dgpParVec = [W;B;A]
scalingType = "FISHER-DIAG"
mod2Samp = GasNetModelDirW1(ones(T,2N), [ zeros(2N) ,0.9  * ones(1), 0.01 * ones(1)], groupsInds,scalingType)
NumberOfGroupsAndABindNodes(mod2Samp,groupsInds)


Y_T, Fitness_TIO = score_driven_filter(mod2Samp,[W;B;A];dgpTN = (T,N))
StrI = sumSq(Y_T,3) ; StrO = sumSq(Y_T,2); StrIO = [StrI  StrO]
Mod2est = GasNetModelDirW1(StrIO,dgpParArr, groupsInds,scalingType)
#test filter
Fitness_TIO_Fil, loglike = score_driven_filter(Mod2est,dgpParVec; groupsInds = groupsInds)
Fitness_TIO - Fitness_TIO_Fil
estimateSnapSeq(Mod2est;print_t = true)

hatBndPar,x,xx = StaticNets.estimate(StaticNets.fooNetModelDirW1;strIO = StrIO[19,:],groupsInds = groupsInds[1] )
hatUbndPar = StaticNets.bndPar2uBndPar(StaticNets.fooNetModelDirW1,hatBndPar)

(hatBndPar[1:N] .+ hatBndPar[N+1:2N]' .>=0)
minInd = indmin(hatBndPar[N+1:2N])
hatBndPar[1+N:2N][minInd] .+ hatBndPar[1:N]
## Test Estimate

estimate(Mod2est)
##plots
1;
    N,Fitness_TI,Fitness_TO = splitMat(Fitness_TIO)
    Fitness_TO[1:5,:]
    using Plots
    plotly()
    pparI = plot(Fitness_TI[:,2:end])
    pparO = plot(Fitness_TO);
    plot(pparI,pparO,layout = (2,1),legend=:none,size = (1200,600))
    pstrI = plot(StrI,title = "  Unc Means = $(round(dgpStrIO[1:N]))");
    pstrO = plot(StrO,title = "  Unc Means = $(round(dgpStrIO[N+1:2N]))");
    plot(pstrI,pstrO,layout = (2,1),legend=:none,size = (1200,600))
    dgpStrIOmat =[dgpStrIO[1:N] dgpStrIO[N+1:2N]]
    meanStrIO =  [meanSq(StrI,1)  meanSq(StrO,1)]
    [meanStrIO  dgpStrIOmat ]
    (meanStrIO .- dgpStrIOmat )./dgpStrIOmat


##
A_prob = 0.00123004
B_prob=  0.840832
W_prob_tmp=    [-Inf         -324.547
                5.33316     -0.274072
               -0.412405    -0.378851
               -0.396107    -0.464246
               -0.357238    -0.401917
               -0.41319     -0.428446
               -0.428308    -0.997147
               -0.496042    -0.361475
               -0.433907    -0.36438
               -0.499994    -0.370849]
W_prob = [W_prob_tmp[1:N];W_prob_tmp[N+1:2N]]
Umean_prob =StaticNets.uBndPar2bndPar(StaMod, W_prob./(1-B_prob))
display(StaticNets.expMatrix(StaMod,Umean_prob))
Y_T, Fitness_TIO =score_driven_filter(Mod2est,[W_prob;B_prob;A_prob]; groupsInds = groupsInds)


##
MeanStrT =  meanSq(Mod2est.obsT,1)
[MeanStrT[1:N] MeanStrT[N+1:2N]]
staMod =  StaticNets.NetModelDirW1((MeanStrT),[zeros(Float64,2N),zeros(Float64,1)],groupsInds[1])
StaticNets.estimate(staMod)
MeanStrT[N+1:end] = -MeanStrT[N+1:end]
sum(MeanStrT)
##










##

##

















#

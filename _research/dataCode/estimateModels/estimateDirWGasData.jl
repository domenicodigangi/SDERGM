#script
## Binary  estimates
    # estSS_1GW,estGas_1GW = 0, 0
    # eMidMod_1GW = DynNets.GasNetModelDirBin1(degsIO_T, [ zeros(2N) ,0.9  * ones(1), 0.01 * ones(1)], [Int.(ones(N)),Int.(ones(1)),Int.(ones(1))],"")
    # @time estSS_1GW =  DynNets.estimateSnapSeq(eMidMod_1GW)
    # @time estGas_1GW = DynNets.estimate(eMidMod_1GW)
    # eMidMod = DynNets.GasNetModelDirBin1(degsIO_T)
    # estSS,estGas,estGasTarg = 0, 0,0
    # estSS =  DynNets.estimateSnapSeq(eMidMod)
    # #estGasTarg = DynNets.estimateTargeting(eMidMod)
    # estGas = DynNets.estimate(eMidMod)
    # # save Bin
    # save_path = save_fold*file_nameStart* "_Bin1_" *file_nameEnd#
    # @save(save_path,estSS,estGas,estGasTarg,estSS_1GW,estGas_1GW,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)
    #
## Weighted Estimates


estSSW1 =  DynNets.estimateSnapSeq(eMidModW1)
estMean = DynNets.estSingSnap(eMidModW1,meanSq(strIO_T,1))

prod(isfinite.(estSSW1[:,2:end]),2)
parest_t=StaticNets.estimate(fooNetModelDirW1; strIO = strIO_T[20,:] ,groupsInds = eMidModW1.groupsInds[1])
strIO_T[20,N+1:2N]
unique(parest_t[1])
prod(parest_t[1][N+1:2N] .> 0)
indmin(parest_t[1][N+1:2N])
log( parest_t[1][N+1]  + parest_t[1][1:N] )
StaticNets.bndPar2uBndPar(fooNetModelDirW1,parest_t[1])

using StatsFuns
smallVal = 1e-15
log(bigConstVal - smallVal)  - log(bigConstVal) + log1p(smallVal/bigConstVal)


##
save_path = save_fold*file_nameStart* "_W1_" *file_nameEnd#
@load(save_path,estSSW1,estGasW1,estGasTargW1,estSSW1_1GW,estGasW1_1GW,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)

# test

B = ones(1)* 0.15999999999999998
A = ones(1)*0.0003000000000000001
uBndPar = StaticNets.bndPar2uBndPar(StaticNets.fooNetModelDirW1, meanSq( DynNets.uBndPar2bndPar_T(eMidModW1,estSSW1),1) )
#StaticNets.bndPar2uBndPar(StaModType(eMidModW1),meanSq(estSSW1,1))
W = uBndPar.*(1-B)
dgpParArr = [W,B,A];dgpParVec = [W;B;A]
filterParTestUbnd,testlike = DynNets.score_driven_filter_or_dgp(eMidModW1,dgpParVec)


filterParTest = DynNets.uBndPar2bndPar_T(eMidModW1,filterParTestUbnd)
 pin = plot(filterParTest[:,1:N][:,meanSq(filterParTest[:,1:N],1).<1e5]);
 pout = plot(filterParTest[:,1+N:2N][:,meanSq(filterParTest[:,1+N:2N],1).<1e5]);
 plot(pin,pout,layout=(2,1),size = (1200,600))

##
targMean =meanSq( estSSW1[1:Tweeks,:],1)# StaticNets.bndPar2uBndPar(StaticNets.fooNetModelDirW1, meanSq( DynNets.uBndPar2bndPar_T(eMidModW1,estSSW1),1) )

 estGasTargW1,x,fVal = DynNets.estimate(eMidModW1;targeting = true, meanEstSS = targMean)
 println(estGasTargW1[2],estGasTargW1[3])
 estGasW1,estSSW1_1GW,estGasW1_1GW = 0,0,0

estParVec = [estGasTargW1[1];estGasTargW1[2];estGasTargW1[3]]
 filterParTestUbndEst,testlike = DynNets.score_driven_filter_or_dgp(eMidModW1,estParVec)
 filterParTest = filterParTestUbndEst#DynNets.uBndPar2bndPar_T(eMidModW1,filterParTestUbndEst)#
 pin = plot(filterParTest[:,1:N][:,meanSq(filterParTest[:,1:N],1).<1e5]);
 pout = plot(filterParTest[:,1+N:2N][:,meanSq(filterParTest[:,1+N:2N],1).<1e5]);
 plot(pin,pout,layout=(2,1),size = (1200,600))

#@save(save_path,estSSW1,estGasW1,estGasTargW1,estSSW1_1GW,estGasW1_1GW,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)

##

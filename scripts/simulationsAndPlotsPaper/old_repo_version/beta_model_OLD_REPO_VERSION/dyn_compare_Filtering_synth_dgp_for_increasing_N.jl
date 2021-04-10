
# The purpose of this script is to compare the filtering capability of the single snapshots
# estimates and the gas filter on a given


#function that Generates the dynamical fitnesses and estimates single snapshot
using StaticNets,JLD, DynNets
 function simWrap1(N,T, NTV;dynType = "SIN", degIOUncMeans =15.*ones(2N),degb = [10, 40] )

    dynFitDgp,indsTVnodes =  StaticNets.dgpDynamic(StaticNets.fooNetModelDirBin1,"SIN",N,T;
                                        NTV = NTV,degIOUncMeans = degIOUncMeans,degb = degb  )
    dynExpDegsDgp = zero(dynFitDgp)
    for t =1:T  dynExpDegsDgp[:,t] = StaticNets.expValStats(StaticNets.fooNetModelDirBin1,dynFitDgp[:,t]) end
    # sanmple dynamical networks from dgp fitnesses
    dynDegsSam = zeros(dynFitDgp)
    for t=1:T
        matSamp = StaticNets.sampl(StaticNets.NetModelDirBin1(dynFitDgp[:,t]),1;  parGroupsIO = dynFitDgp[:,t])
        dynDegsSam[:,t] = [sumSq(matSamp,2);sumSq(matSamp,1)]
    end
    estFitSS =  StaticNets.estimate( StaticNets.SnapSeqNetDirBin1(dynDegsSam); identPost = false,identIter= true )
    return dynFitDgp,dynDegsSam,estFitSS,indsTVnodes
 end

## Simulate estimate and save
Nvals = [ 2000,3000] #[50,60,80,90,100,125,150,200,300,400]# 200##40 #
 minDeg,maxDeg = (10, 40)
 Nsample = 3
 T = 250
    dynFitDgp = Array{Array{Float64,2},3}(length(Nvals),Nsample,2)
    dynDegsSam = Array{Array{Float64,2},3}(length(Nvals),Nsample,2)
    estFitSS = Array{Array{Float64,2},3}(length(Nvals),Nsample,2)
    filFitGas = Array{Array{Float64,2},3}(length(Nvals),Nsample,2)
    indsTVnodes = Array{BitArray{1},3}(length(Nvals),Nsample,2)
    rmseSSandGas = Array{Array{Float64,3},2}(length(Nvals),2)
    storeGasPar = Array{Array{Array{Float64,1},1},3}(length(Nvals),Nsample,2)
    convGasFlag = BitArray{3}(length(Nvals),Nsample,2)
   valDens = [false true]
   for d=1:2
     for indN = 1:length(Nvals)
        N=Nvals[indN]
        NTV = round.(Int,N/2)
        dense = valDens[d] #false#false#true
        if  dense
            unifDeg = round(Int,N*0.1) #15
            degb = [round(Int,N*0.1), round(Int,N*0.8)]
        else
            unifDeg = round(Int,(maxDeg-  minDeg)/2)
             degb = [round(Int,minDeg), round(Int,maxDeg)]#
        end
        tmprmseSSandGas = zeros(2N,2,Nsample)
        prog=0
        @time for s=1:Nsample
            println((d,N,s))
            #sample the dgp
             dynFitDgp[indN,s,d],dynDegsSam[indN,s,d],estFitSS[indN,s,d],indsTVnodes[indN,s,d] =
                   simWrap1(N,T,NTV;degIOUncMeans =unifDeg*ones(2N),degb = degb)
             tmprmseSSandGas[:,1,s] =sqrt.(meanSq((dynFitDgp[indN,s,d] - estFitSS[indN,s,d]).^2,2))
             #estimate single snapshot sequence
              round(s/Nsample,1)>prog? (prog = round(s/Nsample,1);println(prog)):()
              #estimate score_driven_filter_or_dgp
              ## e ora stima il gasPar
              degsIO_T = dynDegsSam[indN,s,d]
              modGasDirBin1 = DynNets.GasNetModelDirBin1(degsIO_T,"FISHER-DIAG")
               parGas,convGasFlag[indN,s,d]  =  DynNets.estimateTarg(modGasDirBin1;SSest = estFitSS[indN,s,d])
               #parGas = [[0.989525],[0.0680809]]
              storeGasPar[indN,s,d] = parGas
              display(parGas[2:3])
              gasFiltFit,~ = DynNets.score_driven_filter_or_dgp(modGasDirBin1,[parGas[1];parGas[2];parGas[3]])
              tmprmseSSandGas[:,2,s] =sqrt.(meanSq((dynFitDgp[indN,s,d] - gasFiltFit).^2,2))
              filFitGas[indN,s,d] = gasFiltFit
              end
              if false # some descriptive plots
                  s=1
               fig, ax = subplots(3,2)
               ax[1,1][:set_title]("N = $(N) const Deg = $(unifDeg)  deg bounds = $(degb)  ")
               ax[1,1][:plot](dynFitDgp[indN,s,d][indsTVnodes[indN,s,d],:]')
               ax[1,2][:plot](dynFitDgp[indN,s,d][.!indsTVnodes[indN,s,d],:]')
               ax[2,1][:plot](estFitSS[indN,s,d][indsTVnodes[indN,s,d],:]',".")
               ax[3,1][:plot](dynFitDgp[indN,s,d][:,:]',estFitSS[indN,s,d][indsTVnodes[indN,s,d],:]',".")
               ax[2,2][:plot](dynDegsSam[indN,s,d][indsTVnodes[indN,s,d],:]')
               ax[3,2][:plot](dynDegsSam[indN,s,d][indmax(maximum(dynDegsSam[indN,s,d],2)),:])
           end
        rmseSSandGas[indN,d] = tmprmseSSandGas
     end
 end

 using JLD2
 save_fold = "./data/estimatesTest/asympTest/"
 @save(save_fold*"FilterTestDynamicIncreasingSize_extraLarge_2_denseAndSparse_$(Nsample).jld",
              dynFitDgp,dynDegsSam,estFitSS,indsTVnodes,storeGasPark,rmseSSandGas,filFitGas,convGasFlag)


##

save_fold = "./data/estimatesTest/asympTest/"
@load(save_fold*"FilterTestDynamicIncreasingSize_extraLarge_2_denseAndSparse_$(Nsample).jld",
             dynFitDgp,dynDegsSam,estFitSS,indsTVnodes,storeGasPark,rmseSSandGas,filFitGas,convGasFlag)

## Test single path sampling and FilterTestDynamicIncreasingSize

T=300
    ## Simulate estimate and save
    minDeg,maxDeg = (10, 160)
    Nsample = 1
    Nvals =  200#[200,300,400]#40 #[40,60,80,90,100,125,150]#
    d=1
    indN = 1
    N= 200
    NTV = round.(Int,N/2)
    dense = false#false#true
    if  dense
        unifDeg = round(Int,N*0.1) #15
        degb = [round(Int,N*0.1), round(Int,N*0.8)]
    else
        unifDeg = round(Int,(maxDeg-  minDeg)/2)
         degb = [round(Int,minDeg), round(Int,maxDeg)]#
    end
    tmprmseSSandGas = zeros(2N,2,Nsample)
    prog=0

    #sample the dgp
    dynFitDgp_test,dynDegsSam_test,estFitSS_test,indsTVnodes_test  =
       simWrap1(N,T,NTV;degIOUncMeans =unifDeg*ones(2N),degb = degb)
    tmprmseSSandGas_test =sqrt.(meanSq((dynFitDgp_test - estFitSS_test).^2,2))
    #estimate single snapshot sequence
    #estimate score_driven_filter_or_dgp
    ## e ora stima il gasPar
    degsIO_T = dynDegsSam_test

 modGasDirBin1 = DynNets.GasNetModelDirBin1(degsIO_T,"")#"")###
 parGas,convGasFlag_test  =  DynNets.estimateTarg(modGasDirBin1;SSest = estFitSS_test)

 modGasDirBin1_2 = DynNets.GasNetModelDirBin1(degsIO_T,"FISHER-DIAG")#"")###
  parGas2,convGasFlag_test  =  DynNets.estimateTarg(modGasDirBin1_2;SSest = estFitSS_test)
 #gasFiltFit,~ = DynNets.score_driven_filter_or_dgp( DynNets.GasNetModelDirBin1(degsIO_T),[parGas[1];parGas[2];0.07])
 println((parGas[2],parGas[3]))
println((parGas2[2],parGas2[3]))
gasFiltFit,~ = DynNets.score_driven_filter_or_dgp( modGasDirBin1,[parGas[1];parGas[2];parGas[3]])
    tmprmseSSandGas_test =sqrt.(meanSq((dynFitDgp_test - gasFiltFit).^2,2))
    filFitGas_test = gasFiltFit
 indsNodesPlot = [1,N]
 #plot(degsIO_T[indsNodesPlot,:]')
 #close("all")
    fig = figure("pyplot_subplot_mixed",figsize=(10,10)) # Create a new blank figure
    subplot(221) # Create the 1st axis of a 2x2 arrax of axes
    title("min max degs =  $((minimum(degsIO_T[indsNodesPlot[1],:]), maximum(degsIO_T[indsNodesPlot[1],:])))") # Give the most recent axis a title
    plot(1:T,dynFitDgp_test[indsNodesPlot[1],:],"--",
       1:T,estFitSS_test[indsNodesPlot[1],:],".")
    plot(1:T,gasFiltFit[indsNodesPlot[1],:],"-")
    suptitle("Compare Fiter Single path 2 nodes")
    subplot(223)
    title("min max degs =  $((minimum(degsIO_T[indsNodesPlot[2],:]), maximum(degsIO_T[indsNodesPlot[2],:])))") # Give the most recent axis a title
    plot(1:T,dynFitDgp_test[indsNodesPlot[2],:],"--",
       1:T,estFitSS_test[indsNodesPlot[2],:],".")
    plot(1:T,gasFiltFit[indsNodesPlot[2],:],"-")

 #gasFiltFit,~ = DynNets.score_driven_filter_or_dgp( DynNets.GasNetModelDirBin1(degsIO_T),[parGas2[1];parGas2[2];0.07])
 gasFiltFit,~ = DynNets.score_driven_filter_or_dgp( modGasDirBin1_2,[parGas2[1];parGas2[2];parGas2[3]])
    tmprmseSSandGas_test =sqrt.(meanSq((dynFitDgp_test - gasFiltFit).^2,2))
    filFitGas_test = gasFiltFit
 indsNodesPlot = [1,N]
 #plot(degsIO_T[indsNodesPlot,:]')
 #close("all")
    fig = figure("pyplot_subplot_mixed",figsize=(10,10)) # Create a new blank figure
    subplot(222) # Create the 1st axis of a 2x2 arrax of axes
    title("min max degs =  $((minimum(degsIO_T[indsNodesPlot[1],:]), maximum(degsIO_T[indsNodesPlot[1],:])))") # Give the most recent axis a title
    plot(1:T,dynFitDgp_test[indsNodesPlot[1],:],"--",
       1:T,estFitSS_test[indsNodesPlot[1],:],".")
    plot(1:T,gasFiltFit[indsNodesPlot[1],:],"-")
    suptitle("Compare Fiter Single path 2 nodes")
    subplot(224)
    title("min max degs =  $((minimum(degsIO_T[indsNodesPlot[2],:]), maximum(degsIO_T[indsNodesPlot[2],:])))") # Give the most recent axis a title
    plot(1:T,dynFitDgp_test[indsNodesPlot[2],:],"--",
       1:T,estFitSS_test[indsNodesPlot[2],:],".")
    plot(1:T,gasFiltFit[indsNodesPlot[2],:],"-")



plot(maximum(degsIO_T,2),".")
plot(mean(degsIO_T,2),".")
plot(minimum(degsIO_T[indsTVnodes_test,:],2))

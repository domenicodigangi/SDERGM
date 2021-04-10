
# The purpose of this script is to compare the filtering capability of the single snapshots
# estimates and the gas filter on a given


#function that Generates the dynamical fitnesses and estimates single snapshot

using Distributed
addprocs(4)

## Simulate estimate and save
@everywhere push!(LOAD_PATH,"../../DynNets/src/")
@everywhere push!(LOAD_PATH,"../../../DynNets/src/")
@everywhere push!(LOAD_PATH,"../DynNets/src/")
@everywhere push!(LOAD_PATH,"./DynNets/src/")


@everywhere using SharedArrays, Statistics, Utilities, DynNets, StaticNets, JLD

@everywhere begin
    function simWrap1(N,T, NTV;dynType = "SIN", degIOUncMeans =0, degb = [10, 40],
                      dynFitDgp = zeros(10,50), indsTVnodes=falses(10))

       if sum(dynFitDgp)==0
           dynFitDgp,indsTVnodes =  StaticNets.dgpDynamic(StaticNets.fooNetModelDirBin1,dynType,N,T;
                                               NTV = NTV,degIOUncMeans = degIOUncMeans,degb = degb  )
       end
       dynExpDegsDgp = zero(dynFitDgp)
       for t =1:T  dynExpDegsDgp[:,t] = StaticNets.expValStats(StaticNets.fooNetModelDirBin1,dynFitDgp[:,t]) end
       # sanmple dynamical networks from dgp fitnesses
       dynDegsSam = zeros(size(dynFitDgp))
       for t=1:T
           matSamp = StaticNets.sampl(StaticNets.NetModelDirBin1(dynFitDgp[:,t]),1;  parGroupsIO = dynFitDgp[:,t])
           dynDegsSam[:,t] = [sumSq(matSamp,2);sumSq(matSamp,1)]
       end
       #println(dynDegsSam[:,4])
       estFitSS =  StaticNets.estimate( StaticNets.SnapSeqNetDirBin1(dynDegsSam); identPost = true,identIter= false,targetErr = 0.001 )
       for t=1:T
           dynFitDgp[:,t] =  StaticNets.identify(StaticNets.NetModelDirBin1(dynDegsSam[:,t]),
                                              dynFitDgp[:,t]; idType = "equalIOsums" )
           estFitSS[:,t] =  StaticNets.identify(StaticNets.NetModelDirBin1(dynDegsSam[:,t]),
                                              estFitSS[:,t]; idType =  "equalIOsums")
       end
       return dynFitDgp,dynDegsSam,estFitSS,indsTVnodes
       end

    N = 10# [500,1000,1500] #[50,60,80,90,100,125,150,200,300,400]# 200##40 #
  minDeg,maxDeg = (3,8)
  degb = [round(Int,minDeg), round(Int,maxDeg)]#
  unifDeg = round(Int,(maxDeg-  minDeg)/2)
  dgpType = "STEPS"
  Nsample = 100
  targetStatic = true
  T=250
    dynDegsSam = SharedArray{Float64,3}(2N,T,Nsample)
    estFitSS = SharedArray{Float64,3}(2N,T,Nsample)
    #dynFitDgp = Array{Array{Float64,2},2}(Nsample,2)
    filFitGas = SharedArray{Float64,3}(2N,T,Nsample)
    rmseSSandGas =  SharedArray{Float64,3}(2N,2,Nsample)
    storeGasPar = SharedArray{Float64,3}(2N,3,Nsample)#  SharedArray{Array{Array{Float64,1},1},2}(Nsample,2)
    # convGasFlag = SharedArrayBitArray{2}(Nsample,2)

     NTV = round.(Int,N/2)

     tmprmseSSandGas = SharedArray{Float64,3}(2N,2,Nsample) #zeros(2N,2,Nsample)
    prog=0
     dynFitDgp,~,~,indsTVnodes =
           simWrap1(N,T,NTV;dynType =dgpType,degIOUncMeans =unifDeg*ones(2N),degb = degb)
end

@sync @distributed for s=1:Nsample
        println((N,s))
        #sample the dgp
        @time ~,degsIO_T,estFitSS_T,~ =
               simWrap1(N,T,NTV;dynType =dgpType,degIOUncMeans =unifDeg*ones(2N),degb = degb, dynFitDgp = dynFitDgp,indsTVnodes = indsTVnodes)
       fooPar =[zeros(2N), zeros(2N) ,zeros(2N)]# [zeros(2N),zeros(1),zeros(1)]#
       indsGroups = [Int.(1:2N),Int.(1:2N)] #[Int.(1:2N),ones(2N)] #
       modGasDirBin1 = DynNets.GasNetModelDirBin1(degsIO_T, fooPar, indsGroups,"FISHER-DIAG")

       rmseSSandGas[:,1,s] =sqrt.(meanSq((dynFitDgp - estFitSS_T).^2,2))
       estFitSS[:,:,s] = estFitSS_T
       if targetStatic
          estFitStatic = DynNets.estimateSnapSeq(modGasDirBin1, degsIO_T = mean(degsIO_T,dims=2) )
          estFitSS_T = repeat(estFitStatic, 1, 5)
        end
       @time        parGas,~  =  DynNets.estimateTarg(modGasDirBin1;SSest =estFitSS_T)


        #estimate single snapshot sequence
        #estimate score_driven_filter_or_dgp
        ## e ora stima il gasPar
        dynDegsSam[:,:,s] = degsIO_T
          storeGasPar[:,1,s] = parGas[1]; storeGasPar[:,2,s] = parGas[2]; storeGasPar[:,3,s] = parGas[3]
          display(parGas[2:3])
          gasFiltFit,~ = DynNets.score_driven_filter_or_dgp(modGasDirBin1,[parGas[1];parGas[2];parGas[3]])
          rmseSSandGas[:,2,s] =sqrt.(meanSq((dynFitDgp - gasFiltFit).^2,2))
          filFitGas[:,:,s] = gasFiltFit

      end

      save_fold = "./data/sdergm_estimates/"
      @save(save_fold*"FilterTestDynamicFixedSize_N_$(N)_Nsample_$(Nsample)_scaling_DGP_" * dgpType  *
        "targetStatic_$(targetStatic)"* ".jld",
              dynFitDgp,dynDegsSam,estFitSS,indsTVnodes,storeGasPar,rmseSSandGas,filFitGas,degb,unifDeg)


##
if false
  using PyCall; pygui(:qt); using PyPlot; pygui(true)

indsPlot = [ 6]
simsToDrop = dropdims(sum(isnan.(abs.(filFitGas)), dims = (1,2)).>0, dims = (1,2))
  Nplot =indsPlot[1]
    plot(1:T,dynFitDgp[Nplot,:], "k")
    plot(1:T,estFitSS[Nplot,:,.!simsToDrop],".b", alpha=0.5)
    plot(1:T, filFitGas[Nplot,:,.!simsToDrop],".r", alpha = 0.3)


end

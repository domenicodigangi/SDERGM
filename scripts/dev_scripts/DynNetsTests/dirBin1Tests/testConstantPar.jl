

# The purpose of this script is to compare the filtering capability of the single snapshots
# estimates and the gas filter on a given


#


using JLD,Utilities,StaticNets,DynNets, Clustering
using PyCall; pygui(:qt); using PyPlot


################## Compute score autororrelation on syntetic data
Nsample=10
save_fold = "./data/estimatesTest/asympTest/"
@load(save_fold*"FilterTestDynamicIncreasingSize_large_denseAndSparse_$(Nsample)_scaling.jld", dynFitDgp,
        dynDegsSam,estFitSS,indsTVnodes,storeGasPar,rmseSSandGas,filFitGas,convGasFlag)

 Nvals = [round(Int,length(dynDegsSam[i,1,1][:,1])/2) for i=1:size(dynDegsSam)[1]]
 Nsample = size(dynDegsSam)[2]

#
#
# indN,s,d = 1,1,1
#  T = length(dynDegsSam[1,1,1][1,:])
#  N = Nvals[indN]
#  close()
#  groupsInds = [Int.(1:2N),Int.(ones(2N))];groupsInds[2][.!indsTVnodes[indN,s,d]] = 0
#  degsIO_T = dynDegsSam[indN,s,d]
#  modGasDirBin1 = DynNets.GasNetModelDirBin1(degsIO_T,[ zeros(2N) ,0.9  * ones(1), 0.01 * ones(1)],groupsInds,"FISHER-DIAG")
#  statPar,~ = StaticNets.estimate(StaModType( modGasDirBin1); degIO = meanSq(degsIO_T,2) )
#  xOrder = 1:2N #maximum(degsIO_T,2) - minimum(degsIO_T,2)# mean(degsIO_T,2)#
#  # Score autocorrelation for static parameters estimates
#  sIO_T,gradIO_T = gasScoreSeries(modGasDirBin1,repmat(statPar,1,T);obsT = degsIO_T)
#  tmp = [autocor(Float64.(sIO_T[i,:]),[1])[1] for i=1:2N]
#  plot(xOrder[indsTVnodes[indN,s,d]],tmp[indsTVnodes[indN,s,d]],"*r")
#  plot(xOrder[!indsTVnodes[indN,s,d]],tmp[.!indsTVnodes[indN,s,d]],"*b")

# tmp1 = Array{Float64,2}(1,100)
# tmp1[1,:] = tmp
# clustTmp = Clustering.kmeans(tmp1,2)
# indTvClust = .!(clustTmp.assignments .==2)


# Score autocorrelation for all TV parameters
groupsInds = [Int.(1:2N),Int.(ones(2N))];
 #groupsInds[2][indTvNodes] = 0
gasPar,~ = estimateTarg(modGasDirBin1;groupsInds = groupsInds)
tmpPar,~ = score_driven_filter(modGasDirBin1,array2VecGasPar(modGasDirBin1,gasPar),groupsInds=groupsInds)
sIO_T,gradIO_T = gasScoreSeries(modGasDirBin1,tmpPar;obsT = degsIO_T)
    tmp = [autocor(Float64.(sIO_T[i,:]),[1])[1] for i=1:2N]
    plot(tmp[indsTVnodes[indN,s,d]],".r")
    plot(tmp[.!indsTVnodes[indN,s,d]],".b")




groupsInds = [Int.(1:2N),Int.(zeros(2N))];
 groupsInds[2][indTvClust] = 1
gasPar,~ = estimateTarg(modGasDirBin1;groupsInds = groupsInds)
tmpPar,~ = score_driven_filter(modGasDirBin1,array2VecGasPar(modGasDirBin1,gasPar),groupsInds = groupsInds)
sIO_T,gradIO_T = gasScoreSeries(modGasDirBin1,tmpPar;obsT = degsIO_T)
    tmp = [autocor(Float64.(sIO_T[i,:]),[1])[1] for i=1:2N]
    plot(tmp[indsTVnodes[indN,s,d]],".r")
    plot(tmp[.!indsTVnodes[indN,s,d]],".b")






################## Compute score autororrelation on real data
using JLD,Utilities,StaticNets,DynNets, Clustering
using PyCall; pygui(:qt); using PyPlot
halfPeriod = false
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)

Ttrain = 100
matA_T = AeMidWeekly_T
N2,T = size(degsIO_T[:,1:Ttrain]);N = round(Int,N2/2)

# autocorrelation wrt constant parameters model
 close()
  groupsInds = [Int.(1:2N),Int.(ones(2N))];
  modGasDirBin1 = DynNets.GasNetModelDirBin1(degsIO_T[:,1:Ttrain],[ zeros(2N) ,0.9  * ones(1), 0.01 * ones(1)],groupsInds,"FISHER-DIAG")
    statPar,~ = StaticNets.estimate(StaModType(fooGasNetModelDirBin1); degIO = meanSq(degsIO_T[:,1:Ttrain],2) )
    xOrder = 1:2N #maximum(degsIO_T[:,1:Ttrain],2) - minimum(degsIO_T[:,1:Ttrain],2)# mean(degsIO_T[:,1:Ttrain],2)#
    # Score autocorrelation for static parameters estimates
    sIO_T,gradIO_T = gasScoreSeries(modGasDirBin1,repmat(statPar,1,T);obsT = degsIO_T[:,1:Ttrain])
    tmp = [autocor(Float64.(sIO_T[i,:]),[1])[1] for i=1:2N]
    plot(xOrder,tmp,".k")
    indConstDegs = isnan(tmp)# falses(2N)#squeeze(std(degsIO_T[:,1:Ttrain],2).<=0.1,2)
    indNonConstDegs = .!indConstDegs
    plot(xOrder[indConstDegs],zeros(sum(indConstDegs)),"*k")
    grid()
    title("Score Autocorrelation emid, Train sample")



# clusters
maxNgroups = 1
 gasforeFitStore = Array{Array{Float64,2},1}(maxNgroups)
 gasParStore = Array{Array{Array{Float64,1},1},1}(maxNgroups)
 groupsIndsABStore = Array{Array{Float64,1},1}(maxNgroups)
 for Nclusters =1: maxNgroups
 #Nclusters =3
 totClustAssign = ones(Int,2N)
 if Nclusters > 1
 tmp1 = Array{Float64,2}(1,length(tmp) - sum(indConstDegs))
 figure()
    tmp1[1,:] = tmp[.!indConstDegs]
    clustTmp = Clustering.kmeans(tmp1,Nclusters)
    nodesClusterNumber = clustTmp.assignments
    constDegsClustInds = indmin(clustTmp.centers)
    for i=1:Nclusters
        inds = nodesClusterNumber .== i
        totClustAssign[find(indNonConstDegs)[inds]] = i
        plot(xOrder[indNonConstDegs][inds],tmp1[inds],".")
    end
    title("Score Autocorrelation WRT static model \n  divided in $(Nclusters) clusters")
 else
     constDegsClustInds = 1
     plot(xOrder,tmp,".")
     title("Score Autocorrelation WRT static model \n  divided in $(Nclusters) clusters")

 end
 groupsInds = [Int.(1:2N),totClustAssign];
 for i in find(indConstDegs)
     groupsInds[2][i] = constDegsClustInds
 end
 println( unique(groupsInds[2]) )
 groupsIndsABStore[Nclusters] = groupsInds[2]
 plot(xOrder[indConstDegs],zeros(sum(indConstDegs)),"*k")
 grid()

 figure()
 # Score autocorrelation for all TV parameters
 gasPar,~ = estimateTarg(modGasDirBin1;groupsInds = groupsInds)
    tmpPar,~ = score_driven_filter(modGasDirBin1,array2VecGasPar(modGasDirBin1,gasPar),groupsInds=groupsInds)
    sIO_T,gradIO_T = gasScoreSeries(modGasDirBin1,tmpPar;obsT = degsIO_T[:,1:Ttrain])
    tmp = [autocor(Float64.(sIO_T[i,:]),[1])[1] for i=1:2N]
    for i in (unique(groupsInds[2]))
        inds = groupsInds[2].==i
        plot(xOrder[inds],tmp[inds],".")
    end
    grid()
 gasParEstOnTrain = gasPar
    GasforeFit,~ = DynNets.score_driven_filter( DynNets.GasNetModelDirBin1(degsIO_T),
            [gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]];groupsInds=groupsInds)
            gasforeFitStore[Nclusters] = Float64.(GasforeFit)
            gasParStore[Nclusters] = gasPar
 end


# Score autocorrelation for all TV parameters
# groupsInds = [Int.(1:2N),Int.(ones(2N))];
#  groupsInds[2][indConstDegs] = 0
#     gasPar,~ = estimateTarg(modGasDirBin1;groupsInds = groupsInds)
#     tmpPar,~ = score_driven_filter(modGasDirBin1,array2VecGasPar(modGasDirBin1,gasPar),groupsInds=groupsInds)
#     sIO_T,gradIO_T = gasScoreSeries(modGasDirBin1,tmpPar;obsT = degsIO_T[:,1:Ttrain])
#     tmp = [autocor(Float64.(sIO_T[i,:]),[1])[1] for i=1:2N]
#     plot(xOrder[indTvClust],tmp[indTvClust],".r")
#     plot(xOrder[.!indTvClust],tmp[.!indTvClust],".b")
#     gasParEstOnTrain = gasPar
#     GasforeFit,~ = DynNets.score_driven_filter( DynNets.GasNetModelDirBin1(degsIO_T),
#             [gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]];groupsInds=groupsInds)
#             gasforeFit = Float64.(GasforeFit)

figure()
 #select links among largest banks in training sample
 remMat = squeeze(prod(.!matA_T,3),3)#falses(N,N)#
 sizeLeg = 15
 i = 1
 tmpTm1 =  StaticNets.nowCastEvalFitNet(  gasforeFitStore[i] ,matA_T,Ttrain;mat2rem = remMat,shift = 0)
    legTex = ["GAS $(i) group t-1 AUC = $(round(tmpTm1[3],3))"]
    # ROC for predictions from mean over Ttrain estimates
 for i=2:maxNgroups
     tmpTm1 =  StaticNets.nowCastEvalFitNet(  gasforeFitStore[i] ,matA_T,Ttrain;mat2rem = remMat,shift = 0)
        legTex = [legTex ; "GAS $(i) groups t-1 AUC = $(round(tmpTm1[3],3))"]
        legend(legTex)
 end

unique(groupsIndsABStore[4])
using GLM
degsIO_T = repmat(degsIO_T_store,20,1)
 N = round(Int,length(degsIO_T[:,1])/2)
 groupsInds = [Int.(1:2N),Int.(ones(2N))];
 modGasDirBin1 = DynNets.GasNetModelDirBin1(degsIO_T[:,1:Ttrain],[ zeros(2N) ,0.9  * ones(1), 0.01 * ones(1)],groupsInds,"FISHER-DIAG")
 @time estimateTarg(modGasDirBin1)


targetErr=1e-5;identPost=false;identIter=false
      statPar = StaticNets.estimate(StaticNets.SnapSeqNetDirBin1(mean(degsIO_T,2)),targetErr = targetErr ,identPost=identPost,identIter= identIter)
     tmpPar = repmat(statPar,1,T)

s_T,grad_T = gasScoreSeries(modGasDirBin1,tmpPar;obsT = degsIO_T)
pval = zeros(2N)
for i=1:2N
 g_T_t = Float64.(grad_T[i,2:end])
 s_T_tm1 = Float64.(s_T[i,1:end-1])
 g_T_tm1 = Float64.(grad_T[i,1:end-1])
 g_delta_t = Float64.(grad_T[(1:2N).!=i,2:end])
 #I_T_t = Float64.(I_T[i,2:end])
 y = ones(T-1)
 ols = lm([   g_T_t  s_T_tm1.*g_T_t],y)
 yHat = predict(ols)
 yMean =1
 RSS = sum((yHat- yMean).^2)
 ESS = T - RSS
 pval[i] =  ccdf(Chisq(1), ESS)
end
plot(pval)
# aggiungere ai parametri statici anche i parametri statici dei parametri che
# sono gia stati battezzati TV.


figure()
 plot(tmp,pval,".")


#

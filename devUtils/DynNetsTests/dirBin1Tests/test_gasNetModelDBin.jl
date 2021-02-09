using Plots, Utilities, ForwardDiff, AverageShiftedHistograms,JLD
plotly()

Model = fooGasNetModelDirBin1
T = 300
N=10
deltaN = 3
GBA = 1
GW = N # N#N # 20
UMdgp,dgpdegs = linSpacedPar(Model,N; NgroupsW = GW,deltaN = deltaN)
groupsIndAB_I = Utilities.distributeAinVecN(Array{Int64}(1:GBA),GW)
groupsIndAB_O = groupsIndAB_I
splitVec(Vector(1:10))
groupsInds = [Utilities.distributeAinVecN(Array{Int64}(1:GW),N), groupsIndAB_I, groupsIndAB_O ]
Bdgp = 0.9*ones(GBA)
Adgp = 0.1*ones(GBA)
Wdgp = UMdgp.*(1-Bdgp)
scoreScal ="FISHER-DIAG"# ""#"FISHER-EWMA"#
# Test sample
Mod = GasNetModelDirBin1(ones(2,2),[Wdgp , Bdgp, Adgp],groupsInds,scoreScal)
Mod,Y_T,Fitness_T = sampl(Mod,T)


## Load Data
save_fold =   "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
file_nameStart = "Weekly_eMid_Estimates"
file_nameEnd = "_from_$(Dates.format(startDate,"yyyy_mm_dd"))"*
        "_to_" *"$(Dates.format(endDate,"yyyy_mm_dd"))"*".jld"
load_path = save_fold*file_nameStart* "_Bin1_" *file_nameEnd#
@load(save_path,estSS,estGas,estGasTarg,estSS_1GW,estGas_1GW,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)

N = round(Int,length(degsIO_T[:,1])/2)
eMidMod = DynNets.GasNetModelDirBin1(degsIO_T,estGas_1GW[1], [Int.(ones(2N)),Int.(ones(1)),Int.(ones(1))],"")
score_driven_filter(eMidMod,[estGas_1GW[1][1];estGas_1GW[1][2];estGas_1GW[1][3]])
estSS =  DynNets.estimateSnapSeq(eMidMod1Snap)'
## test fisher scoring on data
#time = @elapsed  out = estimate(eMidMod)

println((T,N,GW ,round(time),out[1][2],out[1][3]))
Fitness_T_estTarg, tmp = score_driven_filter( Mod, [out[1][1];out[1][2];out[1][3]])
out
#test the pox of autodiff loglikelihood
function loglikeRE(BAvecRe::Vector)
    B = BAvecRe[1]
    A = BAvecRe[2]
    ~, ll = score_driven_filter( eMidMod, [SSMean;B;A])
    return ll
end
loglikeUN(BAvecUn::Vector) = ( ReB = exp(BAvecUn[1])./(1+exp(BAvecUn[1])); ReA = exp(BAvecUn[2]) ; loglikeRE([ReB;ReA]) )

score = x -> ForwardDiff.gradient(loglikeUN,x)
hess = x -> ForwardDiff.hessian(loglikeUN,x)
x0Re = [0.05;0.5]
    α = 1
    x0Un = [log(x0Re[1]/(1-x0Re[1]));log(x0Re[2])]
    x = x0Un
    llVal = loglikeUN(x0Un)
    maxIt = 1000;it= 0
    minScoreNorm = 10*eps(); scoreNorm = 10
    minObjImp = 10*eps();objImp = 10
    @elapsed while (it<maxIt)& (scoreNorm>minScoreNorm)&(objImp>minObjImp)
        it+=1
        score_x = score(x)
        scoreNorm = sum(score_x.^2)
        x = x - α* hess(x)\score_x
        llValNew = loglikeUN(x)
        objImp = abs(llVal - llValNew)
        llVal = llValNew
        println((it,exp(x[1])/(1+exp(x[1])),exp(x[2]),llValNew ,objImp,scoreNorm))
    end


plot(score_driven_filter( eMidMod, [SSMean;x0Re])[1][:,1+N:end])
##
Fitness_T = Float64.(Fitness_T)
corrMat=  squeeze(StatsBase.crosscor( Fitness_T,Fitness_T,[0]),1)
err = corrMat[isfinite(corrMat)]
tmpAsh = ash(err)
plot(tmpAsh)



## Test Target Estimate

estOutTarg = estimateTargeting(Mod)
Mod_estTarg = GasNetModelBin1(Mod.obsT,estOutTarg[1],groupsInd,scoreScal)

vresParTarg = [Mod_estTarg.Par[1];Mod_estTarg.Par[2];Mod_estTarg.Par[3]]
Fitness_T_estTarg, tmp = score_driven_filter( Mod, vresParTarg)

ppar = plot(Fitness_T_estTarg,title = "$(round.(squeeze(mean(Fitness_T_estTarg,1),1),1))  $(round.(UMdgp,1))")
pdegs = plot(degs_T,title = "$(dgpdegs)  $(round.(squeeze(mean(degs_T,1),1),1))")
plot(ppar,pdegs,layout = (2,1),legend=:none,size = (1200,600))

## Test estimate

estOut = estimate(Mod)
Mod_est = GasNetModelBin1(Mod.obsT,estOut[1],groupsInds,scoreScal)

vresPar = [Mod_est.Par[1];Mod_est.Par[2];Mod_est.Par[3]]
Fitness_T_est, tmp = score_driven_filter( Mod, vresPar)

ppar = plot(Fitness_T_est,title = "$(round.(squeeze(mean(Fitness_T_est,1),1),1))  $(round.(UMdgp,1))")
pdegs = plot(degs_T,title = "$(dgpdegs)  $(round.(squeeze(mean(degs_T,1),1),1))")
plot(ppar,pdegs,layout = (2,1),legend=:none,size = (1200,600))


#

##
using JLD

N_est = 300
Tvals =  [1000]
N_ind = length(Tvals)[1]

n  =1
N = 10 ;
GBAdgp =1; GWdgp = N  ;
GBAest =1; GWest = GWdgp; useStartVal = false

useStartVal? startValFlag = "T" : startValFlag = "F";# false #false #
# save
save_fold = "./data/estimatesTest/"
file_name = "EstTest_$(N)_$(GWdgp)_$(GBAdgp)_$(Tvals)_$(N_est)_$(useStartVal)_Miss_$(GWest)_$(GBAest)_rescaling.jld"
save_path = save_fold*file_name#
estDataT = @load(save_path,groupsIndsEst,groupsIndsDgp,allObs,
        simPar_1,simPar_2,simPar_3,estPar_1,estPar_2,estPar_3,est_conv_flag,est_times)

vestPar = [estPar_1[1,n,:]; estPar_2[n];estPar_3[n]]
Fitness_Tfil ,like= score_driven_filter( Mod,vestPar)


    pparfil = plot(Fitness_Tfil)
    pparfil = plot(Fitness_Tfil)
    plot(ppar,pparfil,layout = (2,1),legend=:none,size = (1200,600))
##
vResGasPar = [Mod.Par[1];Mod.Par[2];Mod.Par[3];]
Fitness_Tfil ,like= score_driven_filter( Mod,vResGasPar)


    pparfil = plot(Fitness_Tfil)
    pparfil = plot(Fitness_Tfil)
    plot(ppar,pparfil,layout = (2,1),legend=:none,size = (1200,600))








## test



















#


# using Plots, AverageShiftedHistograms
# plotly()
#
# T = 500
# N=10
# GBA = 1
# GW = 2N#20
#
# UMdgp = identify!(fooGasNetModelDBin1, [linSpacedFitnesses(N,mode = "UNDIRECTED");linSpacedFitnesses(N;mode = "UNDIRECTED")] )
#
# Bdgp = 0.9*ones(GBA)
# Adgp = 0.2*ones(GBA)
# Wdgp = UMdgp.*(1-Bdgp)
# groupsInd = [distributeAinVecN(Array{Int64}(1:GW),2N), distributeAinVecN(Array{Int64}(1:GBA),GW)]
#
#
#
# ##
#
# Mod= GasNetModelDBin1(ones(2,2),[Wdgp , Bdgp, Adgp],groupsInd)
# Mod,Y_T,Fitness_T = sampl(Mod,T)
# degs_T = Mod.obsT
#
# N = round(Int,length(degs_T[1,:])/2)
#     pparI = plot(Fitness_T[:,1:N])
#     pparO = plot(Fitness_T[:,N+1:2N])
#     pdegsI = plot(degs_T[:,1:N])
#     pdegsO = plot(degs_T[:,N+1:2N]);
#     plot(pparI,pdegsI,pparO,pdegsO,layout = (4,1),legend=:none,size = (1200,600))
#
#
# ##
# err =degs_T[:,1]
#
# rng_space = 0:.001:6
# tmpAsh = ash(err;rng = rng_space)
# plot(tmpAsh)
# ##
# vResGasPar = [Mod.Par[1];Mod.Par[2];Mod.Par[3];]
# Fitness_Tfil ,like= score_driven_filterAndLikeliood( Mod,vResGasPar)
#
#
#     pparIfil = plot(Fitness_Tfil[:,1:N])
#     pparOfil = plot(Fitness_Tfil[:,N+1:2N])
#     plot(pparI,pparIfil,pparO,pparOfil,layout = (4,1),legend=:none,size = (1200,600))
#
#
# ## Test estimate
# estimate(Mod)
#
#
# # Mod2 = GasNetModelDBin1(degs_T[1:2,:])
# # Fitness_T_snap = estSnapSeq(Mod2;degsIOT = repmat(Array{Int,1}(1:10)',2,1))
# #
# # plot(Fitness_T_snap')
#
# sampleDynNetsDBin1(T,"CONST";groupsInds=groupsInd, sinMean = UMdgp)
# ##
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# #

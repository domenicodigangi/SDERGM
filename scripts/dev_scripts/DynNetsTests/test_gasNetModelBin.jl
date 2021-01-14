using Plots, BenchmarkTools
plotly()

Model = fooGasNetModelBin1
T = 300
N=10
GBA = 1
GW = 2#N#N#1# N#N # 20
groupsInds = [Utilities.distributeAinVecN(Array{Int64}(1:GW),N), Utilities.distributeAinVecN(Array{Int64}(1:GBA),GW)]

UMdgp,dgpdegs = linSpacedPar(fooGasNetModelBin1,N; NgroupsW = GW)
Bdgp = 0.9*ones(GBA)
Adgp = 0.1*ones(GBA)
Wdgp = UMdgp.*(1-Bdgp[groupsInds[2]])
scoreScal = "FISHER-EWMA"#""#
Mod = GasNetModelBin1(ones(2,2),[Wdgp , Bdgp, Adgp],groupsInds,scoreScal)
Mod,Y_T,Fitness_T = sampl(Mod,T)
degs_T = squeeze(sum(Y_T,2),2)
#time = @elapsed  out = estimate(Mod)
#out
#println((T,N,GW ,round(time),out[1][2],out[1][3]))




Fitness_T = Float64.(Fitness_T)
corrMat=  squeeze(StatsBase.crosscor( Fitness_T,Fitness_T,[0]),1)
# ppar = plot(Fitness_T,title = "$(round.(squeeze(mean(Fitness_T,1),1),1))  $(round.(UMdgp,1))")
# pdegs = plot(degs_T,title = "$(dgpdegs)  $(round.(squeeze(mean(degs_T,1),1),1))")
# plot(ppar,pdegs,layout = (2,1),legend=:none,size = (1200,600))

#
p1 = contour(corrMat,fill=true,yflip = true, title = "Scaling = " *scoreScal )
    p2 = histogram(corrMat[.!Bool.(eye(corrMat))] )
    plot(p1,p2,layout= (2,1),size=(600,600),legend = :none)
    corrMat=  squeeze(StatsBase.crosscor( degs_T,degs_T,[0]),1)
    p1 = contour(corrMat,fill=true,yflip = true, title = " degs Corr with Scaling = " *scoreScal )
    p2 = histogram(corrMat[.!Bool.(eye(corrMat))] )
    plot(p1,p2,layout= (2,1),size=(600,600),legend = :none)    #

# ppar = plot(Fitness_T,title = "$(round.(squeeze(mean(Fitness_T,1),1),1))  $(round.(UMdgp,1))")
# pdegs = plot(degs_T,title = "$(dgpdegs)  $(round.(squeeze(mean(degs_T,1),1),1))")
#     plot(ppar,pdegs,layout = (2,1),legend=:none,size = (1200,600))
#
#     vresPar = [Mod.Par[1];Mod.Par[2];Mod.Par[3]]
#     Fitness_T_filt, tmp = score_driven_filter( Mod, vresPar)
#     ppar_filt = plot(Fitness_T_filt,title = "$(round.(squeeze(mean(Fitness_T_filt,1),1),1))  $(round.(UMdgp,1))")
#     plot(ppar,ppar_filt ,layout = (2,1),legend=:none,size = (1200,600))


## Test Single snapshot estimates
#@time Fitness_T_estSnap = estimateSnapSeq(Mod)
    # UMsingSnap = squeeze(mean(Fitness_T_estSnap,1),1)
    # UMdgp
    # ppar = plot(Fitness_T_estSnap,title = "$(round.(squeeze(mean(Fitness_T_estSnap,1),1),1))  $(round(UMdgp,1))")
    # pdegs = plot(degs_T,title = "$(dgpdegs)  $(round.(squeeze(mean(degs_T,1),1),1))")
    # plot(ppar,pdegs,layout = (2,1),legend=:none,size = (1200,600))

## Test Target Estimate

estOutSS = estimateSnapSeq(Mod)

#println(Mod.obsT[110,:])

estOutTarg = estimateTargeting(Mod)
Mod_estTarg = GasNetModelBin1(Mod.obsT,estOutTarg[1],groupsInd,scoreScal)

vresParTarg = [Mod_estTarg.Par[1];Mod_estTarg.Par[2];Mod_estTarg.Par[3]]
Fitness_T_estTarg, tmp = score_driven_filter( Mod, vresParTarg)

ppar = plot(Fitness_T_estTarg,title = "$(round.(squeeze(mean(Fitness_T_estTarg,1),1),1))  $(round.(UMdgp,1))")
pdegs = plot(degs_T,title = "$(dgpdegs)  $(round.(squeeze(mean(degs_T,1),1),1))")
plot(ppar,pdegs,layout = (2,1),legend=:none,size = (1200,600))

## Test estimate

estOut = estimate(Mod)
Mod_est = GasNetModelBin1(Mod.obsT,estOut[1],groupsInd,scoreScal)

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








##



















#

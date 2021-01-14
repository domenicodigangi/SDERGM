
using JLD, Distributions, Plots, StatPlots, AverageShiftedHistograms, DynNets
 plotly()

 N_est = 5
 NGWpoints = 3
 Tvals =  [  600 400 200 ]#[50 250 1250]#
 N_Tind = length(Tvals)[1]
 scoreRescType = ""#"FISHER-EWMA"#
 useStartVal = false# true# false
 Nmin = 10
 Nspacing =  10
 Nmax = 70
Nvals = Vector(Nmax:-Nspacing:Nmin)
 N_Nind = length(Nvals)



 save_fold = "./data/estimatesTest/compComplexity/"
 file_name = "EstTest_varN_Nest_$(N_est)_$(Nmin)_$(Nspacing)_$(Nmax)_varGW_$(NGWpoints)_T_$(Tvals)_corectSpec"* scoreRescType *".jld"
 save_path = save_fold*file_name#
 @load(save_path,est_conv_flag,est_times,est_flag,Nvals,GWvalStore ,useStartVal)
 
## histograms of static parameters

function Bias(A::Array{<:Real,2},B::Array{<:Real,2})
    return  Bias =    (A.-B)./std(A .-B)
end
Bias(A::SharedArray{<:Real,2},B::SharedArray{<:Real,2}) = BiasRel( convert(Array{Real},A),convert(Array{Real},B))


pW1 = Plots.Plot{Plots.PlotlyBackend};
pB1 = Plots.Plot{Plots.PlotlyBackend};
pA1 = Plots.Plot{Plots.PlotlyBackend};
ind = 0;
for ind = 1:length(Tvals)
#ind=1
errA =  Bias(simPar_3[ind,:,:],estPar_3[ind,:,:])
errB =   Bias(simPar_2[ind,:,:],estPar_2[ind,:,:])
errW =  (simPar_1[ind,:,:].- estPar_1[ind,:,:])./(repmat(std(estPar_1[ind,:,:],1),N_est,1))

rng_space = -5:.1:5

err= errW
tmpAsh = ash(err;rng =rng_space )
pW1 = plot(tmpAsh,
    title = scoreRescType *"  T= $(Tvals[ind]),
    W     ")
plot!(Normal(0,1),linewidth = 5)

err= errB
tmpAsh = ash(err;rng = rng_space)
pB1 = plot(tmpAsh,
title =    " $(N_est) Est   B   ")
plot!(Normal(0,1),linewidth = 5)

err= errA
rng_space = -3:.1:3
tmpAsh = ash(err;rng = rng_space)
pA1 = plot(tmpAsh,
title = "StartVal:" * startValFlag * "  GW: $(GWest), GBA: $(GBAest), N: $(N),   A   ")
plot!(Normal(0,1),linewidth = 5)
exp = "pl_T_$(ind)=  plot(pW1,pB1,pA1,layout=(3,1) )"
eval(parse(exp))
end


exp2 = " " ; for i = 1:ind exp2 *=  "pl_T_$(i) , " end

exp3 = "plot( " * exp2 * " layout = (1,ind),size=(1350,600),legend =:none)"
eval(parse(exp3))
#plot(pl_T_1 , pl_T_2,pl_T_3 , pl_T_3, layout = (1,4),size=(1350,600),legend =:none)

##

RavgMSE(A,B) = sqrt.(mean((A .- B).^2,1))
AvgBias(A,B) = mean((A .- B),1)
tmpModelEst = DynNets.GasNetModelBin1(zeros(10,10),[ones(GWest),ones(GBAest),ones(GBAest)],groupsIndsEst,scoreRescType)
tmpModelDgp = DynNets.GasNetModelBin1(zeros(10,10),[ones(GWdgp),ones(GBAdgp),ones(GBAdgp)],groupsIndsDgp,scoreRescType)

rmseTotGas = zeros(N_est,N,N_ind)
biasTotGas = zeros(N_est,N,N_ind)
rmseTotSnap = zeros(N_est,N,N_ind)
biasTotSnap = zeros(N_est,N,N_ind)
dgpTvParT = Array{Array{Float64,3},1}(N_ind)
filTvParT = Array{Array{Float64,3},1}(N_ind)
ParSnapT =  Array{Array{Float64,3},1}(N_ind)
# Compute MSE Bias and Snap estimates
@time for indT = 1:N_ind
    #indT=1
            T =Tvals[indT]
            tmpParDgp = zeros(T,N,N_est)
            tmpParFil = zeros(T,N,N_est)
            tmpParSnap = zeros(T,N,N_est)
        for n=1:N_est
    #n=100
    obsT = allObs[:,1:T,indT,n]'
    vecParEst =  [estPar_1[indT,n,:]; estPar_2[indT,n];  estPar_3[indT,n]]
    vecParDgp =  [simPar_1[indT,n,:]; simPar_2[indT,n];  simPar_3[indT,n]]
    GroupIndsEst = Array{Array{<:Int,1},1}(2)
    GroupIndsEst[1] = groupsIndsEst[1]; GroupIndsEst[2] = groupsIndsEst[2]
    tmpParFil[:,:,n] =   DynNets.score_driven_filterAndLikeliood(tmpModelEst,vecParEst,
                                                degsT = obsT,
                                                groupsInds = GroupIndsEst )[1]

    fitSnap = DynNets.estSnapSeq(DynNets.GasNetModelBin1(obsT)) # sequence of snapshots estimate
    # remove the infinites
    for i=1:5
        infInd = find(.!isfinite.(fitSnap))

        replInd = infInd[infInd .> 2]
        fitSnap[replInd] = fitSnap[replInd.-1]
        infInd = find(.!isfinite.(fitSnap))
        replInd = infInd[infInd .< (Int(N*T) -2) ]
        fitSnap[replInd] = fitSnap[replInd.+1]
    end
    tmpParSnap[:,:,n] = fitSnap


    GroupIndsDgp = Array{Array{<:Int,1},1}(2)
    GroupIndsDgp[1] = groupsIndsDgp[1]; GroupIndsDgp[2] = groupsIndsDgp[2]
    tmpParDgp[:,:,n] = DynNets.score_driven_filterAndLikeliood(tmpModelDgp,vecParDgp,
                                                degsT = obsT,
                                                groupsInds = GroupIndsDgp )[1]


            rmseTotGas[n,:,indT] = RavgMSE(tmpParFil[:,:,n],tmpParDgp[:,:,n] )
            biasTotGas[n,:,indT] = AvgBias(tmpParFil[:,:,n],tmpParDgp[:,:,n] )
            rmseTotSnap[n,:,indT] = RavgMSE(fitSnap,tmpParDgp[:,:,n] )
            biasTotSnap[n,:,indT] = AvgBias(fitSnap,tmpParDgp[:,:,n] )
        end
    dgpTvParT[indT]  = tmpParDgp
    filTvParT[indT]  = tmpParFil
    ParSnapT[indT] = tmpParSnap
end

## Example single realization plots
Nlines=10
indPl = round(Int,linspace(1,N,Nlines))
C(g::ColorGradient) = RGB[g[z] for z=linspace(0,1,Nlines)]
lineColors = C(cgrad(:rainbow))

indT = 2
n=4
T = Tvals[indT]

uncMeansFit = simPar_1[indT,n,:]./(1-simPar_2[indT,n,:])
fit_T = [uncMeansFit'; dgpTvParT[indT][:,:,n] ]
#Compute the skewness of each degree's distribution
probMat =  1./(1+ Base.exp(uncMeansFit .+ uncMeansFit' ));DynNets.putZeroDiag!(probMat );
    degVar = sum(probMat.*(1-probMat),2);
    degSkew = round(degVar.^(-3/2).*sum( probMat.*(1-probMat).*(1-2.*probMat),2),3)
uncMeansFilFit = estPar_1[indT,n,groupsIndsEst[1]]./(1-estPar_2[indT,n,:])
fil_fit_T = [uncMeansFilFit'; filTvParT[indT][:,:,n] ]
fitSnap = [zeros(1,N); ParSnapT[indT][:,:,n]]
 p=plot()
for i = 1:Nlines p = plot!(fit_T[:,indPl[i]],  color = lineColors[i] , lab = "$(degSkew[i])" ) end #lineColors
for i = 1:Nlines p = plot!(fil_fit_T[:,indPl[i]],linestyle=:dash,  color = lineColors[i]  ) end #lineColors
for i = 1:Nlines p= plot!(fitSnap[:,indPl[i]],linestyle=:dot,markershape=:circle,
    markersize=1,linewidth = 0, markerstrokecolor=  lineColors[i], color = lineColors[i]  ) end

p =plot(p,titlefont = font(16),xlims = (0,50) ,titleloc=:center,
    title = "Gas DGP,N $(N) T = $(T) B = $(simPar_2[indT,n,1]), A = $(simPar_3[indT,n,1])",size=(1350,600))
##
p = plot()
 for i = 1:Nlines p = plot!(ObsAllT[indT][:,indPl[i]],  color = lineColors[i] ) end #lineColors
 plot(p)


## Box Plots
yLim = (0,1)
    pRT1 = boxplot( rmseTotGas[:,:,1]  , ylims = yLim  , legend = false, title = "T = $(Tvals[1]), N = $(N), GWdgp =  $(GWdgp), GWest = $(GWest)")
    pRT2 = boxplot( rmseTotGas[:,:,2]  , ylims = yLim  , legend = false, title = "T = $(Tvals[2])")
      yLim = (-0.8,1.1)
    pBT1 = boxplot( biasTotGas[:,:,1] , ylims = yLim , legend = false, title = "T = $(Tvals[1])")
    pBT2 = boxplot( biasTotGas[:,:,2] , ylims = yLim , legend = false, title = "T = $(Tvals[2])")
     ptot = plot(pRT1,pRT2, pBT1,pBT2, titlefont = font(12) ,titleloc=:center,size=(1350,600) )
##
yLim = (0,1.5)
    pRT1 = boxplot( rmseTotSnap[:,:,1]  , ylims = yLim  , legend = false, title = "T = $(Tvals[1]), N = $(N), GWdgp =  $(GWdgp), Single Snapshot ")
    pRT2 = boxplot( rmseTotSnap[:,:,2]  , ylims = yLim  , legend = false, title = "T = $(Tvals[2])")
          yLim = (-1.4,0.8)
    pBT1 = boxplot( biasTotSnap[:,:,1] , ylims = yLim , legend = false, title = "T = $(Tvals[1])")
    pBT2 = boxplot( biasTotSnap[:,:,2] , ylims = yLim , legend = false, title = "T = $(Tvals[2])")
     ptot = plot(pRT1,pRT2 ,pBT1,pBT2, titlefont = font(12) ,titleloc=:center,size=(1350,600) )


##






##

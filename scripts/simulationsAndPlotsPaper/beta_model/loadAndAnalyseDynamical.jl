

# The purpose of this script is to compare the filtering capability of the single snapshots
# estimates and the gas filter on a given


#function that Generates the dynamical fitnesses and estimates single snapshot
using DynNets,JLD
 using PyCall; pygui(:qt); using PyPlot
 Nsample_2 = 3
 save_fold = "./data/estimatesTest/asympTest/"
 @load(save_fold*"FilterTestDynamicIncreasingSize_extraLarge_denseAndSparse_$(Nsample_2)_scaling.jld", dynFitDgp,
        dynDegsSam,estFitSS,indsTVnodes,storeGasPar,rmseSSandGas,filFitGas,convGasFlag)


 dynDegsSam_2 = dynDegsSam
 rmseSSandGas_2 = rmseSSandGas
 indsTVnodes_2 = indsTVnodes
 Nvals_2 = [round(Int,length(dynDegsSam[i,1,1][:,1])/2) for i=1:size(dynDegsSam)[1]]

Nsample= 10
 @load(save_fold*"FilterTestDynamicIncreasingSize_large_denseAndSparse_$(Nsample)_scaling.jld", dynFitDgp,
       dynDegsSam,estFitSS,indsTVnodes,storeGasPar,rmseSSandGas,filFitGas,convGasFlag)


 Nvals = [round(Int,length(dynDegsSam[i,1,1][:,1])/2) for i=1:size(dynDegsSam)[1]]
 Nsample = size(dynDegsSam)[2]
 T = length(dynDegsSam[1,1,1][1,:])
 Nvals_tot = [Nvals;Nvals_2]

# TWO PLOTS FOR THE PAPER

close("all")
 figure(figsize=(9,7))
 lineW = 5
 ax = axes()
 yPlot= [mean( [sum(dynDegsSam[n,s,1]/(T*Nvals[n]*(Nvals[n]-1))) for s=1:Nsample]) for n=1:length(Nvals)]
 append!(yPlot, [mean( [sum(dynDegsSam_2[n,s,1]/(T*Nvals_2[n]*(Nvals_2[n]-1))) for s=1:Nsample_2]) for n=1:length(Nvals_2)])
 plot(Nvals_tot,yPlot,".-k" ,linewidth=lineW)
 yPlot =  [mean( [sum(dynDegsSam[n,s,2]/(T*Nvals[n]*(Nvals[n]-1))) for s=1:Nsample]) for n=1:length(Nvals)]
 append!(yPlot,[mean( [sum(dynDegsSam_2[n,s,2]/(T*Nvals_2[n]*(Nvals_2[n]-1))) for s=1:Nsample_2]) for n=1:length(Nvals_2)])
 plot(Nvals_tot,yPlot, ".--k",linewidth=lineW)
 xlabel("Number of Nodes ",fontsize=38)
 ylabel("Density",fontsize=38)
 grid()
 ax[:tick_params]("both", labelsize=18)
 legend(["Sparse"; "Dense"],fontsize = 28)
 tight_layout()


close("all")
 figure(figsize=(9,7))
 d=1
 ax = axes()
 plot(Nvals_tot,append!([mean(rmseSSandGas[i,d][:,1,:]) for i=1:length(Nvals)],[mean(rmseSSandGas_2[i,d][:,1,:]) for i=1:length(Nvals_2)]),".-b",linewidth=lineW)
 plot(Nvals_tot,append!([mean(rmseSSandGas[i,d][:,2,:]) for i=1:length(Nvals)],[mean(rmseSSandGas_2[i,d][:,2,:]) for i=1:length(Nvals_2)]),".-r",linewidth=lineW)
 d=2
 ax = axes()
 plot(Nvals_tot,append!([mean(rmseSSandGas[i,d][:,1,:]) for i=1:length(Nvals)],[mean(rmseSSandGas_2[i,d][:,1,:]) for i=1:length(Nvals_2)]),".--b",linewidth=lineW)
 plot(Nvals_tot,append!([mean(rmseSSandGas[i,d][:,2,:]) for i=1:length(Nvals)],[mean(rmseSSandGas_2[i,d][:,2,:]) for i=1:length(Nvals_2)]),".--r",linewidth=lineW)

 #the commented lines plot the rmse only on the TV nodes
 #plot(Nvals_tot,append!([mean(rmseSSandGas[i,d][!indsTVnodes[i,1],1,:]) for i=1:length(Nvals)],[mean(rmseSSandGas_2[i,d][!indsTVnodes_2[i,1],1,:]) for i=1:length(Nvals_2)]),".--b")
 #plot(Nvals_tot,append!([mean(rmseSSandGas[i,d][!indsTVnodes[i,1],2,:]) for i=1:length(Nvals)],[mean(rmseSSandGas_2[i,d][!indsTVnodes_2[i,1],2,:]) for i=1:length(Nvals_2)]),".--r")
 xlabel("Number of Nodes",fontsize=38)
 ylabel("Average RMSE",fontsize=38)
 grid()
 ax[:tick_params]("both", labelsize=18)
 legend(["Beta Model Sparse";"SD-Beta Model Sparse";"Beta Model Dense";"SD-Beta Model Dense"],fontsize = 28)
 tight_layout()


#THE PLOTS FOR THE PAER FINISH HERE

for i=1:3
 println(storeGasPar[1,i,:][1][3])
end




indN,s,d = 1,1,1
 N = Nvals[indN]

 indN = 1
 indsPlot = [1,Nvals[indN]-1]
T =250
close("all")
 ax = axes()
 plot(dynFitDgp[indN,s,d][indsTVnodes[indN,d],:][indsPlot,:]')
 plot(estFitSS[indN,s,d][indsTVnodes[indN,d],:][indsPlot,:]',".")
 ylabel(L"\overleftarrow{\theta}",fontsize = 38)
 xlabel("Time",fontsize=38)
 ylim([-4.5,1.5])
 ylimvals = ylim()
 xlim([0,T])
 grid()
 ax[:tick_params]("both", labelsize=18)
 title("Single Snapshots N = $(Nvals[indN])",fontsize = 38)

close("all")
   ax = axes()
    plot(dynFitDgp[indN,s,d][indsTVnodes[indN,d],:][indsPlot,:]')
   plot(filFitGas[indN,s,d][indsTVnodes[indN,d],:][indsPlot,:]',".")
   ylabel(L"\overleftarrow{\theta}",fontsize = 38)
   xlabel("Time",fontsize=38)
   xlim([0,T])
   ylim(ylimvals)
   grid()
   ylimvals = ylim()
   ax[:tick_params]("both", labelsize=18)
   title("Gas N = $(Nvals[indN])",fontsize = 38)

close("all")
 ax = axes()
 plot(dynFitDgp[indN,s,d][1:N,:][indsTVnodes[indN,d][1:N],:]')
 ylabel(L"\overleftarrow{\theta}",fontsize = 38)
 xlabel("Time",fontsize=38)
 grid()
 xlim([0,T])
  ylim(ylimvals)
 ax[:tick_params]("both", labelsize=18)


close("all")
 ax = axes()
 plot(dynFitDgp[indN,s,d][indsTVnodes[indN,d],:][[1,17],:]')
 ylabel(L"\overleftarrow{\theta}",fontsize = 38)
 xlabel("Time",fontsize=38)
 xlim([0,T])
 ylim(ylimvals)
 grid()
 ylimvals = ylim()
 ax[:tick_params]("both", labelsize=18)






fig, ax = subplots(2,2)
 for f=1:2
     for d=1:2
         ax[1,d][:plot](Nvals,[mean([sum(dynDegsSam[n,s,d]/(T*Nvals[n]*(Nvals[n]-1))) for s=1:Nsample]) for n=1:length(Nvals)] )
         ax[1,d][:set_ylim]([0,1])
         ax[2,d][:plot](Nvals,[mean(rmseSSandGas[i,d][indsTVnodes[i,d],f,:]) for i=1:length(Nvals)])
         ax[2,d][:plot](Nvals,[mean(rmseSSandGas[i,d][.!indsTVnodes[i,d],f,:]) for i=1:length(Nvals)])
          ax[2,d][:set_ylim]([0.1,0.5])
     end
     for i=1:length(ax[:,1]),j=1:length(ax[1,:]) ax[i,j][:grid]() end
 end
ax[2,1][:plot](estFitSS[indsTVnodes,:]',".")


storeGasPar[5,:,:]
rmseSSandGas[i,d][indsTVnodes[i,d],1,:]
rmseSSandGas[i,d][indsTVnodes[i,d],2,:]

indN,s,d = 2,2,2
 parGas,~ =  DynNets.estimateTarg(modGasDirBin1;SSest = estFitSS[indN,s,d])
 degsIO_T = dynDegsSam[indN,s,d]
 modGasDirBin1 = DynNets.GasNetModelDirBin1(degsIO_T)
 parGas,~ =  DynNets.estimateTarg(modGasDirBin1;SSest = estFitSS[indN,s,d])
 gasFiltFit,~ = DynNets.score_driven_filter( DynNets.GasNetModelDirBin1(degsIO_T),[parGas[1];parGas[2];0.05])
 sqrt.(meanSq((dynFitDgp[indN,s,d] - gasFiltFit).^2,2))

plot(gasFiltFit')
plot(dynFitDgp[indN,s,d]')
plot(dynDegsSam[indN,s,d]')
close("all")

for i=[1,20]
    println(i)
end
#
degsDgp = zeros(dynFitDgp)
 for t = 1:T degsDgp[:,t] =  expValStats(fooNetModelDirBin1,dynFitDgp[:,t] ) end
 indPlot = indsTVnodes# [1,round(Int,N/2),N]
 fig, ax = subplots(2,1)
 ax[1,1][:plot](dynFitDgp[indPlot,:]')
 ax[2,1][:plot](degsDgp[indPlot,:]')
sampledNets = sampl(fooSnapSeqNetDirBin1,dynFitDgp;Nsample= Nsample)
sampledDegsIO = zeros(Nsample,2N,T)
estSSFit = zeros(Nsample,2N,T)
# compute the degrees and estimate single snapshot
for n=1:Nsample
    sampledDegsIO[n,:,:] = [sumSq(sampledNets[n,:,:,:],2);sumSq(sampledNets[n,:,:,:],1)]
    estSSFit[n,:,:] = estimate(fooSnapSeqNetDirBin1;degsIO_T= sampledDegsIO[n,:,:] )
end

indPlot = [1,round(Int,N/2),N]
 fig, ax = subplots(2,1)
 ax[1,1][:plot](estSSFit[1,indPlot,:]')
 ax[2,1][:plot](sampledDegsIO[1,indPlot,:]')





##















 ##

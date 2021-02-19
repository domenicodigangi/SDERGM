

using JLD2, StaticNets,  Utilities
using PyCall; pygui(:qt); using PyPlot; pygui(true)

 N = 10#
 dgpType = "SIN"
 Nsample = 100
 targetStatic = true
 T=250
 save_fold = "./data/sdergm_estimates/"
 @load(save_fold*"FilterTestDynamicFixedSize_N_$(N)_Nsample_$(Nsample)_scaling_DGP_" * dgpType  *
      "targetStatic_$(targetStatic)"* ".jld",
      dynFitDgp,dynDegsSam,estFitSS,indsTVnodes,storeGasPar,rmseSSandGas,filFitGas, degb, unifDeg)

T = size(dynDegsSam.s)[2]
  for s=1:Nsample
    for t=1:T
          estFitSS.s[:,t,s] =  StaticNets.identify(StaticNets.NetModelDirBin1(dynDegsSam.s[:,t,s]),
                                                estFitSS.s[:,t,s];idType = "firstZero" )
          filFitGas.s[:,t,s] =  StaticNets.identify(StaticNets.NetModelDirBin1(dynDegsSam.s[:,t,s]),
                                                filFitGas.s[:,t,s];idType = "firstZero" )
        end
      rmseSSandGas.s[:,1,s] = sqrt.(meanSq((dynFitDgp - estFitSS.s[:,:,s]).^2,2))
      rmseSSandGas.s[:,2,s] = sqrt.(meanSq((dynFitDgp - filFitGas.s[:,:,s]).^2,2))
    end
    sum(estFitSS.s[1,:,:])==0 ? () : error()
    # AvgRMSE_SS_SD = (mean(rmseSSandGas.s[:,1,:]) ,   mean(rmseSSandGas.s[:,2,:]))
    #
    #  @save(save_fold*"FilterTestDynamicFixedSize_N_$(N)_Nsample_$(Nsample)_scaling_DGP_" * dgpType  * ".jld",
    #              dynFitDgp,dynDegsSam.s,estFitSS.s,indsTVnodes,storeGasPar,rmseSSandGas.s,filFitGas.s, degb,unifDeg)
    #
    #



 close()
 figure(figsize=(9,7))
 labSize = 26
 tickSize = 20
 legSize = 20
 for n=1:Nsample
    parInd =8
    #subplot(2,1,1);
    plot(1:T,filFitGas.s[ parInd,:,n],"r")
                   plot(1:T,estFitSS.s[parInd,:,n],".b",markersize = 2)
                   plot(1:T,dynFitDgp[parInd,:],"k",linewidth=5)

    #parInd = 17
    # subplot(2,1,2);plot(1:T,filFitGas.s[ parInd,:,n],"r")
    #                plot(1:T,estFitSS.s[parInd,:,n],".b")
    #                plot(1:T,dynFitDgp[parInd,:],"k",linewidth=5)
 end
 legTex = [ "SD- Beta Model "; "Beta Model "; "DGP"]
 legend(legTex,fontsize = legSize,framealpha = 0.95)
 xlabel("Time",fontsize = labSize)
 ylabel("\$ \\overleftarrow{\\theta} \$",fontsize = labSize)

 xticks(fontsize = tickSize)
 yticks(fontsize = tickSize)
 ylim([-2.5,5.5])





#


using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase,CSV, RCall,DynNets,StaticNets
using PyCall; pygui(); using PyPlot

@load("/home/Domenico/Dropbox/Dynamic_Networks/data/congress_covoting_US/juliaEstimates.jld",estParSS_T,obsMat_T,changeStats_T,
        stats_T)

Nterms,T = size(estParSS_T)
for parInd = 1:Nterms
    subplot(Nterms,1,parInd); #plot(1:T,gasFiltPar[parInd,:],"r")
               plot(1:T,estParSS_T[parInd,:],".b")
end

#changeStats_T2 = changeStats_T
t=1
# sortrows(sum(obsMat_T[t],dims = 2),rev=true)
# display(sortrowsUtilities(changeStats_T[t],2))

pvals_SD,~,vConstPar = DynNets.pValStatic_SDERGM(DynNets.GasNetModelDirBinGlobalPseudo(changeStats_T,DynNets.fooGasPar,trues(Nterms),""))
indTvPar = trues(Nterms)#BitArray([true,true])# pvals_SD.< 0.01 #
model = DynNets.GasNetModelDirBinGlobalPseudo(changeStats_T,DynNets.fooGasPar,indTvPar,"")
indTargPar = indTvPar# BitArray([true,true])# indTvPar#
startUM =  meanSq(estParSS_T,2)#zeros(Nterms) #
~,~,ftot_0 = DynNets.pValStatic_SDERGM(DynNets.GasNetModelDirBinGlobalPseudo(changeStats_T[1:5],DynNets.fooGasPar,indTvPar,"") )
#ftot_0 = (meanSq(estParSS_T[:,1:5],2))
estPar,convFlag = DynNets.estimate(model;UM = startUM,indTvPar = indTvPar,indTargPar = indTargPar ,ftot_0 = ftot_0 )#estParSS_T[:,1])
gasParEst = estPar # store gas parameters
gasParVec = zeros(sum(indTvPar)*3); for i=0:(sum(indTvPar)-1) gasParVec[1+3i : 3(i+1)] = estPar[indTvPar][i+1]; end
constParVec = zeros(sum(.!indTvPar)); for i=1:sum(.!indTvPar) constParVec[i] = estPar[.!indTvPar][i][1]; end
vParOptim_0 = gasParVec# [ 0.303867, 0.9, 0.001, 0.1, 0.9, 0.01] #[0.0528737, 0.9, 0.05, -0.0163459, 0.9, 0.05] # [0.0528737, 0.9, 0.0000000005 ] #
 #constParVec = -4.06563 *ones(1)
 gasFiltPar , pseudolike = DynNets.score_driven_filter_or_dgp(model,vParOptim_0,indTvPar;
                              vConstPar = constParVec,ftot_0 = ftot_0)



#static test on single snapshot estimateSlow
constEst,convFlag = DynNets.estimate(model;indTvPar = falses(Nterms),indTargPar = falses(Nterms)  )#estParSS_T[:,1])

@rput(estParSS_T)
@rput(constEst)
@rput(Nterms)
i = 1
mu = constEst[i]
obs = estParSS_T[i,:]
s_mu = mean(obs)

R"
 pVals_R <- array(0,Nterms)
 for (i in 1:Nterms){
  out = t.test(estParSS_T[i,],mu = constEst[[i]],equalVariance = TRUE)
  print(out$p.value)
  pVals_R[i] <- out$p.value
   }
 "
 pVals_SS = @rget(pVals_R)

@show pvals_SD
 @show pVals_SS

 #names = ["Triangles" "Two-stars" "Political Affiliation"]
close()
 figure(figsize=(18,7))
 names = ["Number of Links" "GWESP" "Political Affiliation"]
 labSize = 26
 tickSize = 20
 legSize = 20
 #names = ["Edges" "GWESP" "Political Affiliation"]
 for parInd = 1:Nterms
     subplot(1,Nterms,parInd); plot((1:T).+40,gasFiltPar[parInd,:],"r")
     plot((1:T).+40,estParSS_T[parInd,:],".b")
     plot((1:T).+40,ones(T).*constEst[parInd],"--b")
     #plot(1:T,ones(T).*vConstPar[parInd],"-b")
     title(names[parInd],fontsize = labSize+4)
     ylabel("\$ \\theta_$(parInd) \$",fontsize = labSize)
     xlabel("Congress Number",fontsize = labSize)
     xticks(fontsize = tickSize)
     yticks(fontsize = tickSize)
     legTex = [ "SD-ERGM "; "Cross Sect. ERGM"; "Constant ERGM" ]
     legend(legTex,fontsize = legSize,framealpha = 0.95)
 end

 tight_layout()


##



##--------------------------------------------------------- confidence bounds

#load simulations results divided in chunks
T = 500
      tot_Nsample = 200
      N_each_chunk = 10
      N = 50
      UM = [-3 , -0.5]
      B_Re  = 0.95
      A_Re  = 0.05
      indTvPar = trues(2)
      indTargPar = falses(2)
      Nsample = N_each_chunk
      # Save all Sampled Data
      load_fold = "./data/estimatesTest/sdergmTest/distrib_static_pars/"
      n_chunk = 1
      data = load(load_fold*"sdergm_dgp__$(N)_T_$(T)_Sample_$(Nsample)_UM_$(UM)_intTV_$(Int.(indTvPar[:]))_intTarg_$(Int.(indTargPar[:]))_N_each_chunk_$(N_each_chunk)_chunk_$(n_chunk).jld" )
      obsT_all = data["obsT_all"][1:end-1]
      estPar_all = data["estPar_all"][1:end-1]
      for n_chunk = 2:Int(tot_Nsample/N_each_chunk)
            data = load(load_fold*"sdergm_dgp__$(N)_T_$(T)_Sample_$(Nsample)_UM_$(UM)_intTV_$(Int.(indTvPar[:]))_intTarg_$(Int.(indTargPar[:]))_N_each_chunk_$(N_each_chunk)_chunk_$(n_chunk).jld" )
            global obsT_all = vcat(obsT_all, data["obsT_all"][1:end-1])
            global estPar_all = vcat(estPar_all, data["estPar_all"][1:end-1])
      end

      tot_Nsample = size(estPar_all)[1]
      W_est_all = reduce(vcat, [[estPar_all[i][1][1]  estPar_all[i][2][1]] for i=1:tot_Nsample])
      B_est_all = reduce(vcat, [[estPar_all[i][1][2]  estPar_all[i][2][2]] for i=1:tot_Nsample])
      A_est_all = reduce(vcat, [[estPar_all[i][1][3]  estPar_all[i][2][3]] for i=1:tot_Nsample])
      W_dgp = UM.*(1-B_Re)
      #plt.close()
      fig = figure( figsize=(10,10)) # Create a new blank figure
            indPar = 1
            subplot(321) # Create the 1st axis of a 2x2 arrax of axes
            grid("on") # Create a grid on the axis
            plt.hist(W_est_all[:,indPar])
            legend("W estimates") # Give the most recent axis a titl
            plt.axvline(x=W_dgp[indPar], color=[1,0,0])
            subplot(323) # Create the 1st axis of a 2x2 arrax of axes
            grid("on") # Create a grid on the axis
            plt.hist(B_est_all[:,indPar])
            legend("B estimates") # Give the most recent axis a titl
            plt.axvline(x=B_Re, color=[1,0,0])
            subplot(325) # Create the 4th axis of a 2x2 arrax of axes
            grid("on") # Create a grid on the axis
            plt.hist(A_est_all[:,1])
            plt.axvline(x=A_Re, color=[1,0,0])
            indPar = 2
            subplot(322) # Create the 1st axis of a 2x2 arrax of axes
            grid("on") # Create a grid on the axis
            plt.hist(W_est_all[:,indPar])
            legend("A estimates") # Give the most recent axis a titl
            plt.axvline(x=W_dgp[indPar], color=[1,0,0])
            subplot(324) # Create the 1st axis of a 2x2 arrax of axes
            grid("on") # Create a grid on the axis
            plt.hist(B_est_all[:,indPar])
            plt.axvline(x=B_Re, color=[1,0,0])
            subplot(326) # Create the 4th axis of a 2x2 arrax of axes
            plt.title("A estimates") # Give the most recent axis a title
            grid("on") # Create a grid on the axis
            plt.hist(A_est_all[:,indPar])
            plt.axvline(x=A_Re, color=[1,0,0])
            suptitle("T = $(T), N= $(N), Unc Means = $(UM), targeting = $(indTargPar)")


using Statistics, MLBase, LinearAlgebra, Distributions
# Compute the covariance of different estimates
Nconf_bands = 500
T_obs = size(changeStats_T)[1]
Nsample=Nconf_bands
gasParSampled  = Array{Array{Float64,2}}(undef, Nsample)

par_est_all = [W_est_all[:,1] B_est_all[:,1] A_est_all[:,1] W_est_all[:,2] B_est_all[:,2] A_est_all[:,2]]
cov_hat = cov(par_est_all)
parBootDistr = Distributions.MvNormal(gasParVec,cov_hat)
sampled_gas_par = rand(parBootDistr,Nconf_bands   )
# remove cases with B greater than one
sampled_gas_par = sampled_gas_par[:,((sampled_gas_par[2,:] .>= 1) .+ (sampled_gas_par[5,:] .>= 1)).==0]

N_samp_par = size(sampled_gas_par)[2]

gasFiltPar_conf_all = fill(zeros(N_samp_par,Nterms,T),Nsample)
# Iterate across all SD static pars vectors
 for n=1
     @show(n)
    indTvPar = BitArray([true,true])
    model = DynNets.GasNetModelDirBinGlobalPseudo(changeStats_T,DynNets.fooGasPar,indTvPar,"")

    # Works only when all parameters are time varying
    gasFiltPar_conf = zeros(N_samp_par,Nterms,T_obs)

    for i= 1:N_samp_par
        @show(i)
        samPar = sampled_gas_par[:,i]
        gasParVec = samPar
        tmp, ~ = DynNets.score_driven_filter_or_dgp(model,gasParVec,indTvPar )
            gasFiltPar_conf[i,:,:] = tmp
    end
 end
using Statistics
# Compute confidence bands
quant_vals = [0.95, 0.05]
conf_band = zeros(length(quant_vals),T_obs,2)
 for t=1:T_obs
     for k=1:Nterms
      filt_t = gasFiltPar_conf[:,k,t]
      conf_band[:,t,k] = Statistics.quantile(filt_t[.!isnan.(filt_t)],quant_vals)
    end
 end

# Plot DGP filtered Path and confidence bands
figure()
     legTex = ["DGP";"SDERGM"; "95%"]
     indTvPar = BitArray([true,true])
     #model = DynNets.GasNetModelDirBinGlobalPseudo(changeStats_T,fooGasPar,indTvPar,"")
     estPar = gasParVec # store gas
     parInd = 1
     subplot(1,2,1);
     plot(Int.(1:T_obs).+40, estParSS_T[parInd,:],".b")
     plot(Int.(1:T_obs).+40,ones(T_obs).*constEst[parInd],"--b")
     plot(Int.(1:T_obs) .+ 40,gasFiltPar[parInd,:],"r")
     plt.fill_between(Int.(1:T_obs) .+40, conf_band[1,:,parInd], y2 =conf_band[2,:,parInd],color =(0.9, 0.2 , 0.2, 0.1)  )#, color='b', alpha=.1)
     #N_bands = length(gasFiltPar_conf_all[n])
     # for i=1:N_bands
     #     plot(1:T,gasFiltPar_conf_all[n][i][parInd,:],"r",alpha = 0.02)
     # end
    parInd = 2
    subplot(1,2,2);
    plot((1:T).+40,estParSS_T[parInd,:],".b")
    plot((1:T).+40,ones(T).*constEst[parInd],"--b")
    plot(1:T + 40,gasFiltPar[parInd,:],"r")
    plt.fill_between(1:T, conf_band[1,:,parInd], y2 =conf_band[2,:,parInd],color =(0.9, 0.2 , 0.2, 0.1)  )#, color='b', alpha=.1)
    title(names[parInd],fontsize = labSize+4)
    ylabel("\$ \\theta_$(parInd) \$",fontsize = labSize)
    xlabel("Congress Number",fontsize = labSize)
    xticks(fontsize = tickSize)
    yticks(fontsize = tickSize)
    legTex = [ "SD-ERGM "; "Cross Sect. ERGM"; "Constant ERGM" ]
    legend(legTex,fontsize = legSize,framealpha = 0.95)


# sample sequences of ergms with different parameters' values from R package ergm
# and test the PseudoLikelihoodScoreDrivenERGM filter
using Utilities,AReg,StaticNets,DynNets , JLD,MLBase,StatsBase,CSV, RCall
 using PyCall; pygui(); using PyPlot
 using ErgmRcall

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



# Compute the hessian of the loglikelihood for each sample
using Statistics, ForwardDiff, LinearAlgebra
      Nsample = tot_Nsample
      cov_hat_all = fill( zeros(6,6),Nsample)
      hess_hat_all = fill( zeros(6,6),Nsample)
      for n=1:Nsample
            print(n)
            model = DynNets.GasNetModelDirBinGlobalPseudo(obsT_all[1],DynNets.fooGasPar,indTvPar,"")
            tmpParStat,~ = DynNets.estimate(model;indTvPar = BitArray(undef,2),changeStats_T =obsT_all[n][1:5])
            ftot_0 = [tmpParStat[1][1], tmpParStat[2][1]]
            function likeFun(vecReGasParAll::Array{<:Real,1})
                  oneInADterms  = (maxLargeVal + vecReGasParAll[1])/maxLargeVal
                  foo,loglikelValue = DynNets.score_driven_filter_or_dgp(model,vecReGasParAll,indTvPar;
                                                  obsT = obsT_all[n],ftot_0 = ftot_0 .* oneInADterms)
                  #println(vecReGasPar)
                  return - loglikelValue
            end
            vecAllParGasHat = [W_est_all[n,1], B_est_all[n,1], A_est_all[n,1], W_est_all[n,2], B_est_all[n,2], A_est_all[n,2]]
            likeFun(vecAllParGasHat)
            print( ForwardDiff.gradient(likeFun,vecAllParGasHat))
            hess_hat =  ForwardDiff.hessian(likeFun,vecAllParGasHat)
            cov_hat = round.(pinv(hess_hat),digits = 8)
            cov_hat_all[n] = cov_hat
            hess_hat_all[n] = hess_hat
      end

# Compute the covariance of different estimates
par_est_all = [W_est_all[:,1] B_est_all[:,1] A_est_all[:,1] W_est_all[:,2] B_est_all[:,2] A_est_all[:,2]]
      cov_hat_all2 = cov(par_est_all)
      a =  [cov_hat_all[n][diagind(cov_hat_all[n])] for n =1:Nsample]
      cov_hat_all_diag = permutedims(reshape(hcat(a...), (length(a[1]), length(a))))

      fig, axs = plt.subplots(3,2, figsize=(10,10)) # Create a new blank figure
      for indPar = 1:3
       # Create the 1st axis of a 2x2 arrax of axes
      grid("on") # Create a grid on the axis
      axs[indPar,1].hist(cov_hat_all_diag[cov_hat_all_diag[:,indPar].>0,indPar])
      axs[indPar,1].axvline(x=cov_hat_all2[indPar,indPar], color=[1,0,0])
      grid("on") # Create a grid on the axis
      axs[indPar,2].hist(cov_hat_all_diag[cov_hat_all_diag[:,indPar+3].>0,indPar+3])
      axs[indPar,2].axvline(x=cov_hat_all2[indPar+3,indPar+3], color=[1,0,0])
      end
      suptitle("T = $(T), N= $(N), Unc Means = $(UM), targeting = $(indTargPar)")

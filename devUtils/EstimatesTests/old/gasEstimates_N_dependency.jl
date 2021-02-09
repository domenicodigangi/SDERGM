
# Sample a gas fitness DGP and estimate the parameters for different lenghts of
# time series
#using DynNets
using JLD
using Hwloc
counts = Hwloc.histmap(Hwloc.topology_load())
ncores = counts[:Core]
if nworkers()<ncores
    addprocs(ncores  - nworkers())
end
@everywhere using DynNets
#
N_est = nworkers()
T=200
Nvals = [100  ]#[15,25,35,40,55]
Nmax= maximum(Nvals)
N_N = length(Nvals)
est_times = SharedArray{Float64}(zeros(N_N,N_est))
est_w = SharedArray{Float64}(zeros(N_N,N_est,Nmax))
est_A = SharedArray{Float64}(zeros(N_N,N_est,Nmax))
est_B = SharedArray{Float64}(zeros(N_N,N_est,Nmax))
est_conv_flag = SharedArray{Bool,2}(N_N,N_est)
est_flag = SharedArray{Bool,2}(N_N,N_est)

real_A = SharedArray{Float64}(zeros(N_N,N_est,Nmax))
real_B = SharedArray{Float64}(zeros(N_N,N_est,Nmax))
 @progress "ind" for ind = 1:N_N
    N = Nvals[N_N+1-ind]
    NTV =N
    indTvPars = falses(N);indTvPars[1:NTV] = true

    @sync @parallel for n = 1:N_est
        sim_A = rand((0.02:0.00001:0.8),NTV)
        sim_B = rand((0.4:0.00001:0.999),NTV)
        real_A[ind,n,1:NTV] = sim_A
        real_B[ind,n,1:NTV] = sim_B
        Y_T,fit_T =  DynNets.sampleParDynNetsBin1(T,"Fit-Multi-gas",
                                                    indTvPar = [trues(NTV);falses(N-NTV)],
                                                    α = 5,β=1,
                                                    Agas =sim_A , Bgas = sim_B  )


        degsT = squeeze(sum(Y_T,2),2)

        model_data = DynNets.GasNetModelBin1(degsT,[zeros(i)for i in [N,NTV,NTV]], indTvPars)
        #estimate
        est_times[ind,n] = @elapsed hat_pars,hat_inds,conv_flag = DynNets.estimate(model_data)
        est_w[ind,n,1:N] = hat_pars[1]
        est_A[ind,n,1:NTV] = hat_pars[2]
        est_B[ind,n,1:NTV] = hat_pars[3]
        est_conv_flag[ind,n] = conv_flag
        est_flag[ind,n] = true
        #println([T,n])
    end

    # save
    save_path = "./data/estimatesTest/EstVarN_allTV_T200.jld"
    @save(save_path)
    println([N,sum(est_times[ind,:])])

end

#

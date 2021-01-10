
using AverageShiftedHistograms
using StatPlots
using Plots
using Distributions
plotly()
using JLD


function Bias(A::Array{<:Real,2},B::Array{<:Real,2})
    return  Bias =    (A.-B)./std(A .-B)
end
Bias(A::SharedArray{<:Real,2},B::SharedArray{<:Real,2}) = BiasRel( convert(Array{Real},A),convert(Array{Real},B))

##

N = 10
GBA = 1
GW = 10
N_est = 900

Tvals= [500 1000 5000 10000] #   [100 250 500 1000]#
useStartVal = false
useStartVal? startValFlag = "T" : startValFlag = "F";# false #false #
folderName = "./data/estimatesTest/"
file_name = folderName * "EstTest_$(N)_$(GW)_$(GBA)_$(Tvals)_$(N_est)_$(useStartVal).jld"
estDataT = @load( file_name,simPar_1,simPar_2,simPar_3,estPar_1,estPar_2,estPar_3,est_times)
#

T  = Tvals[1]
ind = find(x->x==T,Tvals)

plot_tot =Vector{Plots.Plot{Plots.PlotlyBackend}}

rng_space = -3:.1:3
pW1 =0
pB1 =0
pA1 =0
for ind = 1:4
#ind = 1

errA = Bias(simPar_3[ind,:,:],estPar_3[ind,:,:])
errB = Bias(simPar_2[ind,:,:],estPar_2[ind,:,:])
errW =  (simPar_1[ind,:,:].- estPar_1[ind,:,:])./(repmat(std(estPar_1[ind,:,:],1),N_est,1))


err= errW
tmpAsh = ash(err;rng =rng_space )
pW1 = plot(tmpAsh,
    title ="  T= $(Tvals[ind]),
    W   mu =  $(round(mean(err[1,:]),2)) ")
plot!(Normal(0,1),linewidth = 5)

err= errB
tmpAsh = ash(err;rng = rng_space)
pB1 = plot(tmpAsh,
title =    " $(N_est) Est   B   mu =  $(round(mean(err[1,:]),2)) ")
plot!(Normal(0,1),linewidth = 5)

err= errA
tmpAsh = ash(err;rng = rng_space)
pA1 = plot(tmpAsh,
title = startValFlag * "  $(GW) / $(GBA) / $(N) A  mu =  $(round(mean(err[1,:]),2)) ")
plot!(Normal(0,1),linewidth = 5)

exp = "pl_T_$(ind)=  plot(pW1,pB1,pA1,layout=(3,1) )"
eval(parse(exp))
end


exp2 = " " ; for i = 1:ind exp2 *=  "pl_T_$(i) , " end

exp3 = "plot( " * exp2 * " layout = (1,ind),size=(1350,600),legend =:none)"
eval(parse(exp3))
#plot(pl_T_1 , pl_T_2,pl_T_3 , pl_T_4, layout = (1,4),size=(1350,600),legend =:none)

##


plot(errW[1,:])




##

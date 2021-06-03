#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\scripts\simulationsAndPlotsPaper\reciprocity_p_star\simulate_AR_and_estimate_variance_edges_GWESP.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\scripts\simulationsAndPlotsPaper\reciprocity_p_star
#Created Date: Wednesday April 28th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Thursday June 3rd 2021 11:28:04 pm
#Modified By:  Domenico Di Gangi
#-----
#Description:
#-----
########


using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"

using ProjUtilities
using Statistics
using DataFrames
using PyPlot
import ScoreDrivenERGM.Utilities.AReg
using ScoreDrivenERGM


@elapsed df = collect_results( datadir("sims", "dgpAR_Fil_SS")) 
df["modelTag"] = string.(df["model"]) 

# @elapsed df = df[typeof.(df.model) .== ScoreDrivenERGM.DynNets.SdErgmPml, :]
@elapsed df = df[ df.T .== 73, :]
# viewtab(df[[:model, :T, :N]])

model = df.model[1]
N = df.N[1]
df.allfVecT_filt[1]
@elapsed begin 

    # df = transform(df, :allfVecT_filt_SS => ByRow(y-> dropdims(mapslices( (x-> AReg.fitARp(Float64.(x), 1)[2] ), y, dims=(2) ), dims=2)) => :std_SS )

    df = transform(df, :allfVecT_filt => ByRow(y-> dropdims(mapslices( (x-> AReg.ar1_variance(x) ), y, dims=(2) ), dims=2)) => :std_SD )

    df = transform(df, :allParDgpT => ByRow(y-> dropdims(mapslices( (x-> AReg.fitARp(Float64.(x), 1)[2] ), y, dims=(2) ), dims=2)) => :std_DGP_est )

    df = transform(df, :dgpSettings => ByRow(x-> getfield(getfield(x,:opt), :sigma)) => :sigma )
end


viewtab(df[[:model, :T, :N, :std_SD]])

begin
df = sort(df, :sigma)
qVals = getfield.(getfield.(df.dgpSettings,:opt), :sigma)
fig, ax = plt.subplots(nrows=2, ncols=length(qVals), figsize=(12, 8))#, sharey = "row")

labNames = ["edges", "gwesp"]
for i=1:nrow(df)
    for p=1:2
        all_data = ([df.std_SD[i][p,:] df.std_DGP_est[i][p,:]] .-df.sigma[i][p])./df.sigma[i][p]
        # plot violin plot
        ax[p, i].boxplot(all_data,
                        showmeans=true, showfliers=false)
        xlims = ax[p, i].get_xlim()
        ax[p, i].hlines(0, xlims[1], xlims[2], linestyle=":" , colors = "r")

        ax[p,i].set_xticklabels([ "SD", "DGP"])
        ax[p,i].grid()
        ax[p,1].set_ylabel(labNames[p])
        ax[1,i].set_title("var dgp = $(qVals[i][2]) ")
    end
end
   
fig.suptitle("Estimate Variance of Latent AR1 in $(model.staticModel.ergmTermsString)  N = $N")
plt.tight_layout()

# plt.savefig("C:\\Users\\digan\\Dropbox\\Dynamic_Networks\\repos\\ScoreDrivenExponentialRandomGraphs\\_research\\analysis_for_paper_revision\\plots_for_revision\\estimate_var_latent_AR1_N_$(N)_$model")
end
#endregion


# test stationarity of filtered parameters
using RCall
R"library('aTSA')"


i = 1
SDFiltPar_T = df.allfVecT_filt[1][2,:,i]
@rput  SDFiltPar_T
R"""
    x = as.numeric(SDFiltPar_T) 
    adf.test(x)
    """
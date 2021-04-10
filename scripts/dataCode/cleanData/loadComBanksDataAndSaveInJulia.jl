using JLD2,CSV
#Load CsV files and save as julia data
T = 59

assetHeader = "names"
for a=2:21
    assetHeader = [assetHeader;"$(a)"]
end

matrT = [ zeros(Int,2,2) for t=1:T]
equityT = [ zeros(Int,1) for t=1:T]
IDsT = [ zeros(Int,1) for t=1:T]
for t=1:T
    matfilename = "/home/Domenico/Dropbox/Network_Ensembles/data/Real_networks_raw_data/bipartite/dataArtAssessing/US_com_banks/matr_$(t).csv"
    matframe = CSV.read(matfilename;header = assetHeader)
    mat = Array{Int,2}(matframe)
    matrT[t] = mat[:,2:end]
    IDsT[t] = mat[:,1] 
    eqfilename = "/home/Domenico/Dropbox/Network_Ensembles/data/Real_networks_raw_data/bipartite/dataArtAssessing/US_com_banks/equities_$(t).csv"
    equframe = CSV.read(eqfilename; header = ["Equity"])
    eqvec = Array{Int,2}(equframe)
    equityT[t] = squeeze(eqvec,2)
end

save_path = "/home/Domenico/Dropbox/Network_Ensembles/data/Real_networks_raw_data/bipartite/dataArtAssessing/USComBanksMatsAndEquities.jld"
@save(save_path,equityT,IDsT,matrT)
##



































#

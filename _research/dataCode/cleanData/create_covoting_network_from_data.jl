
# Load raw congress voting data and create a covoting matrix as described in  Roy,  S.,  Y.  Atchad ́e,  and  G.  Michailidis  (2017).   Change  point  estimation  in  high  dimensionalmarkov  random-field  models.Journal  of  the  Royal  Statistical  Society:   Series  B  (StatisticalMethodology) 79(4), 1187–1206



using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"
DrWatson.greet()

import ScoreDrivenERGM:StaticNets, Utilities

using JLD2,StatsBase,CSV, RCall
using DelimitedFiles

loadPath = datadir("US_congr_covoting", "raw_data") * "/Senate_all_votes.csv"
tmp = readdlm(loadPath,',')
rowNames = tmp[1,:]


ids_members = Int.(tmp[2:end,4])
congNumb = Int.(tmp[2:end,1])
votes = Int.(tmp[2:end,5])
rollNumbAll = Int.(tmp[2:end,3])
[sum(votes.==sort(unique(votes))[i]) for i =1:9]

allCongrs = unique(congNumb)
yesVotesCodes = [1,2,3]
nayVotesCodes = [4,5,6]
Ncongresses = length(allCongrs)
votingMatrices  = Array{Array{Int,2},1}(undef,Ncongresses)
for c=1:Ncongresses
    cong = allCongrs[c]
    inds_c = congNumb.==cong
    votesCong_c = votes[inds_c]
    ids_cong_c = ids_members[inds_c]
    senatorsPres = unique(ids_cong_c)
    N = length(senatorsPres)
    rollNumb = rollNumbAll[inds_c]
    println(c)
    rolls_c = unique(rollNumb)
    Nrolls_c = length(rolls_c)
    #matrix that stores for each roll (row) the vote of each member present(column)
    tmpMat = zeros(Int,Nrolls_c,N)
    for i=1:N
        # select the rows that refer to votes casted by i
        indsId_i = ids_cong_c.==senatorsPres[i]
        votes_i = votesCong_c[indsId_i]
        #indices of votes
        roll_i = rollNumb[indsId_i]
        #esclude multi voting of same member
        length(unique(roll_i)) == length((roll_i)) ? Nrolls_i = length(unique(roll_i)) : error()
        #for each roll of i
        for r=1:Nrolls_i
            rollInd = roll_i[r]
            if  sum(votes_i[r] .== yesVotesCodes)>0
                tmpMat[rollInd,i] = 1
            elseif  sum(votes_i[r] .== nayVotesCodes)>0
                tmpMat[rollInd,i] = 2
            end
        end
        votingMatrices[c] = tmpMat
    end
end

coVotingMatrices  = Array{Array{Float64,2},1}(undef, Ncongresses)
for c = 1:Ncongresses
    println(c)
    mat = votingMatrices[c]
    N = size(mat)[2]
    tmpCovot = zeros(N,N)
    for i=1:N
        vot_i = mat[:,i]
        for j=i+1:N
            vot_j = mat[:,j]
            ind_copresent = (vot_i.!=0).&(vot_j.!=0)
            covot = sum((vot_i.==vot_j).&ind_copresent)
            copresent = sum(  ind_copresent)
            if copresent>0
                tmpCovot[i,j] = covot/copresent
            end
        end
    end
    coVotingMatrices[c] = tmpCovot
end




#identify republicans and democrats. Not used in the paper, hence not checked and not saved
# loadPath = datadir("US_congr_covoting", "raw_data") * "/members_Info.csv"
# tmp = readdlm(loadPath,',')
# rowNames = tmp[1,:]

# ids_members = Int.(tmp[2:end,4])
# congNumb = Int.(tmp[2:end,1])
# partyCode = Int.(tmp[2:end,7])
# Ncong = length(unique(congNumb))
# codesBigParties = zeros(Ncong,2)
# affiliationVecs = Array{Array{String,1},1}(undef, Ncong)
# congressMenPresentVecs = Array{Array{Int,1},1}(undef, Ncong)
# for c = 1:Ncong
#     indsCong_c = congNumb .== c
#     congressMansVec = unique(ids_members[indsCong_c])
#     congressMansPartiesVec = unique(ids_members[indsCong_c])
#     congressMenPresentVecs[c] = congressMansVec
#     Ncongressmas_c = length(congressMansVec)
#     partyCode_c = partyCode[indsCong_c]
#     parties_c  = unique(partyCode_c)
#     Nparties = length(parties_c)
#     countPartiesMembers = zeros(length(parties_c))
#     for i=1:Nparties
#          countPartiesMembers[i] =  sum(partyCode_c.==parties_c[i])

#     end
#     println(length(parties_c))
#     println(countPartiesMembers)
#     if Nparties>2
#         indsBigParties = .!(countPartiesMembers/mean(countPartiesMembers).<0.5)
#         sum(indsBigParties) != 2 ? println("warning more than 2 large parties") : ()#error if more than two parties are flagged as big ones
#         tmp = [countPartiesMembers collect(1:Nparties)]
#         tmp = sortslices(tmp,rev = true, dims=1)
#         indsBigParties = tmp[1:2,2]
#     else
#         indsBigParties = [1 2]
#     end
#     codesBigParties[c,:] = parties_c[Int.(indsBigParties)]
#     tmpVecAffil = ["" for i=1:Ncongressmas_c]
#     for i=1:Ncongressmas_c
#         inds_man_i = ids_members[indsCong_c] .==congressMansVec[i]
#         partyCongressMan_i = unique( partyCode_c[inds_man_i])
#         length(partyCongressMan_i)>1 ? println("warning one man with more parties") : ()
#         partyCongressMan_i=partyCongressMan_i[1]
#         # it does not matter that Ds are democrats and R are reps. it matters only that all Reps have the same letter etc
#         independent = true
#         partyCongressMan_i.== codesBigParties[c,1] ? (tmpVecAffil[i] = "D";independent=false) : ()
#         partyCongressMan_i.== codesBigParties[c,2] ? (tmpVecAffil[i] = "R"; independent = false) : ()
#         independent ? tmpVecAffil[i] = "I" : ()
#     end
#      affiliationVecs[c]=tmpVecAffil
#      congressMenPresentVecs[c] = congressMansVec
# end

#make matrix symmetric and remove inactive nodes
using LinearAlgebra
for c=1:Ncong
        tmp = coVotingMatrices[c]
        indsActive = .!((dropdims(sum(tmp,dims=1), dims=1).==0).&(dropdims(sum(tmp,dims=2), dims=2).==0))
        Nactive = sum(indsActive)
        coVotingMatrices[c] = Symmetric(tmp[indsActive,indsActive])

        # affiliationVecs[c] = affiliationVecs[c][indsActive]
end

# save Data
@save(datadir("US_congr_covoting", "clean_data") * "/coVotingMatrices_julia.jld2", coVotingMatrices,votingMatrices,congressMenPresentVecs)




module Scalings


n_pox_pairs(N) = N*(N-1)/2
export n_pox_pairs


n_pox_dir_links(N) = N*(N-1)
export n_pox_dir_links


"""
Number of links for various network size, in dense regime
"""
denseLScal_DirBin(N, C)  = C * n_pox_dir_links(N)


denseLScal_DirBin(N)  = denseLScal_DirBin(N, 0.2)
export denseLScal_DirBin

"""
Number of links for various network size, in semiDense regime
"""
semiDenseLScal_DirBin(N, C)  = C * n_pox_dir_links(N)/sqrt(N)

semiDenseLScal_DirBin(N)  = semiDenseLScal_DirBin(N, 0.1)
export semiDenseLScal_DirBin


"""
Number of links for various network size, in sparse regime
"""
sparseLScal_DirBin(N, C) =  C * N 

sparseLScal_DirBin(N) = sparseLScal_DirBin(N, 4)
export sparseLScal_DirBin


"""
Average number of reciprocal pairs in Erdos Reny Model 
 - avgL is the expected number of links, related to the probability p of a single link being present by expL = (N^2-N) * p  

"""
erdosRenyRecScal_DirBin(expL, N) = expL^2/(2*(N^2-N)) 
export erdosRenyRecScal_DirBin



end

using PhyloNetworks

# Convert binary vectors to binary strings, then parse binary strings as
# integers for readability (at least in terms of identifying common clusters
# across networks)
function clust2int(clust::Vector{Int64})
    parse(Int,foldl((x,y) -> x*string(y),clust;init="");base=2)
end

# Applies `clust2int` above to the output of `hardwiredClusters`
function clust2int(clust::Matrix{Int64})
    res = zeros(Int64,size(clust,1))
    # keep track of edge numbers for plotting
    edgenum = zeros(Int64,size(clust,1))
    count = 1
    for i = 1:size(clust,1)
        binvec = clust[i,2:(end-1)]
        if sum(binvec) > 1 # remove singleton clusters
            res[count] = clust2int(binvec)
            edgenum[count] = clust[i,1]
            count += 1
        end
    end
    (res[1:(count-1)],edgenum[1:(count-1)])
end

###############################################################################

# Computes RF distance between two networks with the same tips, based on their
# output from `hardwiredClusters`
function rfDist(clust1::Matrix{Int64},clust2::Matrix{Int64})
    # Count the no. of unique clusters represented in one but not the other
    length(symdiff(unique(clust2int(clust1)[1]),unique(clust2int(clust2)[1])))
end

# Applies `rfDist` above to two ::HybridNetwork values with the same tips.
function rfDist(net1::HybridNetwork,net2::HybridNetwork)
    taxon_labels = tipLabels(net1);
    rfDist(hardwiredClusters(net1,taxon_labels),
           hardwiredClusters(net2,taxon_labels))
end

###############################################################################

# Print hardwired clusters that are found in one network but not the other
# Supply (1) the integer-coded cluster, (2) the no. of taxa in the network, and
# (3) the taxon order for representing hardwired clusters
function printCluster(val::Int64,ntaxa::Int64,taxon_labels::Vector{String})
    # Convert integer to bitstring, to binary vector
    binvec = map(x -> parse(Int,x),
                 split(bitstring(val)[(end-ntaxa+1):end],""))
    # Extract taxa corresponding to 1s
    taxon_labels[findall(!iszero,binvec)] |> print; print("\n")
end
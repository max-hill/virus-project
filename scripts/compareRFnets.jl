using PhyloNetworks

l2r = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick");
r2l = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick");

# l2r, 7-reticulation net
l2r_7ret = l2r[8];
# r2l, 6-reticulation net
r2l_6ret = r2l[7];

# Standardize taxa order for representing hardwired clusters
taxon_labels = tipLabels(l2r_7ret);
# ["B589", "SP1777", "K22", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33", "216_II", "BHV5"]

# Matrices describing clusters, which are uniquely represented as binary vectors
l2r_7ret_hcmat = hardwiredClusters(l2r_7ret,taxon_labels);
r2l_6ret_hcmat = hardwiredClusters(r2l_6ret,taxon_labels);

# Convert binary vectors to binary strings, then parse binary strings as
# integers for readability (at least in terms of identifying common clusters
# across networks)
function clust2int(clust::Vector{Int64})
    parse(Int,foldl((x,y) -> x*string(y),
        clust[2:(end-1)];init="");base=2)
end

# Applies `clust2int` above to the output of `hardwiredClusters`
function clust2int(clust::Matrix{Int64})
    map(i -> clust2int(clust[i,:]),1:size(clust,2))
end

clust2int(l2r_7ret_hcmat) 
# [510, 510, 508, 508, 444, 256, 256, 256, 188, 60, 60]
clust2int(r2l_6ret_hcmat)
# [510, 510, 510, 446, 318, 62, 62, 60, 60, 56, 56]

# Computes RF distance between two networks with the same tips, based on their
# output from `hardwiredClusters`
function rfDist(clust1::Matrix{Int64},clust2::Matrix{Int64})
    # Count the no. of unique clusters represented in one but not the other
    length(symdiff(unique(clust2int(clust1)),unique(clust2int(clust2))))
end

# Applies `rfDist` above to two ::HybridNetwork values with the same tips.
function rfDist(net1::HybridNetwork,net2::HybridNetwork)
    taxon_labels = tipLabels(net1);
    rfDist(hardwiredClusters(net1,taxon_labels),
           hardwiredClusters(net2,taxon_labels))
end

rfDist(l2r_7ret,r2l_6ret) # 8
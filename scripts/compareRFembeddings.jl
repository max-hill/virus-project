using PhyloNetworks
include("RFutilities.jl")

# RF nets for L2R
l2r_7ret = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick");
taxon_labels = tipLabels(l2r_7ret);

# Gene tree embeddings for L2R and R2L
l2r_gte = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick.embeddings.tre");
r2l_gte = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick.embeddings.tre");

# Union of (integer-coded) clusters/splits among all embeddings
l2r_gte_splits = foldl(union,
    map(i -> clust2int(hardwiredClusters(l2r_gte[i],taxon_labels))[1],
        1:length(l2r_gte));init=Int64[]) |> sort;
# [6, 17, 24, 25, 49, 56, 57, 59, 130, 132, 134, 187, 191, 260, 313, 315, 317, 319, 447]
r2l_gte_splits = foldl(union,
    map(i -> clust2int(hardwiredClusters(r2l_gte[i],taxon_labels))[1],
        1:length(r2l_gte));init=Int64[]) |> sort;
# [6, 9, 17, 24, 25, 49, 56, 57, 130, 132, 134, 191, 313, 315, 319, 441, 443, 447]

# Clusters/splits in one of L2R or R2L, but not the other
symdiff(l2r_gte_splits,r2l_gte_splits);
# [59, 187, 260, 317, 9, 441, 443]
# In L2R: [59, 187, 260, 317]
# In R2L: [9, 441, 443]

# Print these clusters
for i in [59, 187, 260, 317, 9, 441, 443]
    printCluster(i,9,taxon_labels)
end
# 59: ["C33", "C46", "Cooper", "SP1777", "Titanium_IBR_MLV_vaccine"]
# 187: ["B589", "C33", "C46", "Cooper", "SP1777", "Titanium_IBR_MLV_vaccine"]
# 260: ["216_II", "K22"]
# 317: ["216_II", "C33", "C46", "Cooper", "K22", "Titanium_IBR_MLV_vaccine"]

# 9: ["Cooper", "Titanium_IBR_MLV_vaccine"]
# 441: ["216_II", "B589", "C33", "C46", "Cooper", "Titanium_IBR_MLV_vaccine"]
# 443: ["216_II", "B589", "C33", "C46", "Cooper", "SP1777", "Titanium_IBR_MLV_vaccine"]

###############################################################################

# Gene trees for L2R and R2L
l2r_gt = readMultiTopology("../data/bhv/genetrees/ogrooted/15genes_blksize10000_set1c.treefile");
r2l_gt = readMultiTopology("../data/bhv/genetrees/ogrooted/15genes_blksize10000_rev_set1c.treefile");

# Union of (integer-coded) clusters/splits among all gene trees
l2r_gt_splits = foldl(union,
    map(i -> clust2int(hardwiredClusters(l2r_gt[i],taxon_labels))[1],
        1:length(l2r_gt));init=Int64[]) |> sort;
# [6, 17, 24, 25, 49, 56, 57, 59, 130, 132, 134, 187, 191, 260, 313, 315, 317, 319, 447]
r2l_gt_splits = foldl(union,
    map(i -> clust2int(hardwiredClusters(r2l_gt[i],taxon_labels))[1],
        1:length(r2l_gt));init=Int64[]) |> sort;
# [6, 9, 17, 24, 25, 49, 56, 57, 130, 132, 134, 191, 313, 315, 319, 441, 443, 447]
using PhyloNetworks
include("RFutilities.jl")

# RF nets for L2R
l2r_7ret = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick")[8];
taxon_labels = tipLabels(l2r_7ret);

# Gene tree embeddings for L2R and R2L
l2r_gte = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick.embeddings.tre");
r2l_gte = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick.embeddings.tre");

# Union of (integer-coded) clusters/splits among all embeddings
l2r_gte_splits = foldl(union,
    map(i -> clust2int(hardwiredClusters(l2r_gte[i],taxon_labels))[1],
        1:length(l2r_gte));init=Int64[]) |> sort;
# [24, 28, 40, 44, 56, 60, 62, 66, 126, 188, 190, 192, 254, 320, 384, 444, 448, 508, 510]
r2l_gte_splits = foldl(union,
    map(i -> clust2int(hardwiredClusters(r2l_gte[i],taxon_labels))[1],
        1:length(r2l_gte));init=Int64[]) |> sort;
# [24, 28, 40, 44, 48, 56, 60, 62, 190, 192, 254, 318, 320, 384, 446, 448, 508, 510]

# Clusters/splits in one of L2R or R2L, but not the other
symdiff(l2r_gte_splits,r2l_gte_splits)
# [66, 126, 188, 444, 48, 318, 446]
# In L2R: [66, 126, 188, 444]
# In R2L: [48, 318, 446]

# Print these clusters
for i in [66, 126, 188, 444, 48, 318, 446]
    printCluster(i,9,taxon_labels)
end
# 66: ["K22", "216_II"]
# 126: ["K22", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33", "216_II"]
# 188: ["SP1777", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33"]
# 444: ["B589", "SP1777", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33"]

# 48: ["Titanium_IBR_MLV_vaccine", "Cooper"]
# 318: ["B589", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33", "216_II"]
# 446: ["B589", "SP1777", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33", "216_II"]

###############################################################################

# Gene trees for L2R and R2L
l2r_gt = readMultiTopology("../data/bhv/genetrees/ogrooted/15genes_blksize10000_set1c.treefile");
r2l_gt = readMultiTopology("../data/bhv/genetrees/ogrooted/15genes_blksize10000_rev_set1c.treefile");

# Union of (integer-coded) clusters/splits among all gene trees
l2r_gt_splits = foldl(union,
    map(i -> clust2int(hardwiredClusters(l2r_gt[i],taxon_labels))[1],
        1:length(l2r_gt));init=Int64[]) |> sort;
# [24, 28, 40, 44, 56, 60, 62, 66, 126, 188, 190, 192, 254, 320, 384, 444, 448, 508, 510]
r2l_gt_splits = foldl(union,
    map(i -> clust2int(hardwiredClusters(r2l_gt[i],taxon_labels))[1],
        1:length(r2l_gt));init=Int64[]) |> sort;
# [24, 28, 40, 44, 48, 56, 60, 62, 190, 192, 254, 318, 320, 384, 446, 448, 508, 510]
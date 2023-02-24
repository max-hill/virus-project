using PhyloNetworks
include("RFutilities.jl")

l2r = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick");
r2l = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick");

# l2r, 7-reticulation net
l2r_7ret = l2r[8];
# r2l, 6-reticulation net
r2l_6ret = r2l[7];

# Standardize taxon order for representing hardwired clusters
taxon_labels = tipLabels(l2r_7ret);
# ["B589", "SP1777", "K22", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33", "216_II", "BHV5"]

# Matrices describing clusters, which are uniquely represented as binary vectors
l2r_7ret_hcmat = hardwiredClusters(l2r_7ret,taxon_labels);
r2l_6ret_hcmat = hardwiredClusters(r2l_6ret,taxon_labels);

(l2r_7ret_hcvec,l2r_7ret_edgenum) = clust2int(l2r_7ret_hcmat) 
# Clusters: [510, 510, 508, 508, 444, 188, 60, 60, 60, 60, 56, 24, 320, 320, 382, 62]
# Edge numbers (wrt l2r_7ret): [36, 34, 15, 13, 8, 7, 6, 28, 26, 24, 22, 21, 12, 11, 33, 31]
(r2l_6ret_hcvec,r2l_6ret_edgenum) = clust2int(r2l_6ret_hcmat)
# Clusters: [510, 510, 510, 446, 318, 62, 62, 60, 60, 56, 56, 56, 24, 192]
# Edge numbers (wrt r2l_6ret): [33, 13, 11, 7, 3, 2, 31, 29, 27, 26, 24, 21, 20, 10]

symdiff(unique(l2r_7ret_hcvec),unique(r2l_6ret_hcvec))
# Clusters not shared: [508, 444, 188, 320, 382, 446, 318, 192]
# Clusters from l2r_7ret: [508, 444, 188, 320, 382]
# Edge numbers (wrt l2r_7ret): [15,13, 8, 7, 12,11 33]
# Clusters from r2l_6ret: [446, 318, 192]
# Edge numbers (wrt r2l_6ret): [7, 3, 10]

rfDist(l2r_7ret,r2l_6ret) # 8

###############################################################################
using PhyloPlots
using RCall

# Thicken edges whose hardwired clusters are found in one network but not the other
R"layout(matrix(1:2,nrow=1,ncol=2,byrow=T))"
l2r_7ret.node[13].name = "Ti";
r2l_6ret.node[12].name = "Ti";
plot(l2r_7ret,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     edgewidth=Dict(map(x -> (x,3),[15,13, 8, 7, 12,11, 33])),showedgenumber=true)
plot(r2l_6ret,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     edgewidth=Dict(map(x -> (x,3),[7,3,10])),showedgenumber=true)

# l2r_7ret
for i in [508,444,188,320,382]
    printCluster(i,9,taxon_labels)
end
# ["B589", "SP1777", "K22", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33"]
# ["B589", "SP1777", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33"]
# ["SP1777", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33"]
# ["B589", "K22"]
# ["B589", "K22", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33", "216_II"]

# r2l_6ret
for i in [446,318,192]
    printCluster(i,9,taxon_labels)
end
# ["B589", "SP1777", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33", "216_II"]
# ["B589", "Titanium_IBR_MLV_vaccine", "Cooper", "C46", "C33", "216_II"]
# ["SP1777", "K22"]
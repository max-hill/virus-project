using PhyloNetworks
include("RFutilities.jl")

############## SOFTWIRED CLUSTERS OF THE INFERRED NETRAX NETWORKS ##############
nr_nets = Vector{HybridNetwork}(undef,15);
for i in 1:15
    runi_str = readlines("../results/bhv/netrax/multistart-titanium_subruns/run_$(i-1)_inferred_network.nw")[1];
    # Pre-process newick string so that hybrid node names have the form #H<name>
    # E.g. #H1 instead of #1
    nr_nets[i] = replace(runi_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
        readTopology;
end

taxon_labels = tipLabels(nr_nets[1]);
# Standardize taxon ordering with the RF-Net analyses
# taxon_labels = tipLabels(readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick"));

# Store softwired clusters for each inferred network
clusters = Vector{Vector{Int64}}(undef,15)
for i in 1:15
    # displayedTrees(nr_nets[i],0.0): lists displayed trees after suppressing degree-2
    # nodes except for the root
    clusters[i] = foldl(union, # union of softwired clusters (encoded as integers)
        map(t -> clust2int(hardwiredClusters(t,taxon_labels))[1],
            displayedTrees(nr_nets[i],0.0));init=Int[]) |> sort
end

# show(stdout,"text/plain",clusters)
# [12, 14, 15, 31, 47, 63, 80, 95, 96, 111, 112, 127, 143, 144, 159, 160, 175, 191, 192, 208, 223, 224, 239, 240]
# [12, 14, 15, 31, 96, 111, 112, 127]
# [12, 14, 15, 31, 96, 111, 127, 144, 159]
# [12, 14, 15, 31, 63, 96, 127, 144, 159, 191]
# [12, 14, 15, 31, 96, 127, 144, 159]
# [12, 14, 15, 31, 96, 111, 112, 127]
# [12, 14, 15, 31, 47, 63, 96, 111, 127, 159, 191, 192, 224]
# [12, 14, 15, 31, 47, 63, 96, 111, 127]
# [12, 14, 15, 31, 79, 95, 96, 111, 127, 143, 159, 207, 223, 239]
# [12, 14, 15, 31, 79, 95, 96, 111, 127, 144, 159, 224, 239]
# [12, 14, 15, 31, 96, 111, 127, 144, 224, 239]
# [12, 14, 15, 31, 96, 143, 144, 159, 191]
# [12, 14, 96, 144, 192, 224, 240, 241]
# [12, 14, 15, 31, 79, 95, 96, 111, 127, 144, 159]
# [12, 14, 15, 31, 96, 111, 127, 144, 224, 239]

# Union of all softwired clusters across the inferred networks
foldl(union,clusters;init=Int[]) |> sort
# [12, 14, 15, 31, 47, 63, 79, 80, 95, 96, 111, 112, 127, 143, 144, 159, 160, 175, 191, 192, 207, 208, 223, 224, 239, 240, 241]
for c in [12, 14, 15, 31, 47, 63, 79, 80, 95, 96, 111, 112, 127, 143, 144, 159, 160, 175, 191, 192, 207, 208, 223, 224, 239, 240, 241]
    printCluster(c,8,taxon_labels)
end
# ["C46", "Cooper"]
# ["C46", "Cooper", "Titanium_IBR_MLV_vaccine"]
# ["C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["SP1777", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["SP1777", "216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["B589", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["B589", "216_II"]
# ["B589", "216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["B589", "SP1777"]
# ["B589", "SP1777", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["B589", "SP1777", "216_II"]
# ["B589", "SP1777", "216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["BHV5", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["BHV5", "216_II"]
# ["BHV5", "216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["BHV5", "SP1777"]
# ["BHV5", "SP1777", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["BHV5", "SP1777", "216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["BHV5", "B589"]
# ["BHV5", "B589", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["BHV5", "B589", "216_II"]
# ["BHV5", "B589", "216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["BHV5", "B589", "SP1777"]
# ["BHV5", "B589", "SP1777", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["BHV5", "B589", "SP1777", "216_II"]
# ["BHV5", "B589", "SP1777", "216_II", "C33"]

# Intersection of all softwired clusters across the inferred networks
foldl(intersect,clusters;init=clusters[1]) # [12, 14, 96]
for c in [12, 14, 96]
    printCluster(c,8,taxon_labels)
end
# ["C46", "Cooper"]
# ["C46", "Cooper", "Titanium_IBR_MLV_vaccine"]
# ["B589", "SP1777"]

########################### COMPARE WITH START TREES ###########################

# Start trees for respective runs
start_trees = Vector{HybridNetwork}(undef,15);
for i in 1:15
    start_trees[i] = readTopology("../results/bhv/netrax/multistart-titanium_subruns/run_$(i-1)_start_network.nw");
end

# Store softwired clusters for each start tree
clusters_start = Vector{Vector{Int64}}(undef,15)
for i in 1:15
    clusters_start[i] = foldl(union, # union of softwired clusters (encoded as integers)
        map(t -> clust2int(hardwiredClusters(t,taxon_labels))[1],
            displayedTrees(start_trees[i],0.0));init=Int[]) |> sort
end
# [12, 14, 15, 31, 96]
# [12, 14, 15, 31, 96]
# [12, 14, 15, 31, 96, 127]
# [10, 14, 15, 31, 96]
# [10, 14, 15, 31, 96]
# [12, 14, 15, 31, 96]
# [12, 13, 15, 31, 96]
# [12, 14, 15, 31, 96]
# [10, 11, 15, 96, 111, 127]
# [12, 14, 15, 31, 96]
# [12, 14, 15, 31, 63]
# [10, 11, 15, 31, 63]
# [12, 14, 15, 31, 96]
# [12, 14, 15, 31, 96]
# [10, 11, 15, 31, 63]

# Union of all softwired clusters across the start trees
foldl(union,clusters_start;init=Int[]) |> sort
# [10, 11, 12, 13, 14, 15, 31, 63, 96, 111, 127]

for c in [10, 11, 12, 13, 14, 15, 31, 63, 96, 111, 127]
    printCluster(c,8,taxon_labels)
end
# ["C46", "Titanium_IBR_MLV_vaccine"]
# ["C46", "Titanium_IBR_MLV_vaccine", "C33"]
# ["C46", "Cooper"]
# ["C46", "Cooper", "C33"]
# ["C46", "Cooper", "Titanium_IBR_MLV_vaccine"]
# ["C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["SP1777", "216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["B589", "SP1777"]
# ["B589", "SP1777", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
# ["B589", "SP1777", "216_II", "C46", "Cooper", "Titanium_IBR_MLV_vaccine", "C33"]
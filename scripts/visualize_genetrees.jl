using PhyloNetworks
include("RFutilities.jl")

# gene tree embeddings
l2r_7ret = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick");
taxon_labels = tipLabels(l2r_7ret);

l2r_gte = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick.embeddings.tre");
r2l_gte = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick.embeddings.tre");

for i = 1:length(l2r_gte)
    clust2int(hardwiredClusters(l2r_gte[i],taxon_labels))[1] |> sort |> print
    print("\n")
end
foldl(union,
    map(i -> clust2int(hardwiredClusters(l2r_gte[i],taxon_labels))[1],1:length(l2r_gte));
    init=[]) |> sort
# [6, 17, 24, 25, 49, 56, 57, 59, 130, 132, 134, 187, 191, 260, 313, 315, 317, 319, 447]

foldl(union,
    map(i -> clust2int(hardwiredClusters(r2l_gte[i],taxon_labels))[1],1:length(r2l_gte));
    init=[]) |> sort
# [6, 9, 17, 24, 25, 49, 56, 57, 130, 132, 134, 191, 313, 315, 319, 441, 443, 447]

symdiff([6, 17, 24, 25, 49, 56, 57, 59, 130, 132, 134, 187, 191, 260, 313, 315, 317, 319, 447],
        [6, 9, 17, 24, 25, 49, 56, 57, 130, 132, 134, 191, 313, 315, 319, 441, 443, 447])
# [59, 187, 260, 317, 9, 441, 443]
# [59, 187, 260, 317]
# [9, 441, 443]

for i in [59, 187, 260, 317, 9, 441, 443]
    printCluster(i,9,taxon_labels)
end
# ["C33", "C46", "Cooper", "SP1777", "Titanium_IBR_MLV_vaccine"]
# ["B589", "C33", "C46", "Cooper", "SP1777", "Titanium_IBR_MLV_vaccine"]
# ["216_II", "K22"]
# ["216_II", "C33", "C46", "Cooper", "K22", "Titanium_IBR_MLV_vaccine"]
# ["Cooper", "Titanium_IBR_MLV_vaccine"]
# ["216_II", "B589", "C33", "C46", "Cooper", "Titanium_IBR_MLV_vaccine"]
# ["216_II", "B589", "C33", "C46", "Cooper", "SP1777", "Titanium_IBR_MLV_vaccine"]

###############################################################################

genetrees = readMultiTopology("rfnet_results/15genes_blksize10000_set1c_rooted.treefile");
ec_genetrees = readMultiTopology("rfnet_results/15genes_blksize10000_set1c_rooted_r10.newick.embeddings.tre");
net = readMultiTopology("rfnet_results/15genes_blksize10000_set1c_rooted_r10.newick")[8]
supertree = readMultiTopology("rfnet_results/15genes_blksize10000_set1c_rooted_r10.newick")[1]

writeTopology(net)

for e in net.edge
    try
        setGamma!(e,0.5)
    catch
        @warn "Cannot change gamma in a tree edge"
    end
end
temp = displayedTrees(net,0.4) # 128 displayed trees

@assert length(genetrees) == length(ec_genetrees) "unequal lengths"
discordant = []
for i in 1:length(genetrees)
    d = hardwiredClusterDistance(genetrees[i],ec_genetrees[i],true)
    if d > 0
        push!(discordant,(i,d))
    end
end
# Any[(5, 1), (7, 1), (9, 1), (11, 3), (12, 1), (14, 1)]
# Any[(9, 1), (11, 2)]

R"layout(matrix(1:2,nrow=1,ncol=2,byrow=T))"
R"layout.show(n=2)"
for i in [11]
    plot(genetrees[i],useedgelength=false,
         minorhybridedgecolor="blue",majorhybridedgecolor="red")
    plot(ec_genetrees[i],useedgelength=false,
         minorhybridedgecolor="blue",majorhybridedgecolor="red")
end


genetrees = readMultiTopology("rfnet_results/15genes_blksize10000_rev_set1c_rooted.treefile");
ec_genetrees = readMultiTopology("rfnet_results/15genes_blksize10000_rev_set1c_rooted_r10.newick.embeddings.tre");
net = readMultiTopology("15genes_blksize10000_rev_set1c_rooted_r10.newick")[8]

@assert length(genetrees) == length(ec_genetrees) "unequal lengths"
discordant = []
for i in 1:length(genetrees)
    d = hardwiredClusterDistance(genetrees[i],ec_genetrees[i],true)
    if d > 0
        push!(discordant,(i,d))
    end
end
# Any[(1, 3), (2, 1), (3, 1), (5, 1), (7, 1), (8, 1), (9, 2), (10, 1), (14, 1), (15, 1)]
# Any[(1, 2), (9, 2)]


gt9 = genetrees[9]
gt9.node[5].name = "216-II"
gt9.node[6].name = "Ti"
ec_gt9 = ec_genetrees[9]
ec_gt9.node[5].name = "216-II"
ec_gt9.node[7].name = "Ti"

gt1 = genetrees[1]
gt1.node[5].name = "216-II"
gt1.node[9].name = "Ti"
ec_gt1 = ec_genetrees[1]
ec_gt1.node[5].name = "216-II"
ec_gt1.node[9].name = "Ti"

R"layout(matrix(1:2,nrow=1,ncol=2,byrow=T))"
R"layout.show(n=2)"
for i in [1]
    plot(genetrees[i],useedgelength=false,
         minorhybridedgecolor="blue",majorhybridedgecolor="red")
    R"mtext"("Gene tree 1")
    plot(ec_genetrees[i],useedgelength=false,
         minorhybridedgecolor="blue",majorhybridedgecolor="red")
    R"mtext"("Error-corrected gene tree 1")
end

R"layout(matrix(1:15,nrow=5,ncol=3,byrow=T))"
R"layout.show(n=15)"

for i in 1:length(genetrees)
    plot(genetrees[i],useedgelength=false)
    R"mtext"("Gene tree $i")
end

net = readMultiTopology("rfnet_results/15genes_blksize10000_set1c_rooted_r10.newick")[8]
rotate!(net,-7); rotate!(net,-11); rotate!(net,-14)
rotate!(net,-17); rotate!(net,-22); rotate!(net,-18)

net.node[6].name = "Ti"
net.node[14].name = "216-II"

plot(net,minorhybridedgecolor="blue",majorhybridedgecolor="red",
     shownodelabel=false,
     nodelabel=DataFrame(nodenumber=[3,4,7,9,10,12,13],
                         label=["H18","H29","H31","H26","H22","H20","H24"]),
     nodecex=0.9,nodelabelcolor="red")

net_rev = readMultiTopology("rfnet_results/15genes_blksize10000_rev_set1c_rooted_r10.newick")[8]
rotate!(net_rev,-4)
rotate!(net_rev,-19)
rotate!(net_rev,-8)

net_rev.node[6].name = "Ti"
net_rev.node[14].name = "216-II"

plot(net_rev,minorhybridedgecolor="blue",majorhybridedgecolor="red",
     shownodelabel=false,
     nodelabel=DataFrame(nodenumber=[3,4,6,9,10,12,13],
                         label=["H24","H22","H28","H30","H20","H18","H26"]),
     nodecex=0.9,nodelabelcolor="red")

R"layout(matrix(1:2,nrow=1,ncol=2,byrow=T))"
R"layout.show(n=2)"
plot(net,minorhybridedgecolor="blue",majorhybridedgecolor="red")
plot(net_rev,minorhybridedgecolor="blue",majorhybridedgecolor="red")


hardwiredClusterDistance(net,net_rev,true)

R"layout(matrix(1:9,nrow=3,ncol=3,byrow=T))"
R"layout.show(n=9)"

for i in 1:length(speciesnetworks)
    plot(speciesnetworks[i],useedgelength=false,
    minorhybridedgecolor="blue",majorhybridedgecolor="red")
    R"mtext"("$(i-1) reticulations")
end

R"layout(matrix(1:4,nrow=2,ncol=2,byrow=F))"
R"layout.show(n=4)"

for i in [1,3,5,7]
    plot(speciesnetworks[i],useedgelength=false,
    minorhybridedgecolor="blue",majorhybridedgecolor="red")
    R"mtext"("$(i-1) reticulations")
end

################################################################################

unrooted_genetrees = readMultiTopology("data/15genes_blksize10000_set1c.treefile");
rooted_genetrees = readMultiTopology("rfnet_results/15genes_blksize10000_set1c_rooted.treefile")

# not clock-like: 2
# slower-rate: 4, 6, 10, 15
R"layout(matrix(1:4,nrow=2,ncol=2,byrow=T))"
for i in [4,6,10,15]
    plot(rooted_genetrees[i],useedgelength=true)
    R"mtext"("gene $i")
end

# faster-rate: 1, 3, 5, 7
# 8, 9, 11, 12
# 13, 14
R"layout(matrix(1:4,nrow=2,ncol=2,byrow=T))"
for i in [1,3,5,7]
    plot(rooted_genetrees[i],useedgelength=true)
    R"mtext"("gene $i")
end
for i in [8,9,11,12]
    plot(rooted_genetrees[i],useedgelength=true)
    R"mtext"("gene $i")
end
for i in [13,14,2]
    plot(rooted_genetrees[i],useedgelength=true)
    R"mtext"("gene $i")
end
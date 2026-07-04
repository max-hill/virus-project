#_______________________________________________________________________________
#
# SNP-bootstrap.jl -- Run bootstrap analysis on the snaq networks that were
# obtained from SNP data, and plot the bootstrap values on the networks
#_______________________________________________________________________________

# INSTRUCTIONS: Run the code in this file line-by-line in a Julia REPL in the
# directory /scrips/snaq/
using DataFrames, CSV, SNaQ, PhyloNetworks, PhyloPlots, RCall
dir="../../analysis/snaq/fromSNPs/bootstrap-results/"
mkpath(dir)

#   [336ed68f] CSV v0.10.16
#   [13f3f980] CairoMakie v0.15.11
# ⌅ [861a8166] Combinatorics v1.0.2
#   [a93c6f00] DataFrames v1.8.2
# ⌃ [31c24e10] Distributions v0.25.125
#   [5790dc55] FourLeafMLE v1.0.1
# ⌃ [14197337] GenericLinearAlgebra v0.3.19
# ⌃ [f213a82b] HomotopyContinuation v2.18.2
# ⌃ [033835bb] JLD2 v0.5.15
#   [b964fa9f] LaTeXStrings v1.4.0
#   [23fbe1c1] Latexify v0.16.10
# ⌃ [f1435218] Oscar v1.3.1
#   [33ad39ac] PhyloNetworks v1.3.1
#   [c0d5b6db] PhyloPlots v2.1.0
#   [b98c9c47] Pipe v1.3.0
#   [91a5bcdd] Plots v1.41.6
#   [49802e3a] ProgressBars v1.5.1
#   [92933f4c] ProgressMeter v1.11.0
#   [274fc56d] PythonPlot v1.0.6
#   [6f49c342] RCall v0.14.13
#   [c2bf7a07] SNaQ v1.1.2
#   [2913bbd2] StatsBase v0.34.11
#   [24249f21] SymPy v2.3.3
# ⌃ [bc8888f7] SymPyPythonCall v0.4.1
# ⌃ [0c5d862f] Symbolics v6.46.0
#   [37e2e46d] LinearAlgebra

#_______________________________________________________________________________
#
# Load utility functions
#_______________________________________________________________________________

function rootaboveoutgroup!(net::HybridNetwork, outgroup)
  good = true
  # direct rooting at the outgroup is possible, or
  # the direct parent of the outgroup is a hybrid
  try
    rootatnode!(net, outgroup)
  catch e1
    isa(e1, PhyloNetworks.RootMismatch) || rethrow(e1)
    directEdges!(net)
    oi = findfirst(n -> n.name == outgroup, net.node)
    pa = PhyloNetworks.getparent(net.node[oi].edge[1])
    for ntrials in eachindex(net.node)
      pa.hybrid && break
      good = false
      oi = findfirst(n -> n === pa, net.node)
      pa = PhyloNetworks.getparent(net.node[oi])
    end
    pae = PhyloNetworks.getparentedge(pa)
    try
      rootonedge!(net, pae)
    catch e2
      isa(e2, PhyloNetworks.RootMismatch) || rethrow(e2)
      directEdges!(net)
      msg = "can't root above $outgroup's 1st ancestral hybrid: the major edge is still below a reticulation"
      throw(PhyloNetworks.RootMismatch(msg)) # e2.msg * msg
    end
  end
  msg = (good ? "" : "can't root above $outgroup directly: its parent isn't a hybrid")
  return(msg)
end

function breakedge!(net::HybridNetwork, lab::String)
    it = findfirst(x -> x.name==lab, net.node)
    isnothing(it) && error("taxon $lab not found in network")
    e = getparentedge(net.node[it])
    PhyloNetworks.breakedge!(e, net)
end


#_______________________________________________________________________________
#
# Generate and save 100 bootstrap networks -- THIS HAS ALREADY BEEN DONE
#_______________________________________________________________________________

# taxonset="set1b" # or run with "set1c"
# startnetwork = readTopology(taxonset=="set1b" ? "((((((C46,Cooper),Titanium_IBR_MLV_vaccine),C33),216_II),(B589,SP1777)),BHV5);" : "((((((C46,Cooper),Titanium_IBR_MLV_vaccine),C33),216_II),((SP1777,B589),K22)),BHV5);")

# COMMENTARY: For the boostrap analysis, we chose as our start network the
# best RFnet tree. For set1b, is from the file
# /osf/rfnet/networks/15genes_blksize10000/15genes_blksize10000_set1b_rooted_r10.newick
# and for set1c, the best RFnet tree is this is from
# /osf/rfnet/networks/15genes_blksize10000/15genes_blksize10000_set1c_rooted_r10.newick

# read in the concordance factor data
# df = CSV.read("../../analysis/snaq/fromSNPs/$(taxonset)/tableCF.csv", DataFrame)

# to redo the bootstrap analysis, uncomment the following and run the following
# code. It will do the bootstrap analysis, saving the output to files to
# gitrepo/analysis/snaq/fromSNPs/bootstrap-reults/bootsnaq-snp-{set1b,set1c}.{log,out}.

# bootnet = bootsnaq(startnetwork, df, hmax=2, nrep=100, seed=34417,  filename="bootsnaq-snp-$(taxonset)")
# "virus-project/analysis/snaq/fromSNPs/bootsnaq-snp-set1c"
# for set in ("set1b", "set1c"), ext in ("log", "out")
#     file = "bootsnaq-snp-$set.$ext"
#     mv(file, joinpath(dir, file))
# end

# For ten runs, total runtime was about 13 minutes per bootstrap replicate (10
# replicates = 2 hours 13 minutes) (running with 12 threads).


#_____________________________________________________________________
#
# Plot bootstrap support network for set1b
#_____________________________________________________________________

taxonset="set1b"

# We will annotate the following network by adding bootstrap supports. The
# network is the best network from SNaQ with h=2 using SNP data for set1b. It
# is located in the file /snaq/fromSNPS/set1b/net2.out
bestNetwork_h2_set1b = readnewick("((C33,(#H9:::0.31830799181252717,(((B589,SP1777):0.2279294082463688,(BHV5)#H12:::0.6700282110815909):2.8143275311865317,(216_II,#H12:::0.329971788918409):0.6282723945427968):1.5424108018558045):1.4473827233513228):1.3189378365443265,(C46,(Cooper)#H9:::0.6816920081874729):3.2407330804204553,Titanium_IBR_MLV_vaccine);") 
rootonedge!(bestNetwork_h2_set1b, 7)

# preview network
plot(bestNetwork_h2_set1b,showedgenumber=true,shownodenumber=true)

# load the bootstrap results
bootnet = readmultinewick(joinpath(dir,"bootsnaq-snp-$(taxonset).out"))

# we need to re-root all the bootstrap networks before calculating bootstrap
# supports. This is because root location affects what counts as a clade.
for i in [1:17;19:20;22:65;67:69;71:100]
    rootaboveoutgroup!(bootnet[i],"BHV5")
    print(i)
    print("\n")
    i=i+1
end

# COMMENTARY: There were four bootstrap networks (18, 21, 66, 70) that we
# couldn't reroot in an obvious way. So we didn't re-root them. You can see
# them here:
# plot(bootnet[18],showgamma=true)
# plot(bootnet[21],showgamma=true)
# plot(bootnet[66],showgamma=true)
# plot(bootnet[70],showgamma=true)




# summarize bootstrap supports for tree and hybrid edges
BSe_tree, major_tree_1b = treeedges_support(bootnet,bestNetwork_h2_set1b);
BSn, BSe, BSc, BSgam, BSedgenum = hybridclades_support(bootnet, bestNetwork_h2_set1b);

# plot(bestNetwork_h2_set1b, edgelabel=BSe[:,[:edge,:BS_hybrid_edge]]);

# make table of labels
tmp = filter(row -> !ismissing(row[:edge]), BSe) # filter rows
select!(tmp, [:edge,:BS_hybrid_edge])            # select 2 columns only
rename!(tmp, :BS_hybrid_edge => :proportion)     # rename those columns, to match names in BSe_tree
rename!(tmp, :edge => :edgeNumber)
tmp = vcat(BSe_tree, tmp)
# preview plot:
# plot(bestNetwork_h2_set1b, edgelabel=tmp, nodelabel=BSn[:, [:hybridnode,:BS_hybrid_samesisters]])

# COMMENTARY: The following two networks show the bootstrap support for donor and reciepient of gene flow
plot(bestNetwork_h2_set1b,
     edgelabel=filter(r->r[:BS_minor_sister]>0,BSn)[!,[:edge,:BS_minor_sister]]); # donor support
plot(bestNetwork_h2_set1b,
     edgelabel=filter(row->row[:BS_hybrid]>0, BSn)[!,[:edge,:BS_hybrid]]); # recipient support

# this is the same as bestNetwork_h2_set1b, but with some cosmetic changes
net_1b=readnewick("((C33,(#H9:::0.31830799181252717,(((B589,SP1777):0.2279294082463688,(BHV5)#H12:::0.6700282110815909):2.8143275311865317,(216_II,#H12:::0.329971788918409):0.6282723945427968):1.5424108018558045):1.4473827233513228):1.3189378365443265,(C46,(Cooper)#H9:::0.6816920081874729):3.2407330804204553,Ti);")
rootonedge!(net_1b, 7)
rotate!(net_1b, -9)
rotate!(net_1b, -3)
rotate!(net_1b, -10)
rotate!(net_1b, -7)
breakedge!(net_1b, "B589") # add degree-2 node: to move parent node to the left

# preview the plot
plot(net_1b, edgelabel=tmp, nodelabel=BSn[:, [:hybridnode,:BS_hybrid_samesisters]], edgelabeladj=[.5,-.2], showgamma=true)

# save the plot
outfile=joinpath(dir,"boostrap-annotated-network-SNaQ-SNP-h2-$(taxonset).pdf")
R"pdf"("$outfile", height=4, width=6);
R"par"(mar=[0,0,0,0]); R"layout"([1])
plot(net_1b, edgelabel=tmp, nodelabel=BSn[:, [:hybridnode,:BS_hybrid_samesisters]], edgelabeladj=[.5,-.2], nodelabeladj=[-.7,-.3], showgamma=true)
R"dev.off"()




#_____________________________________________________________________
#
# Plot bootstrap support network for set1c -- new
#_____________________________________________________________________

taxonset="set1c"

# We will annotate the following network by adding bootstrap supports. The
# network is the best network from SNaQ with h=2 using SNP data for set1c. It
# is located in the file /snaq/fromSNPS/set1c/net2.out

# from the file /snaq/fromSNPS/set1c/net2.out
bestNetwork_h2_set1c = readnewick("(C46,(Titanium_IBR_MLV_vaccine,(C33,(#H10:::0.32326003580298523,((((SP1777,K22):0.038648065141035505,B589):0.28823486443461865,(BHV5)#H11:::0.6864639147113899):2.2944605344476803,(216_II,#H11:::0.3135360852886101):0.5581012926269834):1.6572111385089179):1.3722081419095815):1.4032650601214929):3.4819492825051683,(Cooper)#H10:::0.6767399641970148);")
rootonedge!(bestNetwork_h2_set1c, 11)

# load the bootstrap results
bootnet = readmultinewick(joinpath(dir,"bootsnaq-snp-$(taxonset).out"))

# we need to re-root all the bootstrap networks before calculating bootstrap
# supports. This is because root location affects what counts as a clade.
for i in 2:100
    rootaboveoutgroup!(bootnet[i],"BHV5")
    print(i)
    print("\n")
    i=i+1
end
# COMMENTARY: There was one bootstrap networks (1) that we couldn't reroot in
# an obvious way. So we didn't re-root them. You can see it here:
# plot(bootnet[1],showgamma=true)

# preview network
plot(bestNetwork_h2_set1c,showedgenumber=true,shownodenumber=true)

# summarize bootstrap supports for tree and hybrid edges
BSe_tree, major_tree_1c = treeedges_support(bootnet,bestNetwork_h2_set1c);
BSn, BSe, BSc, BSgam, BSedgenum = hybridclades_support(bootnet, bestNetwork_h2_set1c);

# plot(bestNetwork_h2_set1c, edgelabel=BSe[:,[:edge,:BS_hybrid_edge]]);

# make table of labels
tmp = filter(row -> !ismissing(row[:edge]), BSe) # filter rows
select!(tmp, [:edge,:BS_hybrid_edge])            # select 2 columns only
rename!(tmp, :BS_hybrid_edge => :proportion)     # rename those columns, to match names in BSe_tree
rename!(tmp, :edge => :edgeNumber)
tmp = vcat(BSe_tree, tmp)

# preview plot:
# plot(bestNetwork_h2_set1c, edgelabel=tmp, nodelabel=BSn[:, [:hybridnode,:BS_hybrid_samesisters]])

# # COMMENTARY: The following two networks show the bootstrap support for donor and reciepient of gene flow
plot(bestNetwork_h2_set1c, edgelabel=filter(r->r[:BS_minor_sister]>0,BSn)[!,[:edge,:BS_minor_sister]]); # donor support
plot(bestNetwork_h2_set1c, edgelabel=filter(row->row[:BS_hybrid]>0, BSn)[!,[:edge,:BS_hybrid]]); # recipient support

# this is the same as bestNetwork_h2_set1c, but with some cosmetic changes
net_1c=readnewick("(C46,(Ti,(C33,(#H10:::0.32326003580298523,((((SP1777,K22):0.038648065141035505,B589):0.28823486443461865,(BHV5)#H11:::0.6864639147113899):2.2944605344476803,(216_II,#H11:::0.3135360852886101):0.5581012926269834):1.6572111385089179):1.3722081419095815):1.4032650601214929):3.4819492825051683,(Cooper)#H10:::0.6767399641970148);")
rootonedge!(net_1c, 11)
rotate!(net_1c, -11)
rotate!(net_1c, -4)
rotate!(net_1c, -3)
rotate!(net_1c, -2)
rotate!(net_1c, -8)
breakedge!(net_1c, "K22") # add degree-2 node: to move parent node to the left
plot(net_1c,showedgenumber=true,shownodenumber=true)

# # preview the plot
# plot(net_1c, edgelabel=tmp, nodelabel=BSn[:, [:hybridnode,:BS_hybrid_samesisters]], edgelabeladj=[.5,-.2], showgamma=true)


# save the plot
outfile=joinpath(dir,"boostrap-annotated-network-SNaQ-SNP-h2-$(taxonset).pdf")
R"pdf"("$outfile", height=4, width=6);
R"par"(mar=[0,0,0,0]); R"layout"([1])
plot(net_1c, edgelabel=tmp, nodelabel=BSn[:, [:hybridnode,:BS_hybrid_samesisters]], edgelabeladj=[.5,-.2], nodelabeladj=[-.7,-.3], showgamma=true)
# plot(bestNetwork_h2_set1c, edgelabel=tmp, nodelabel=BSn[:, [:hybridnode,:BS_hybrid_samesisters]], edgelabeladj=[.5,-.2], nodelabeladj=[-.3,-.3], showgamma=true)
R"dev.off"()





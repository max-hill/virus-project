using PhyloNetworks

# Read in network using PhyloNetworks
l2r_7ret = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick")[8];
taxon_labels = tipLabels(l2r_7ret);

# Preserves metadata (e.g. gene embeddings) that PhyloNetworks ignores
net_string = readlines("../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick")[8];

# Tuples of hybrid node names and major hybrid edge numbers
# h.edge[2:end]: because for each internal node, the "first" edge is its child edge
mjedges = map(h -> (h.name,
                    filter(e -> e.isMajor,h.edge[2:end])[1].number),
              l2r_7ret.hybrid);
#  ("H20", 2)
#  ("H22", 3)
#  ("H26", 4)
#  ("H30", 12)
#  ("H24", 17)
#  ("H18", 19)
#  ("H28", 29)
d_mjedges = Dict(mjedges);

mnedges = map(h -> (h.name,
                    filter(e -> e.hybrid && !e.isMajor,h.edge)[1].number),
              l2r_7ret.hybrid);
#  ("H20", 14)
#  ("H22", 10)
#  ("H26", 35)
#  ("H30", 32)
#  ("H24", 27)
#  ("H18", 25)
#  ("H28", 6)
d_mnedges = Dict(mnedges);

# Regex expression to capture hybrid name and suppporting genes
rx = r"#(?<hybrid>H[\d]+)\[&genes={(?<genes>[\d,]+)}\]";

# Search for all matches of `rx` and return an iterator over the matches
hname2genes = eachmatch(rx,net_string) |> collect;

d_gene2hedges = Dict() # Maps genes to hybrid edges used in its embedding
for i = 1:length(hname2genes)
    hname = hname2genes[i].captures[1] # hybrid name
    genes = parse.(Int64,split(hname2genes[i].captures[2],",")) # gene numbers
    if length(genes) < ceil(length(taxon_labels)/2) # if minor edge
        d_hedges = d_mnedges
    else # if major edge
        d_hedges = d_mjedges
    end
    hedge = d_hedges[hname]
    for g in genes
        if haskey(d_gene2hedges,g)
            # Collect the hybrid edges used to embed each gene
            d_gene2hedges[g] = push!(d_gene2hedges[g],hedge)
        else
            d_gene2hedges[g] = [hedge]
        end
    end    
end

###############################################################################

# Example: Suppose we are interested in gene 11
d_gene2hedges[11] # embedded with hybrid edges: [2, 3, 17, 19, 29, 32, 35]

# Make the above hybrid edges major, then plot the major tree
l2r_7ret_copy = deepcopy(l2r_7ret);
setGamma!.(l2r_7ret_copy.edge[[2,3,17,19,29,32,35]],1.0)

using PhyloPlots
using RCall
R"layout(matrix(1:2,nrow=1,ncol=2,byrow=T))"
plot(l2r_7ret) # 7-reticulation RF net
plot(majorTree(l2r_7ret_copy)) # embedding of gene 11 in the 7-reticulation RF net
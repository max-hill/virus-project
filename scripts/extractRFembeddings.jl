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
                    filter(e -> !e.isMajor,h.edge[2:end])[1].number),
              l2r_7ret.hybrid);
#  ("H20", 14)
#  ("H22", 10)
#  ("H26", 35)
#  ("H30", 32)
#  ("H24", 27)
#  ("H18", 25)
#  ("H28", 6)
d_mnedges = Dict(mnedges);

# Regex expression to capture hybrid name and supporting genes
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

# Example: Suppose we are interested in gene 12
d_gene2hedges[12] # embedded with hybrid edges: [2, 3, 12, 17, 25, 29, 35]

# Make the above hybrid edges major, then plot the major tree
l2r_7ret_copy = deepcopy(l2r_7ret);
setGamma!.(l2r_7ret_copy.edge[[2, 3, 12, 17, 25, 29, 35]],1.0)
l2r_ec_g12 = majorTree(l2r_7ret_copy)
rotate!(l2r_ec_g12,-15)
l2r_ec_g12.node[5].name = "Ti"

l2r_ogrooted_g12 = readMultiTopology("../data/bhv/genetrees/ogrooted/15genes_blksize10000_set1c.treefile")[12];
rotate!(l2r_ogrooted_g12,-8)
l2r_ogrooted_g12.node[3].name = "Ti"
writeTopology(l2r_ec_g12,"../results/bhv/rfnet/temp/l2r_ec_g12.newick")
writeTopology(l2r_ogrooted_g12,"../results/bhv/rfnet/temp/l2r_ogrooted_g12.newick")

###############################################################################

# Read in network using PhyloNetworks
r2l_6ret = readMultiTopology("../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick")[7];
taxon_labels = tipLabels(r2l_6ret);

# Preserves metadata (e.g. gene embeddings) that PhyloNetworks ignores
net_string = readlines("../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick")[7];

# Tuples of hybrid node names and major hybrid edge numbers
# h.edge[2:end]: because for each internal node, the "first" edge is its child edge
mjedges = map(h -> (h.name,
                    filter(e -> e.isMajor,h.edge[2:end])[1].number),
              r2l_6ret.hybrid);
#  ("H20", 5)
#  ("H26", 6)
#  ("H24", 17)
#  ("H18", 19)
#  ("H28", 23)
#  ("H22", 32)
d_mjedges = Dict(mjedges);

mnedges = map(h -> (h.name,
                    filter(e -> !e.isMajor,h.edge[2:end])[1].number),
              r2l_6ret.hybrid);
#  ("H20", 12)
#  ("H26", 9)
#  ("H24", 25)
#  ("H18", 22)
#  ("H28", 28)
#  ("H22", 2)
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

# Example: Suppose we are interested in gene 9
d_gene2hedges[9] # embedded with hybrid edges: [2, 5, 6, 17, 19, 23]

# Make the above hybrid edges major, then plot the major tree
r2l_6ret_copy = deepcopy(r2l_6ret);
setGamma!.(r2l_6ret_copy.edge[[2, 5, 6, 17, 19, 23]],1.0)
r2l_ec_g9 = majorTree(r2l_6ret_copy)
rotate!(r2l_ec_g9,-7); rotate!(r2l_ec_g9,-2);
rotate!(r2l_ec_g9,-12); rotate!(r2l_ec_g9,-14);
rotate!(r2l_ec_g9,-17); rotate!(r2l_ec_g9,-18);
r2l_ec_g9.node[8].name = "Ti"

r2l_ogrooted_g9 = readMultiTopology("../data/bhv/genetrees/ogrooted/15genes_blksize10000_rev_set1c.treefile")[9];
rotate!(r2l_ogrooted_g9,-8);
r2l_ogrooted_g9.node[6].name = "Ti"
writeTopology(r2l_ec_g9,"../results/bhv/rfnet/temp/r2l_ec_g9.newick")
writeTopology(r2l_ogrooted_g9,"../results/bhv/rfnet/temp/r2l_ogrooted_g9.newick")
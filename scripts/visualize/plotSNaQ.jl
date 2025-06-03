#= 2025-06
- re-use code from `scripts/snaq/{snaq.jl,snaq-97.jl}`
- plot networks with h=2 only

package versions:
  [336ed68f] CSV v0.10.15
  [a93c6f00] DataFrames v1.7.0
  [33ad39ac] PhyloNetworks v1.1.0
  [c0d5b6db] PhyloPlots v2.0.1
  [6f49c342] RCall v0.14.8
=#

using DataFrames
using RCall
# @rlibrary ggplot2 # ERROR: Normalized names are no longer unique: is_facet, ...
R"library(ggplot2)"

# pseudo-likelihood profile versus h
dat = DataFrame(
  score = # data from snaq.jl:
    # from SNPs: lines 209-211, for set1b (8 taxa) and set1c (9 taxa)
    [501.7797654686652, 295.145876646058, 233.87695612845857, 233.87695612806655, # b
     660.7009851999568, 400.83826191436407, 315.80346871955686, 312.09545236068, # c
     # from 15 gene trees, 10,000-bp L2R partition: lines 119-120
     365.1587301188808, 153.94391099485654, 77.30057399899863, 77.3005739986875, # b
     513.0663771591117, 312.7311634534902, 187.5802997215012, 187.5802997214814], # c
  taxonset = repeat(["- K22","+ K22","- K22","+ K22"], inner=4),
  input = repeat(["SNPs", "15 gene trees"], inner=8), # L2R
  h = repeat([0,1,2,3], outer=4)
)
R"""ggplot(data=$dat, aes(x=h, y=score, color=taxonset, linetype=input, shape=input)) +
  geom_line() + geom_point() +
  xlab("h: maximum number of reticulations") +
  ylab("negative pseudolikelihood score") +
  labs(color = "taxon set", linetype = "input from", shape = "input from") +
  guides(color = guide_legend(override.aes=list(shape=NA))) +
  theme_minimal()
""";
R"ggsave"("../../figures/snaq_scores.pdf", width=5, height=3.5, units="in");

# SNaQ networks with h=2

using PhyloNetworks, PhyloPlots
function breakedge!(net::HybridNetwork, lab::String)
    it = findfirst(x -> x.name==lab, net.node)
    isnothing(it) && error("taxon $lab not found in network")
    e = getparentedge(net.node[it])
    PhyloNetworks.breakedge!(e, net)
end
function renameTitatium!(net::HybridNetwork)
    it = findfirst(x -> startswith(x.name, "Titanium"), net.node)
    if isnothing(it)
        @error("Titanium not found in network")
        return nothing
    end
    net.node[it].name = "Ti"
    return net
end

outdir_root = "../../analysis/snaq"
n_8tax_SNP_file = joinpath(outdir_root,"fromSNPs", "set1b", "snaq_net_0123.out");
n_9tax_SNP_file = joinpath(outdir_root,"fromSNPs", "set1c", "snaq_net_0123.out");
n_8tax_15gL2R_file = joinpath(outdir_root,"fromIQtrees", "set1b", "snaq_net_0123.out");
n_9tax_15gL2R_file = joinpath(outdir_root,"fromIQtrees", "set1c", "snaq_net_0123.out");
n_9tax_97gL2R_file = joinpath(outdir_root,"fromIQtrees", "set1c-97genes-L2R", "snaq_net_0123.out");
n_9tax_97gR2L_file = joinpath(outdir_root,"fromIQtrees", "set1c-97genes-R2L", "snaq_net_0123.out");

n_8tax_SNP = readmultinewick(n_8tax_SNP_file)[3] # 3rd network is for h=2
n_9tax_SNP = readmultinewick(n_9tax_SNP_file)[3]
n_8tax_15gL2R = readmultinewick(n_8tax_15gL2R_file)[3]
n_9tax_15gL2R = readmultinewick(n_9tax_15gL2R_file)[3]
n_9tax_97gL2R = readmultinewick(n_9tax_97gL2R_file)[3]
n_9tax_97gR2L_h1, n_9tax_97gR2L = readmultinewick(n_9tax_97gR2L_file)[[2,3]]
# both incorrect: cannot be rooted with BHV5
rootonedge!(n_9tax_97gR2L, 19) # to make BHV1.1 closer to being a clade

nets = [n_8tax_SNP, n_9tax_SNP, n_8tax_15gL2R, n_9tax_15gL2R,
    n_9tax_97gL2R, n_9tax_97gR2L]
renameTitatium!.(nets)

cfinput = ["SNPs","SNPs","15 gene trees","15 gene trees",
           "97 gene trees","97 gene trees"]
details = ["- K22", "+ K22", "- K22", "+ K22", "(L2R)", "(R2L)"]

# root with BHV5, if possible
rootonedge!(nets[1], 7)
rootonedge!(nets[2], 11)
rootonedge!(nets[3], 9)
rootatnode!(nets[4], "BHV5")
rootonedge!(nets[5], 8)
# reverse ismajor above BHV5: to show ML tree has a displayed tree
nets[5].edge[ 8].ismajor = true
nets[5].edge[11].ismajor = false

# rotate edges to untangle them when plotting
for i in [-9,-3,-10] rotate!(nets[1],i); end
breakedge!(nets[1], "B589") # add degree-2 node: to move parent node to the left
for i in [-11,-4,-3,-2,-8] rotate!(nets[2],i); end
breakedge!(nets[2], "K22")
for i in [-8,-9, -4,-3,-2] rotate!(nets[3],i); end
for i in [-8,-9,-10, -4] rotate!(nets[4],i); end
for i in [-5,-3, -11] rotate!(nets[5],i); end
for i in [-5,-3,-4, -7,-8,-9] rotate!(nets[6],i); end

xmax = [10,10,10,13.5,10.1,13.5]
R"pdf"(file="../../figures/snaq_networks_h2.pdf", width=9, height=4.5);
R"layout"(reshape(1:6, 2,3));
R"par"(mar=[0,0,0,0]);
for (i,n) in enumerate(nets)
    plot(n; showgamma=true, xlim=[0.4,xmax[i]])
    R"mtext"(cfinput[i], side=3, line=-1, adj=0, cex=0.8)
    isempty(details[i]) ||
    R"mtext"(details[i], side=3, line=-2.3, adj=0, cex=0.8)
end
R"dev.off"();

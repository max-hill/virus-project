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
n_9tax_97gR2L = readmultinewick(n_9tax_97gR2L_file)[3]

nets = [n_8tax_SNP, n_9tax_SNP, n_8tax_15gL2R, n_9tax_15gL2R,
    n_9tax_97gL2R, n_9tax_97gR2L]
# root with BHV5, if possible
rootonedge!(nets[1], 3)
rootonedge!(nets[2], 11)
rootonedge!(nets[3], 9)
rootatnode!(nets[4], "BHV5")
rootonedge!(nets[5], 8)
# reverse ismajor above BHV5: to show ML tree has a displayed tree
nets[5].edge[ 8].ismajor = true
nets[5].edge[11].ismajor = false
rootonedge!(nets[6], 19) # to make BHV1.1 closer to being a clade

# rotate edges to untangle them when plotting

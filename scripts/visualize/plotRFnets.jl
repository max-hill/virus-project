using PhyloNetworks, PhyloPlots
using RCall

# RF nets for different partition schemes (L2R, R2L) and RF-Net modes (regular,
# regular fast, tree-child, tree-child fast)
l2r = readMultiTopology("../../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick");
l2r_f = readMultiTopology("../../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted_f.newick");
l2r_tc = readMultiTopology("../../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted_tc.newick");
l2r_tc_f = readMultiTopology("../../results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted_tc_f.newick");

r2l = readMultiTopology("../../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick");
r2l_f = readMultiTopology("../../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_f.newick");
r2l_tc = readMultiTopology("../../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_tc.newick");
r2l_tc_f = readMultiTopology("../../results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_tc_f.newick");

###############################################################################

# l2r, 7-reticulation net, 15 genes
l2r_7ret = deepcopy(l2r[8]);
rotate!(l2r_7ret,-4); rotate!(l2r_7ret,-7); 
rotate!(l2r_7ret,-11); rotate!(l2r_7ret,-13);
rotate!(l2r_7ret,-19); rotate!(l2r_7ret,-17);
rotate!(l2r_7ret,-22); rotate!(l2r_7ret,-15);
rotate!(l2r_7ret,-5);
l2r_7ret.node[13].name = "Ti"; # abbreviate
R"pdf"(file="../../figures/l2r_7ret.pdf", width=8, height=8);
plot(l2r_7ret, majorhybridedgecolor="red", minorhybridedgecolor="blue",
     shownodelabel=true, edgewidth=2);
R"dev.off"();

# r2l, 6-reticulation net, 15 genes
r2l_6ret = deepcopy(r2l[7]);
rotate!(r2l_6ret,-3); rotate!(r2l_6ret,-7);
rotate!(r2l_6ret,-10); rotate!(r2l_6ret,-15);
rotate!(r2l_6ret,-12);
r2l_6ret.node[12].name = "Ti"; # abbreviate
R"pdf(file='../../figures/r2l_6ret.pdf',width=8,height=8)"
plot(r2l_6ret,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true,edgewidth=2)
R"dev.off()"

################################################################################

# l2r, 7-reticulation net, 97 genes
l2r_7ret_97g = readMultiTopology("../../results/bhv/rfnet/97genes_blksize1500_set1c_ogrooted.newick")[8];
l2r_7ret_97g.node[17].name = "Ti"
for i in [-3,-5,-9] rotate!(l2r_7ret_97g, i); end
# r2l, 7-reticulation net, 97 genes
r2l_7ret_97g = readMultiTopology("../../results/bhv/rfnet/97genes_blksize1500_rev_set1c_ogrooted.newick")[8];
r2l_7ret_97g.node[16].name = "Ti";
for i in [-6,-14,-16,-17,-22] rotate!(r2l_7ret_97g, i); end

# R"pdf"(file="../../figures/l2r_7ret_97g.pdf", width=8, height=8);
R"pdf"(file="../../figures/rfnet_7ret_97g.pdf", width=12, height=6);
R"layout"([1 2]);
R"par"(mar=[0,0,0,0]);
plot(l2r_7ret_97g, majorhybridedgecolor="red", minorhybridedgecolor="blue",
     shownodelabel=true, edgewidth=2, xlim=[1,21]);
# R"dev.off"();

# R"pdf"(file="../../figures/r2l_7ret_97g.pdf", width=8, height=8)
plot(r2l_7ret_97g, majorhybridedgecolor="red", minorhybridedgecolor="blue",
     shownodelabel=true, edgewidth=2, xlim=[1,23]);
R"dev.off"();

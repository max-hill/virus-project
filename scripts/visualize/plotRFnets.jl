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

# l2r, 7-reticulation net
l2r_7ret = deepcopy(l2r[8]);
rotate!(l2r_7ret,-4); rotate!(l2r_7ret,-7); 
rotate!(l2r_7ret,-11); rotate!(l2r_7ret,-13);
rotate!(l2r_7ret,-19); rotate!(l2r_7ret,-17);
rotate!(l2r_7ret,-22); rotate!(l2r_7ret,-15);
rotate!(l2r_7ret,-5);
l2r_7ret.node[13].name = "Ti"; # abbreviate
R"pdf(file='../../figures/l2r_7ret.pdf',width=8,height=8)"
plot(l2r_7ret,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true)
R"dev.off()"

# r2l, 6-reticulation net
r2l_6ret = deepcopy(r2l[7]);
rotate!(r2l_6ret,-3); rotate!(r2l_6ret,-7);
rotate!(r2l_6ret,-10); rotate!(r2l_6ret,-15);
rotate!(r2l_6ret,-12);
r2l_6ret.node[12].name = "Ti"; # abbreviate
R"pdf(file='../../figures/r2l_6ret.pdf',width=8,height=8)"
plot(r2l_6ret,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true)
R"dev.off()"
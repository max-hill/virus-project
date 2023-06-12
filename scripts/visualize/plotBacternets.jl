using DataFrames
using PhyloNetworks, PhyloPlots
using RCall

function readTopology_nwstr(nw_str)
    nw_str_1 = replace(nw_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:");
    nw_str_2 = replace(nw_str_1,r"\\\"" => s"");
    nw_str_3 = replace(nw_str_2,r"\)[\d]+:" => s"):");
    return(readTopology(nw_str_3))
end

# experimentMC3
net = readTopology_nwstr(readline("../results/bhv/bacter/experimentMC3/tree.newick",keep=false));
R"pdf(file='../figures/bacter/experimentMC3_50.pdf',width=8,height=8)"
plot(net,majorhybridedgecolor="red",minorhybridedgecolor="blue",style=:fulltree,
    edgelabel=DataFrame(edgenum=[12],edgelab=["[11364,15577]"]),
    shownodelabel=true,
    edgecex=0.7,edgelabelcolor="blue")
R"dev.off()"

# experimentMC3, rep2
net = readTopology_nwstr(readline("../results/bhv/bacter/experimentMC3_rep2/tree.newick",keep=false));
R"pdf(file='../figures/bacter/experimentMC3_rep2_50.pdf',width=8,height=8)"
plot(net,majorhybridedgecolor="red",minorhybridedgecolor="blue",style=:fulltree,
    edgelabel=DataFrame(edgenum=[20,8],edgelab=["[11363,15594]","[84167,87738]"]),
    shownodelabel=true,
    edgecex=0.7,edgelabelcolor="blue")
R"dev.off()"

# experimentMC3, rep3
net = readTopology_nwstr(readline("../results/bhv/bacter/experimentMC3_rep3/tree.newick",keep=false));
R"pdf(file='../figures/bacter/experimentMC3_rep3_50.pdf',width=8,height=8)"
plot(net,majorhybridedgecolor="red",minorhybridedgecolor="blue",
    edgelabel=DataFrame(edgenum=[3,9],edgelab=["[81192,82528]","[84165,87739]"]),
    shownodelabel=true,
    edgecex=0.7,edgelabelcolor="blue")
R"dev.off()"

# experimentMC3, rep4
net = readTopology_nwstr(readline("../results/bhv/bacter/experimentMC3_rep4/tree.newick",keep=false));
R"pdf(file='../figures/bacter/experimentMC3_rep4_50.pdf',width=8,height=8)"
plot(net,majorhybridedgecolor="red",minorhybridedgecolor="blue",
    edgelabel=DataFrame(edgenum=[3,16,22],edgelab=["[81192,82531]","[84168,87737]",
    "[11372,15447]"]),
    shownodelabel=true,
    edgecex=0.7,edgelabelcolor="blue")
R"dev.off()"

# experimentMC3, rep4a
net = readTopology_nwstr(readline("../results/bhv/bacter/experimentMC3_rep4a/tree.newick",keep=false));
R"pdf(file='../figures/bacter/experimentMC3_rep4a_50.pdf',width=8,height=8)"
plot(net,majorhybridedgecolor="red",minorhybridedgecolor="blue",
    edgelabel=DataFrame(edgenum=[17],edgelab=["[11362,15610]"]),
    shownodelabel=true,
    edgecex=0.7,edgelabelcolor="blue")
R"dev.off()"

# experimentMC3, rep5
net = readTopology_nwstr(readline("../results/bhv/bacter/experimentMC3_rep5/tree.newick",keep=false));
R"pdf(file='../figures/bacter/experimentMC3_rep5_50.pdf',width=8,height=8)"
plot(net,majorhybridedgecolor="red",minorhybridedgecolor="blue",
    edgelabel=DataFrame(edgenum=[10,21,25],edgelab=["[118003,122044]","[11383,15363]",
    "[17932,134958]"]),
    shownodelabel=true,
    edgecex=0.7,edgelabelcolor="blue")
R"dev.off()"

# experimentMC3, rep5a
net = readTopology_nwstr(readline("../results/bhv/bacter/experimentMC3_rep5a/tree.newick",keep=false));
R"pdf(file='../figures/bacter/experimentMC3_rep5a_50.pdf',width=8,height=8)"
plot(net,majorhybridedgecolor="red",minorhybridedgecolor="blue",
    edgelabel=DataFrame(edgenum=[9,16,24],edgelab=["[84173,87736]","[96873,99717]",
    "[73218,77382]"]),
    shownodelabel=true,
    edgecex=0.7,edgelabelcolor="blue")
R"dev.off()"

# experimentMC3, rep5b
net = readTopology_nwstr(readline("../results/bhv/bacter/experimentMC3_rep5b/tree.newick",keep=false));
R"pdf(file='../figures/bacter/experimentMC3_rep5b_50.pdf',width=8,height=8)"
plot(net,majorhybridedgecolor="red",minorhybridedgecolor="blue",
    edgelabel=DataFrame(edgenum=[2,15],edgelab=["[81194,82525]","[11366,15533]"]),
    shownodelabel=true,
    edgecex=0.7,edgelabelcolor="blue")
R"dev.off()"
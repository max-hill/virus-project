using PhyloNetworks, PhyloPlots
# include("RFutilities.jl")
using RCall

## Multistart (15 genetrees), set1c, L2R, 15 genes, best network
l2r_str = readline("../analysis/netrax/multistart-set1cL2R15genes/multistart-set1cL2R15genes_inferred_network.nw",keep=false);
l2r_infnet = replace(l2r_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
    readTopology; # l2r_infnet.numHybrids: 8
l2r_infnet_no3cyc = deepcopy(l2r_infnet);
shrink3cycles!(l2r_infnet_no3cyc); # l2r_infnet_no3cyc.numHybrids: 7

# findfirst("Titanium_IBR_MLV_vaccine" .== map(x -> x.name,l2r_infnet_no3cyc.node)) # 8
l2r_infnet_no3cyc.node[8].name = "Ti";

R"pdf(file='../figures/netrax/multistart-set1c-L2R-15genes-no3cyc.pdf',width=8,height=8)"
plot(l2r_infnet_no3cyc,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true,edgewidth=2)
R"dev.off()"

## Multistart (15 genetrees), set1c, R2L, 15 genes, best network
r2l_str = readline("../analysis/netrax/multistart-set1cR2L15genes/multistart-set1cR2L15genes_inferred_network.nw",keep=false);
r2l_infnet = replace(r2l_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
    readTopology; # r2l_net.numHybrids: 8
r2l_infnet_no3cyc = deepcopy(r2l_infnet);
shrink3cycles!(r2l_infnet_no3cyc); # r2l_infnet_no3cyc.numHybrids: 7

findfirst("Titanium_IBR_MLV_vaccine" .== map(x -> x.name,r2l_infnet_no3cyc.node)) # 5
r2l_infnet_no3cyc.node[5].name = "Ti";

R"pdf(file='../figures/netrax/multistart-set1c-R2L-15genes-no3cyc.pdf',width=8,height=8)"
rotate!(r2l_infnet_no3cyc,-23);
rotate!(r2l_infnet_no3cyc,18);
rotate!(r2l_infnet_no3cyc,7);
rotate!(r2l_infnet_no3cyc,23);
rotate!(r2l_infnet_no3cyc,-10);
plot(r2l_infnet_no3cyc,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true,edgewidth=2)
R"dev.off()"
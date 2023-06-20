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

###############################################################################

## Multistart (15 genetrees), set1b, L2R, 15 genes, best network
l2r_str = readline("../analysis/netrax/multistart-titanium/multistart-titanium_inferred_network.nw",keep=false);
l2r_infnet = replace(l2r_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
    (x -> replace(x,r"NODE_[\d.]+" => s"")) |> readTopology; # l2r_infnet.numHybrids: 6
l2r_infnet_no3cyc = deepcopy(l2r_infnet);
shrink3cycles!(l2r_infnet_no3cyc); # l2r_infnet_no3cyc.numHybrids: 2

findfirst("Titanium_IBR_MLV_vaccine" .== map(x -> x.name,l2r_infnet.node)) # 25
l2r_infnet.node[25].name = "Ti"

findfirst("Titanium_IBR_MLV_vaccine" .== map(x -> x.name,l2r_infnet_no3cyc.node)) # 17
l2r_infnet_no3cyc.node[17].name = "Ti";

R"pdf(file='../figures/netrax/multistart-set1b-L2R-15genes-no3cyc.pdf',width=8,height=8)"
rootonedge!(l2r_infnet_no3cyc,13);
rotate!(l2r_infnet_no3cyc,-12);
rotate!(l2r_infnet_no3cyc,-4);
rotate!(l2r_infnet_no3cyc,-8);
rotate!(l2r_infnet_no3cyc,-17);
rotate!(l2r_infnet_no3cyc,-3);
plot(l2r_infnet_no3cyc,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true,edgewidth=2)
R"dev.off()"

## Multistart (15 genetrees), set1b, L2R, 15 genes, best network, 3-cycles not removed
R"pdf(file='../figures/netrax/multistart-set1b-L2R-15genes-with3cyc.pdf',width=8,height=8)"
rootonedge!(l2r_infnet,20);
rotate!(l2r_infnet,-17);
rotate!(l2r_infnet,-18);
rotate!(l2r_infnet,-19);
rotate!(l2r_infnet,-3);
rotate!(l2r_infnet,-6);
rotate!(l2r_infnet,-8);
rotate!(l2r_infnet,-16);
plot(l2r_infnet,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true,edgewidth=2)
R"dev.off()"
using PhyloNetworks, PhyloPlots
# include("RFutilities.jl")
using RCall
using DataFrames

## Multistart (15 genetrees), set1c, L2R, 15 genes, best network
l2r_str = readline("../analysis/netrax/multistart-set1cL2R15genes/multistart-set1cL2R15genes_inferred_network.nw",keep=false);
l2r_infnet = replace(l2r_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
    readTopology; # l2r_infnet.numHybrids: 8
l2r_infnet_no3cyc = deepcopy(l2r_infnet);
shrink3cycles!(l2r_infnet_no3cyc); # l2r_infnet_no3cyc.numHybrids: 7
# findfirst("Titanium_IBR_MLV_vaccine" .== map(x -> x.name,l2r_infnet_no3cyc.node)) # 8
l2r_infnet_no3cyc.node[8].name = "Ti";

## Multistart (15 genetrees), set1c, R2L, 15 genes, best network
r2l_str = readline("../analysis/netrax/multistart-set1cR2L15genes/multistart-set1cR2L15genes_inferred_network.nw",keep=false);
r2l_infnet = replace(r2l_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
    readTopology; # r2l_net.numHybrids: 8
r2l_infnet_no3cyc = deepcopy(r2l_infnet);
shrink3cycles!(r2l_infnet_no3cyc); # r2l_infnet_no3cyc.numHybrids: 7
findfirst("Titanium_IBR_MLV_vaccine" .== map(x -> x.name,r2l_infnet_no3cyc.node)) # 5
r2l_infnet_no3cyc.node[5].name = "Ti";
for i in [-23,18,7,23,-10] rotate!(r2l_infnet_no3cyc,i); end

#= # now plot them, 2024
R"pdf"(file="../figures/netrax/multistart-set1c-L2R-15genes-no3cyc.pdf", width=8, height=8)
plot(l2r_infnet_no3cyc, majorhybridedgecolor="red", minorhybridedgecolor="blue",
     shownodelabel=true, edgewidth=2);
R"dev.off"()

R"pdf(file='../figures/netrax/multistart-set1c-R2L-15genes-no3cyc.pdf',width=8,height=8)"
plot(r2l_infnet_no3cyc,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true,edgewidth=2)
R"dev.off()"
=#

# new plots, 2025-05
for n in l2r_infnet_no3cyc.node # remove names from internal tree nodes
     (n.hybrid || n.leaf) && continue
     n.name = ""
end
for i in [-14,21,-8] rotate!(l2r_infnet_no3cyc, i); end
for n in r2l_infnet_no3cyc.node # remove names from internal tree nodes
     (n.hybrid || n.leaf) && continue
     n.name = ""
end
rootonedge!(r2l_infnet_no3cyc, 34);
shrink3cycles!(r2l_infnet_no3cyc); # 1 fewer reticulation
for i in [24,-24,23,11,18,-5] rotate!(r2l_infnet_no3cyc, i); end
#l2r_edat = DataFrame(
#  num = [e.number for e in l2r_infnet_no3cyc.edge if e.hybrid],
#  lab = [round(e.gamma, digits=2) for e in l2r_infnet_no3cyc.edge if e.hybrid])
R"pdf"(file="../figures/netrax_15genes_no3cyc.pdf", width=12, height=6);
R"layout"([1 2]);
R"par"(mar=[0,0,0,0]);
plot(l2r_infnet_no3cyc, majorhybridedgecolor="red", minorhybridedgecolor="blue",
     shownodelabel=true, edgewidth=2, showgamma=true);
plot(r2l_infnet_no3cyc, majorhybridedgecolor="red", minorhybridedgecolor="blue",
     shownodelabel=true, edgewidth=2, showgamma=true);
R"dev.off()"


###############################################################################

## Multistart (15 genetrees), set1b, L2R, 15 genes, best network
l2r_str = readline("../analysis/netrax/multistart-titanium/multistart-titanium_inferred_network.nw",keep=false);
l2r_infnet = replace(l2r_str,r"#(?<hybrid>[\d]+):" => s"#H\g<hybrid>:") |>
    (x -> replace(x,r"NODE_[\d.]+" => s"")) |> readTopology; # l2r_infnet.numHybrids: 6
findfirst("Titanium_IBR_MLV_vaccine" .== map(x -> x.name,l2r_infnet.node)) # 25
l2r_infnet.node[25].name = "Ti"

## Multistart (15 genetrees), set1b, L2R, 15 genes, best network, 3-cycles not removed
l2r_infnet_no3cyc = deepcopy(l2r_infnet);
shrink3cycles!(l2r_infnet_no3cyc); # l2r_infnet_no3cyc.numHybrids: 2
rootonedge!(l2r_infnet_no3cyc,13);
for i in [-13,-4,-8,-17,-3,-16] # [-12,-4,-8,-17,-3]
     rotate!(l2r_infnet_no3cyc,i)
end

rootonedge!(l2r_infnet,20);
for i in [-17,-18,-19,-3,-6,-8,-16]
     rotate!(l2r_infnet,i)
end
#= # now plot them
R"pdf(file='../figures/netrax/multistart-set1b-L2R-15genes-no3cyc.pdf',width=8,height=8)"
plot(l2r_infnet_no3cyc,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true,edgewidth=2)
R"dev.off()"

R"pdf(file='../figures/netrax/multistart-set1b-L2R-15genes-with3cyc.pdf',width=8,height=8)"
plot(l2r_infnet,majorhybridedgecolor="red",minorhybridedgecolor="blue",
     shownodelabel=true,edgewidth=2)
R"dev.off()"
=#

# new plots, 2025-05
l2r_infnet_no3cyc.node[findfirst(x -> x.name == "H0", l2r_infnet_no3cyc.node)].name = "H01"
R"pdf"(file="../figures/netrax_15genes_8tax_l2r.pdf", width=12, height=6);
R"layout"([1 2]);
R"par"(mar=[0,0,0,0]);
plot(l2r_infnet, majorhybridedgecolor="red", minorhybridedgecolor="blue",
     shownodelabel=true, edgewidth=2, showgamma=true);
plot(l2r_infnet_no3cyc, majorhybridedgecolor="red", minorhybridedgecolor="blue",
     shownodelabel=true, edgewidth=2, showgamma=true);
R"dev.off()"

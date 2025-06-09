using DataFrames
using RCall
using PhyloNetworks
using Revise # to use local version of PhyloPlots
using PhyloPlots # v2.0.1 + code to fix bug: `edgewidth` was not used for arrows

net0 = readnewick("(((BHV5:1.6)#H11:1.2,#H10:0):3.2,((SP1777:4,(B589:3.7,((K22:2,#H20:0):0.2)#H21:1.5):0.3):1,(((#H21:0,((#H11:0,216_II:1.6):0.4)#H20:0.2):0.6)#H10:1.2,(((Cooper:1,C46:1):1,Ti:2):1,C33:3):1):1):1);");
# arguments for plot(net0)
args = (style=:majortree, useedgelength=true, arrowlen=0.12, tipoffset=0.15,
        majorhybridedgecolor="black", edgewidth=edgewidth, xlim=[1.1,7.8])
# to add minor edges manually after plot()
minor_num = [3,14, # H10 and H11
             8,13] # H20 and H21
edgewidth = Dict(e.number => 1 for e in net0.edge)
for i in eachindex(minor_num)
    edgewidth[minor_num[i]] = 0 # to hide the minor edges when calling 'plot'
end
minor_width = [2,1.5,2,1.5]
minor_type = ["solid","longdash","solid","longdash"]
minor_color = ["deepskyblue","deepskyblue","orangered","orangered"]
# from pp[:edge_x_lo][minor_num] etc. with pp defined below
edge_x    = [4.2, 5.4, 5.0, 4.8] # :edge_x_lo and :edge_x_hi are equal
edge_y_lo = [1, 5, 4, 5]
edge_y_hi = [5, 1, 5, 4]

R"pdf"(file="figures/bootscanevents_possiblenetwork.pdf", width=5, height=4);
R"par"(mar=[0,0,0,0]);
pp = plot(net0; args...);
for i in eachindex(minor_num)
    R"arrows"(edge_x[i], edge_y_lo[i], edge_x[i], edge_y_hi[i],
              length=args[:arrowlen], angle=20, col=minor_color[i],
              lty=minor_type[i], lwd=minor_width[i])
end
R"dev.off()";


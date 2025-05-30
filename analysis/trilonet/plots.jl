#=
using Pkg # install packages, if not done ealier
Pkg.add("PhyloNetworks")
Pkg.add("PhyloPlots")
Pkg.add("Plots") # not needed here
Pkg.add("RCall") # Pkg.build("RCall")
=#
using PhyloNetworks
using PhyloPlots
# using Plots # not needed, causes to require PhyloPlots.plot() instead of plot()
using RCall
###

# Plot experiment 1
net01 = readTopology("(BHV5,((216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1)),(B589,(SP1777,K22))))root;")
R"png(filename='trilonet-experiment-01-best-network-k6.5-set1c.png', width=36, height=24, units='cm', res=300)"
PhyloPlots.plot(net01, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 2
net02 = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
R"png(filename='trilonet-experiment-02-best-network-k6.5-set1b.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net02,-3)
PhyloPlots.plot(net02, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 3
net03 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
PhyloNetworks.rotate!(net03,-3)
PhyloNetworks.rotate!(net03,-7)
PhyloNetworks.rotate!(net03,-8)
PhyloNetworks.rotate!(net03,-9)
R"png(filename='trilonet-experiment-03-best-network-k4.0-set1c.png', width=36, height=24, units='cm', res=300)"
PhyloPlots.plot(net03, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 4
net04 = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
R"png(filename='trilonet-experiment-04-best-network-k4.0-set1b.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net04,-3)
PhyloNetworks.rotate!(net04,-6)
PhyloNetworks.rotate!(net04,-7)
PhyloNetworks.rotate!(net04,-8)
PhyloPlots.plot(net04, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 5
net05 = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))root;")
R"png(filename='trilonet-experiment-05-best-network-k6.5-set1c-with-breakpoints.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net05,-3)
PhyloNetworks.rotate!(net05,-7)
PhyloPlots.plot(net05, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 6
net06 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
R"png(filename='trilonet-experiment-06-best-network-k1.0-set1c.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net06,-3)
PhyloNetworks.rotate!(net06,-7)
PhyloNetworks.rotate!(net06,-9)
PhyloPlots.plot(net06, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 7
net07 = readTopology("(BHV5,((B589,((Cooper,(Titanium,(C33,(216_II)#H1))),(SP1777,K22))),(C46,#H1)))root;")
R"png(filename='trilonet-experiment-07-best-network-k10-set1c.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net07,-3)
PhyloNetworks.rotate!(net07,-4)
PhyloNetworks.rotate!(net07,-6)
PhyloNetworks.rotate!(net07,-7)
PhyloNetworks.rotate!(net07,-8)
PhyloPlots.plot(net07, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 8
net08 = readTopology("(BHV5,((K22,(B589,(216_II,(SP1777,(C33,(Titanium,(C46,#H1))))))),(Cooper)#H1))root;")
R"png(filename='trilonet-experiment-08-best-network-k20-set1c.png', width=36, height=24, units='cm', res=300)"
for i in [-3,-4,-5,-6,-7,-8,-9,-10]
    PhyloNetworks.rotate!(net08,i)
end
PhyloPlots.plot(net08, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 9
net09 = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((Cooper,((C46,Titanium))#H2),(C33,#H2)))))root;")
R"png(filename='trilonet-experiment-09-best-network-k4.0-set1c-with-breakpoints.png', width=36, height=24, units='cm', res=300)"
# PhyloNetworks.rotate!(net09,-3)
PhyloNetworks.rotate!(net09,-7)
PhyloNetworks.rotate!(net09,-9)
PhyloNetworks.rotate!(net09,-10)
# PhyloNetworks.rotate!(net09,-12)
PhyloPlots.plot(net09, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 10
net10 = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
R"png(filename='trilonet-experiment-10-best-network-k4.0-set1b-with-breakpoints.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net10,-6)
PhyloNetworks.rotate!(net10,-9)
PhyloNetworks.rotate!(net10,-7)
PhyloPlots.plot(net10, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 11
net11 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
R"png(filename='trilonet-experiment-11-best-network-k4.0-set1c-with-major-breakpoints.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net11,-3)
PhyloNetworks.rotate!(net11,-7)
PhyloNetworks.rotate!(net11,-8)
PhyloNetworks.rotate!(net11,-9)
PhyloPlots.plot(net11, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 12
net12 = readTopology("(BHV5,((216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1)),(B589,(SP1777,K22))))root;")
R"png(filename='trilonet-experiment-12-best-network-k6.5-set1c-with-major-breakpoints.png', width=36, height=24, units='cm', res=300)"
PhyloPlots.plot(net12, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 13
net13 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
R"png(filename='trilonet-experiment-13-best-network-k2.0-set1c.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net13,-7)
PhyloNetworks.rotate!(net13,-8)
PhyloNetworks.rotate!(net13,-3)
PhyloNetworks.rotate!(net13,-10)
PhyloPlots.plot(net13, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"


# Plot experiment 14
net14 = readTopology("(BHV5,((216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1))),(B589,(SP1777,K22))))root;")
R"png(filename='trilonet-experiment-14-best-network-k8.0-set1c.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net14,-6)
PhyloNetworks.rotate!(net14,-5)
PhyloNetworks.rotate!(net14,-7)
PhyloPlots.plot(net14, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 15
net15 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
R"png(filename='trilonet-experiment-15-best-network-k0.5-set1c.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net15,-3)
PhyloNetworks.rotate!(net15,-7)
PhyloNetworks.rotate!(net15,-9)
PhyloPlots.plot(net15, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"

# Plot experiment 16
net16 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
R"png(filename='trilonet-experiment-16-best-network-k5-set1c.png', width=36, height=24, units='cm', res=300)"
PhyloNetworks.rotate!(net16,-3)
PhyloNetworks.rotate!(net16,-7)
PhyloNetworks.rotate!(net16,-8)
PhyloNetworks.rotate!(net16,-9)
PhyloPlots.plot(net16, tipcex=1.5, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true, arrowlen=0.1, style = :fulltree)
R"dev.off()"


# Plot networks for set1c with kappa =  0.5,1,2,4,5,6.5,8,10,20
# old: re-done below
R"pdf"("trilonet-effect-of-varying-kappa-on-set1c.pdf", height=10, width=10);
R"layout(matrix(1:9, nrow=3, ncol=3, byrow=TRUE))"
R"par"(mar=[1,0,1.4,.1]); 
nets =  [net15, net06, net13, net03, net16, net01, net14, net07, net08]
kappa_values = [0.5,1,2,4,5,6.5,8,10,20]
limits=[[0,9],[0,9],[0,9.8],[0,10.2],[0,10.2],[0,11],[0,10.2],[0,11.4],[0,13.4]]
for i in 1:9
    net=nets[i]
    k = kappa_values[i]
    PhyloPlots.plot(net, xlim=limits[i], tipcex=1.3, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true)
    R"mtext"("kappa = $k")
end
R"dev.off()"




# Plot networks for set1b with kappa =  0.5,1,2,4,5,6.5,8,10,20
net_k0_5 = readTopology("(BHV5,((B589,SP1777),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
net_k1 = readTopology("(BHV5,((B589,SP1777),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
net_k2 = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
net_k4 = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
net_k5 = readTopology("(BHV5,((216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1))),(SP1777,B589)))root;")
net_k6_5 = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
net_k8 = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
net_k10 = readTopology("(BHV5,(((Cooper,(Titanium,(C33,(216_II)#H1))),(B589,SP1777)),(C46,#H1)))root;")
net_k20 = readTopology("(BHV5,((SP1777,(B589,(216_II,(C33,(Titanium,(C46,#H1)))))),(Cooper)#H1))root;")
nets = [net_k0_5, net_k1, net_k2, net_k4, net_k5, net_k6_5, net_k8, net_k10, net_k20]
for j in [1,2,3,4] rotate!(nets[j],-4); end
rotate!(nets[3],-6)
rotate!(nets[4],-10)
for i in [-3,-9]  rotate!(nets[5],i); end
setgamma!(nets[6].edge[2], 0.501) # to make it major: plots 216II in better position
rootonedge!(nets[6], 17)
for i in [-3,-4,-6,-8,-9,11,-11] rotate!(nets[6],i); end # [-3,-4,-6,-8,-9]
setgamma!(nets[7].edge[2], 0.501)
rootonedge!(nets[7], 17)
for i in [-3,-4,-6,-8,-9,11,-11] rotate!(nets[7],i); end # [-3,-4,-6,-8,-9]
setgamma!(nets[9].edge[8], 0.501)
for i in [-4,-10,-9]  rotate!(nets[8],i); end # [-4,-10]
kappa_values = ["0.5","1","2","4","5","6.5","8","10","20"]
xmin = 0.7; ymin=0.9
xmax=[9,9,10.4,10.3,10.3,12,12,11.4,13.4]
ymax=[8.5,8.5,9.5,9.5,9.5,10.5,10.5,9.5,9.8]
#= # commented-out: plot combining both 8-taxon and 9-taxon networks later
R"cairo_pdf"("trilonet-effect-of-varying-kappa-on-set1b.pdf", height=5, width=5); # height=10, width=10
R"layout"(transpose(reshape(1:9,3,3)))
R"par"(mar=[0,0,0,0]); # [1,0,1.4,.1]
for (i,net) in enumerate(nets)
    k = kappa_values[i]
    plot(net, xlim=[xmin,xmax[i]], ylim=[ymin,ymax[i]]) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"()
=#

# Plot networks for set1c with kappa =  0.5,1,2,4,5,6.5,8,10,20
net_k0_5 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
net_k1 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
net_k2 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
net_k4 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
net_k5 = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
net_k6_5 = readTopology("(BHV5,((216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1)),(B589,(SP1777,K22))))root;")
net_k8 = readTopology("(BHV5,((216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1))),(B589,(SP1777,K22))))root;")
net_k10 = readTopology("(BHV5,((B589,((Cooper,(Titanium,(C33,(216_II)#H1))),(SP1777,K22))),(C46,#H1)))root;")
net_k20 = readTopology("(BHV5,((K22,(B589,(216_II,(SP1777,(C33,(Titanium,(C46,#H1))))))),(Cooper)#H1))root;")
nets_9tax = [net_k0_5, net_k1, net_k2, net_k4, net_k5, net_k6_5, net_k8, net_k10, net_k20]

for j in [1,2,3,4,5] rotate!(nets_9tax[j],-4); end
for j in [6,7] rotate!(nets_9tax[j],-10); end
setgamma!(nets_9tax[4].edge[15], 0.501)
setgamma!(nets_9tax[5].edge[15], 0.501)
setgamma!(nets_9tax[7].edge[10], 0.501)
setgamma!(nets_9tax[9].edge[9], 0.501)
rotate!(nets_9tax[3],-7)
rotate!(nets_9tax[4],-11)
rotate!(nets_9tax[5],-11)
for i in [-5,-7,-8,-3] rotate!(nets_9tax[6],i); end
for i in [-3,-9] rotate!(nets_9tax[7],i); end
for i in [-5,-11]  rotate!(nets_9tax[8],i); end # [-5,-11]
kappa_values = ["0.5","1","2","4","5","6.5","8","10","20"]
xmin = 0.7; ymin=0.9
xmax_9tax=[9,9,10.4,10.3,10.3,12,11.2,12,14.4]
ymax_9tax=[9.5,9.5,10.5,10.5,10.5,10.5,10.5,10.5,10.8]
#=
R"cairo_pdf"("trilonet-effect-of-varying-kappa-on-set1b.pdf", height=5, width=5); # height=10, width=10
R"layout"(transpose(reshape(1:9,3,3)))
R"par"(mar=[0,0,0,0]); # [1,0,1.4,.1]
for (i,net) in enumerate(nets_9tax)
    k = kappa_values[i]
    plot(net, xlim=[xmin,xmax_9tax[i]], ylim=[ymin,ymax_9tax[i]]) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off()"
=#

# combine the 2 plots above into 1 plot, with 2 x 9 panels
R"cairo_pdf"("trilonet-effect-of-varying-kappa.pdf", height=5, width=10.5);
R"layout"([1 2 3 0 10 11 12; 4 5 6 0 13 14 15; 7 8 9 0 16 17 18],
    widths=[2 2 2 1 2 2 2]);
R"par"(mar=[0,0,0,0]);
for (i,net) in enumerate(nets)
    k = kappa_values[i]
    plot(net, xlim=[xmin,xmax[i]], ylim=[ymin,ymax[i]]) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
for (i,net) in enumerate(nets_9tax)
    k = kappa_values[i]
    plot(net, xlim=[xmin,xmax_9tax[i]], ylim=[ymin,ymax_9tax[i]]) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"();

# Plot networks for experiment 19  (k=4,6.5,8, both sets 1b and 1c, with breakpoints)
R"pdf"("trilonet-with-breakpoints.pdf", height=10, width=10);
R"layout(matrix(1:6, nrow=3, ncol=2, byrow=TRUE))"
R"par"(mar=[1,0,1.4,.1]);
net_1b_k4 = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
net_1b_k65 = readTopology("(BHV5,((B589,SP1777),(216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1))))root;")
net_1b_k8 = readTopology("(BHV5,((B589,SP1777),(216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1))))root;")
net_1c_k4 = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((Cooper,((C46,Titanium))#H2),(C33,#H2)))))root;")
net_1c_k65 = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))root;")
net_1c_k8 = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))root;")
nets = [net_1b_k4, net_1c_k4, net_1b_k65, net_1c_k65, net_1b_k8, net_1c_k8]
kappa_values = [4,4,6.5,6.5,8,8]
limits=[[0,9],[0,9],[0,10],[0,10],[0,10],[0,10]]

PhyloNetworks.rotate!(nets[1],-10)
PhyloNetworks.rotate!(nets[2],-13)
PhyloNetworks.rotate!(nets[2],-7)
PhyloNetworks.rotate!(nets[4],-7)
PhyloNetworks.rotate!(nets[6],-7)

for i in 1:6
    net=nets[i]
    k = kappa_values[i]
    PhyloPlots.plot(net, xlim=limits[i], tipcex=1.3, useedgelength=false, showedgenumber=false, shownodenumber=false, showgamma=true)
    R"mtext"("kappa = $k")
end
R"dev.off()"


# Plot networks for experiment 20  (k=4,6.5,8, both sets 1b and 1c, with only major breakpoints)
R"pdf"("trilonet-with-only-major-breakpoints.pdf", height=10, width=10);
R"layout(matrix(1:6, nrow=3, ncol=2, byrow=TRUE))"
R"par"(mar=[1,0,1.4,.1]);
net_1b_k4 = readTopology("(BHV5,(((#H1,((Cooper,((C46,Titanium))#H2),(C33,#H2))),(B589,SP1777)),(216_II)#H1))root;")
net_1b_k65 = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
net_1b_k8 = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
net_1c_k4 = readTopology("(BHV5,(((#H1,((Cooper,((C46,Titanium))#H3),(C33,#H3))),((K22,(B589)#H2),(SP1777,#H2))),(216_II)#H1))root;")
net_1c_k65 = readTopology("(BHV5,(((#H1,((Cooper,(C33,(C46)#H3)),(Titanium,#H3))),((K22,(B589)#H2),(SP1777,#H2))),(216_II)#H1))root;")
net_1c_k8 = readTopology("(BHV5,(((#H1,((Cooper,(C33,(C46)#H3)),(Titanium,#H3))),((K22,(B589,#H2)),(SP1777)#H2)),(216_II)#H1))root;")
nets = [net_1b_k4, net_1c_k4, net_1b_k65, net_1c_k65, net_1b_k8, net_1c_k8]
kappa_values = [4,4,6.5,6.5,8,8]
limits=[[0,11],[0,11],[0,11],[0,11],[0,11],[0,11]]

[PhyloNetworks.rotate!(nets[1],i) for i in [-3, -6, -7]]
[PhyloNetworks.rotate!(nets[2],i) for i in [-3,-14,-6,-7]]
PhyloNetworks.rotate!(nets[3],-3)
[PhyloNetworks.rotate!(nets[4],i) for i in [-3, -7, -6, -8, -14]]
PhyloNetworks.rotate!(nets[5],-3)
PhyloNetworks.rotate!(nets[6],-3)
PhyloNetworks.rotate!(nets[6],-6)
PhyloNetworks.rotate!(nets[6],-7)
PhyloNetworks.rotate!(nets[6],-8)

for i in 1:6
    net=nets[i]
    k = kappa_values[i]
    PhyloPlots.plot(net, xlim=limits[i], tipcex=1.3, useedgelength=false, showedgenumber=false, shownodenumber=true, showgamma=true)
    R"mtext"("kappa = $k")
end
R"dev.off()"

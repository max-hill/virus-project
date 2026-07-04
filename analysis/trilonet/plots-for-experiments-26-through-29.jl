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

#_______________________________________________________________________________
#
# Trilonet experiments 26, 27, 28, and 29
#________________________________________________________________________________

# The following function extracts networks from the output files of Trilonet
function get_networks(taxonset, dir, experiment_number)
    pattern = "finalNetwork eNewick Short:"
    println("\n\n# Networks for ", taxonset, " from directory ", dir)
    results = []  # will store (kvalue, ksnippet, network)
    for file in readdir(dir)
        if occursin(taxonset, file) && endswith(file, ".txt")
            m = match(r"-k([0-9]+(?:\.[0-9]+)?)\-", file)
            m === nothing && continue
            kvalue = parse(Float64, m.captures[1])
            ksnippet = "ex" * string(experiment_number) * "_net_k" * replace(m.captures[1], "." => "_")
            path = joinpath(dir, file)
            open(path) do io
                lines = readlines(io)
                for i in 1:length(lines)-1
                    if strip(lines[i]) == pattern
                        network = strip(lines[i+1])
                        network = replace(network, "Titanium_IBR_MLV_vaccine" => "Titanium")
                        push!(results, (kvalue, ksnippet, network))
                        break
                    end
                end
            end
        end
    end
    sort!(results, by = x -> x[1])
    for (_, ksnippet, network) in results
        println(ksnippet, "_",taxonset," = readTopology(\"", network, "\")")
    end
end

# EXPERIMENT 21 NETWORKS (major breakpoint only)

# we run the following, which gives the code for the networks listed below
get_networks("set1b","experiment-26-with-major-breakpoint-only",26)

# Ex26_Networks for set1b from directory experiment-26-with-major-breakpoint-only
ex26_net_k5_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(SP1777,B589)),(216_II)#H1))root;")
ex26_net_k6_5_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(SP1777,B589)),(216_II)#H1))root;")
ex26_net_k8_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(SP1777,B589)),(216_II)#H1))root;")



# Todo: rotate and re-root the ex26_networks to make them pretty

# Save lists of the ex26_networks
ex26_nets_set1b = [ex26_net_k5_set1b, ex26_net_k6_5_set1b, ex26_net_k8_set1b]

# Make dimensions for the plotting
xmin_set1b = 0.7
ymin_set1b = 0.9
xmax_set1b = [13,13,13]
ymax_set1b = [11,11,11]
xpad = 0
kappa_values = ["5","6.5","8"]

# Untangle the networks (for visualization purposes)
for k in [-4,-5,-6,-7,-8,-9] rotate!(ex26_nets_set1b[1],k); end
for k in [-4,-5,-6,-7,-8,-9] rotate!(ex26_nets_set1b[2],k); end
for k in [-4,-5,-6,-7,-8,-9] rotate!(ex26_nets_set1b[3],k); end
for k in [-3,-4,-5] rotate!(ex26_nets_set1b[1],k); end
for k in [-3,-4,-5] rotate!(ex26_nets_set1b[2],k); end
for k in [-3,-4,-5] rotate!(ex26_nets_set1b[3],k); end

# For displaying reticulations, we may choose which edge is major and which is
# minor, since TriLoNet doesn't determine this. For reticulations involving
# 216_II, the evidence suggests that the major reticulation edge should be the
# one closest to BHV1.1, since the genome of 216_II is overwhelmingly more
# simiar to strains in that clade, so we set a value for the reticulation
# parameter value greater than .5 to reflect that fact. For other
# reticulations, we don't have clear evidence one way or the other, so the
# minor/major assigment is arbitrary in accordance with the default convention
# in PhyloPlot.
for i in [1,2,3] setgamma!(ex26_nets_set1b[i].edge[2], 0.989); end

# Plot the networks in two 9x9 panels
R"cairo_pdf"("trilonet-major-breakpoint-only-PERMUTED.pdf", height=5, width=12);
R"layout"([1 2 3],
    widths=[2 2 2]);
R"par"(mar=[0,.5,0,.5]);
for (i,net) in enumerate(ex26_nets_set1b)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1b-xpad,xmax_set1b[i]+xpad], ylim=[ymin_set1b,ymax_set1b[i]],
         showedgenumber=false, shownodenumber=false, minorlinetype="solid", arrowlen=0,
         majorhybridedgecolor="deepskyblue", minorhybridedgecolor="deepskyblue", tipcex=.8)
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"();


# EXPERIMENT 22 NETWORKS (all breakpoints)

# we run the following, which gives the code for the networks listed below
get_networks("set1b","experiment-27-with-all-conjectured-breakpoints",27)

# Networks for set1b from directory experiment-27-with-all-conjectured-breakpoints
ex27_net_k5_set1b = readTopology("(BHV5,((216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1)),(SP1777,B589)))root;")
ex27_net_k6_5_set1b = readTopology("(BHV5,((216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1)),(SP1777,B589)))root;")
ex27_net_k8_set1b = readTopology("(BHV5,((216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1)),(SP1777,B589)))root;")

# Save lists of the networks
ex27_nets_set1b = [ex27_net_k5_set1b,ex27_net_k6_5_set1b,ex27_net_k8_set1b]


# Make dimensions for the plotting
xmin_set1b = 0.7
ymin_set1b = 0.9
xmax_set1b = [13,13,13]
ymax_set1b = [11,11,11]
xpad = 0
kappa_values = ["5","6.5","8"]
                    
# Plot the networks in two 9x9 panels
R"cairo_pdf"("trilonet-with-all-breakpoints-PERMUTED.pdf", height=5, width=12);
R"layout"([1 2 3],
    widths=[2 2 2]);
R"par"(mar=[0,.5,0,.5]);
for (i,net) in enumerate(ex27_nets_set1b)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1b-xpad, xmax_set1b[i]+xpad], ylim=[ymin_set1b,ymax_set1b[i]],
         showedgenumber=false, shownodenumber=false, minorlinetype="solid", arrowlen=0,
         majorhybridedgecolor="deepskyblue", minorhybridedgecolor="deepskyblue",tipcex=.8) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"();


# EXPERIMENT 28 NETWORKS (no breakpoints)

# we run the following, which gives the code for the networks listed below
get_networks("set1b","experiment-28-with-no-breakpoints",28)
# Networks for set1b from directory experiment-28-with-no-breakpoints
ex28_net_k5_set1b = readTopology("(BHV5,(((#H1,((Cooper,(C33,(C46)#H2)),(Titanium,#H2))),(SP1777,B589)),(216_II)#H1))root;")
ex28_net_k6_5_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(SP1777,B589)),(216_II)#H1))root;")
ex28_net_k8_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(SP1777,B589)),(216_II)#H1))root;")

# Save lists of the networks
ex28_nets_set1b = [ex28_net_k5_set1b, ex28_net_k6_5_set1b, ex28_net_k8_set1b]

# Untangle the networks
rotate!(ex28_nets_set1b[1],-5)
rotate!(ex28_nets_set1b[1],-10)
rotate!(ex28_nets_set1b[1],-4)

# Make dimensions for the plotting
xmin_set1b = 0.7
ymin_set1b = 0.9
xmax_set1b = [13,13,13]
ymax_set1b = [11,11,11]
xpad = 0

                  

# Plot the networks in two 9x9 panels
kappa_values = ["5","6.5","8"]
R"cairo_pdf"("trilonet-with-no-breakpoints-PERMUTED.pdf", height=5, width=12);
R"layout"([1 2 3],
    widths=[2 2 2]);
R"par"(mar=[0,.5,0,.5]);
for (i,net) in enumerate(ex28_nets_set1b)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1b-xpad,xmax_set1b[i]+xpad],
         ylim=[ymin_set1b,ymax_set1b[i]],
         showedgenumber=true, shownodenumber=true,  minorlinetype="solid", arrowlen=0,
         majorhybridedgecolor="deepskyblue", minorhybridedgecolor="deepskyblue", tipcex=.8) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"();




# EXPERIMENT 24 NETWORKS (uniform 1500bp breakpoint mesh)

# we run the following, which gives the code for the networks listed below
get_networks("set1b","experiment-24-with-uniform-1500-bp-mesh",24)
get_networks("set1c","experiment-24-with-uniform-1500-bp-mesh",24)

# Networks for set1b from directory experiment-24-with-uniform-1500-bp-mesh
ex24_net_k0_5_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex24_net_k1_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex24_net_k2_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex24_net_k4_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
ex24_net_k5_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,(C46,(Titanium)#H1)),(C33,#H1)))))root;")
ex24_net_k6_5_set1b = readTopology("(BHV5,(((C33,(Titanium,(C46,(Cooper)#H1))),(216_II,#H1)),(B589,SP1777)))root;")
ex24_net_k8_set1b = readTopology("(BHV5,(((Cooper,(Titanium,(C33,(216_II)#H1))),(B589,SP1777)),(C46,#H1)))root;")
ex24_net_k10_set1b = readTopology("(BHV5,(((Cooper,(Titanium,(C33,(216_II)#H1))),(B589,SP1777)),(C46,#H1)))root;")
ex24_net_k20_set1b = readTopology("(BHV5,((SP1777,(B589,(216_II,(C33,(Titanium,(C46,#H1)))))),(Cooper)#H1))root;")

# Networks for set1c from directory experiment-24-with-uniform-1500-bp-mesh
ex24_net_k0_5_set1c = readTopology("(BHV5,(((B589,(K22)#H1),(SP1777,#H1)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex24_net_k1_set1c = readTopology("(BHV5,(((B589,(K22)#H1),(SP1777,#H1)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex24_net_k2_set1c = readTopology("(BHV5,(((B589,(K22)#H1),(SP1777,#H1)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex24_net_k4_set1c = readTopology("(BHV5,(((B589,(K22)#H1),(SP1777,#H1)),(216_II,((Cooper,(C33,(C46)#H2)),(Titanium,#H2)))))root;")
ex24_net_k5_set1c = readTopology("(BHV5,(((B589,(K22)#H1),(SP1777,#H1)),(216_II,((Cooper,(C46,(Titanium)#H2)),(C33,#H2)))))root;")
ex24_net_k6_5_set1c = readTopology("(BHV5,(((Cooper,(Titanium,(C33,(216_II)#H2))),(C46,#H2)),((B589,(K22)#H1),(SP1777,#H1))))root;")
ex24_net_k8_set1c = readTopology("(BHV5,(((Cooper,(Titanium,(C33,(216_II)#H2))),(C46,#H2)),((SP1777,(K22,#H1)),(B589)#H1)))root;")
ex24_net_k10_set1c = readTopology("(BHV5,(((Cooper,(Titanium,(C33,(216_II)#H1))),((SP1777,(K22,#H2)),(B589)#H2)),(C46,#H1)))root;")
ex24_net_k20_set1c = readTopology("(BHV5,((K22,(B589,(216_II,(SP1777,(C33,(Titanium,(C46,#H1))))))),(Cooper)#H1))root;")

# Save lists of the networks
ex24_nets_set1b = [ex24_net_k0_5_set1b, ex24_net_k1_set1b, ex24_net_k2_set1b, ex24_net_k4_set1b, ex24_net_k5_set1b, ex24_net_k6_5_set1b, ex24_net_k8_set1b, ex24_net_k10_set1b, ex24_net_k20_set1b]
ex24_nets_set1c = [ex24_net_k0_5_set1c, ex24_net_k1_set1c, ex24_net_k2_set1c, ex24_net_k4_set1c, ex24_net_k5_set1c, ex24_net_k6_5_set1c, ex24_net_k8_set1c, ex24_net_k10_set1c, ex24_net_k20_set1c]


# TODO: untangle networks
for i in [4,5,7,8] rotate!(ex24_nets_set1b[i],-10) end;
for i in [7,8] rotate!(ex24_nets_set1b[i],-4) end;
rotate!(ex24_nets_set1b[6],-9) end;
for i in [1,2,3,4,5] rotate!(ex24_nets_set1c[i],-7) end;
for i in [4,5,6] rotate!(ex24_nets_set1c[i],-13) end;
for i in [6,7] rotate!(ex24_nets_set1c[i],-9) end;
rotate!(ex24_nets_set1c[8],-4)
rotate!(ex24_nets_set1c[8],-13)

# Make dimensions for the plotting
xmin_set1b = 0.7
ymin_set1b = 0.9
xmax_set1b = [12,12,12,12,13,13,13,13,13] #[11,11,11.4,12.3,11.3,13,13,12.4,13.4]
ymax_set1b = [11,11,11,11,11,11,11,11,11] #[9.5,9.5,10.5,10.5,10.5,11.5,11.5,10.5,10.8]
xmin_set1c = 0.7
ymin_set1c = 0.9
xmax_set1c = [12,12,12,12,13,13,13,15.5,15.5] #[9,9,10.4,10.3,10.3,12,11.2,12,14.4]
ymax_set1c = [11,11,11,12,12,12,12,12,12] # [9.5,9.5,10.5,10.5,10.5,10.5,10.5,10.5,10.8]
xpad = 0

# Plot the networks in two 9x9 panels
kappa_values = ["0.5","1","2","4","5","6.5","8","10","20"]
R"cairo_pdf"("trilonet-with-uniform-1500-bp-mesh.pdf", height=5, width=12);
R"layout"([1 2 3 0 10 11 12;
           4 5 6 0 13 14 15;
           7 8 9 0 16 17 18],
    widths=[2 2 2 1 2 2 2]);
R"par"(mar=[0,.5,0,.5]);
for (i,net) in enumerate(ex24_nets_set1b)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1b-xpad,xmax_set1b[i]+xpad],
         ylim=[ymin_set1b,ymax_set1b[i]],
         showedgenumber=false, shownodenumber=false, minorlinetype="solid", arrowlen=0,
         majorhybridedgecolor="deepskyblue", minorhybridedgecolor="deepskyblue", tipcex=.8) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
for (i,net) in enumerate(ex24_nets_set1c)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1c-xpad,xmax_set1c[i]+xpad], ylim=[ymin_set1c,ymax_set1c[i]],
         showedgenumber=false, shownodenumber=false, minorlinetype="solid", arrowlen=0,
         majorhybridedgecolor="deepskyblue", minorhybridedgecolor="deepskyblue", tipcex=.8)# tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"();




#_______________________________________________________________________________
#
# Trilonet experiments 24
#________________________________________________________________________________

# The following function extracts networks from the output files of Trilonet
function get_networks(taxonset, dir, experiment_number)
    pattern = "finalNetwork eNewick Short:"
    println("\n\n# Networks for ", taxonset, " from directory ", dir)
    results = []  # will store (kvalue, ksnippet, network)
    for file in readdir(dir)
        if occursin(taxonset, file) && endswith(file, ".txt")
            m = match(r"-k([0-9]+(?:\.[0-9]+)?)\-", file)
            m === nothing && continue
            kvalue = parse(Float64, m.captures[1])
            ksnippet = "ex" * string(experiment_number) * "_net_k" * replace(m.captures[1], "." => "_")
            path = joinpath(dir, file)
            open(path) do io
                lines = readlines(io)
                for i in 1:length(lines)-1
                    if strip(lines[i]) == pattern
                        network = strip(lines[i+1])
                        network = replace(network, "Titanium_IBR_MLV_vaccine" => "Titanium")
                        push!(results, (kvalue, ksnippet, network))
                        break
                    end
                end
            end
        end
    end
    sort!(results, by = x -> x[1])
    for (_, ksnippet, network) in results
        println(ksnippet, "_",taxonset," = readTopology(\"", network, "\")")
    end
end
get_networks
get_networks("set1b","experiment-25-with-permuted-rows-no-breakpoints",25)

ex25_net_k5_set1b = readTopology("(BHV5,((216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1))),(SP1777,B589)))root;")
ex25_net_k6_5_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
ex25_net_k8_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")

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
# Trilonet experiments 21, 22, and 23
#________________________________________________________________________________

# The following function extracts networks from the output files of Trilonet
function get_networks(taxonset, dir)
    pattern = "finalNetwork eNewick Short:"
    println("\n\n# Networks for ", taxonset, " from directory ", dir)

    results = []  # will store (kvalue, ksnippet, network)

    for file in readdir(dir)
        if occursin(taxonset, file) && endswith(file, ".txt")

            m = match(r"-k([0-9]+(?:\.[0-9]+)?)\-", file)
            m === nothing && continue

            kvalue = parse(Float64, m.captures[1])
            ksnippet = "net_k" * replace(m.captures[1], "." => "_")

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
get_networks("set1b","experiment-21-with-major-breakpoint-only")
get_networks("set1c","experiment-21-with-major-breakpoint-only")

# Ex21_Networks for set1b from directory experiment-21-with-major-breakpoint-only
ex21_net_k0_5_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex21_net_k1_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex21_net_k2_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex21_net_k4_set1b = readTopology("(BHV5,(((#H1,((Cooper,((C46,Titanium))#H2),(C33,#H2))),(B589,SP1777)),(216_II)#H1))root;")
ex21_net_k5_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
ex21_net_k6_5_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
ex21_net_k8_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
ex21_net_k10_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
ex21_net_k20_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")


# Networks for set1c from directory experiment-21-with-major-breakpoint-only
ex21_net_k0_5_set1c = readTopology("(BHV5,(((K22,(B589)#H1),(SP1777,#H1)),(216_II,((Cooper,((C46,Titanium))#H2),(C33,#H2)))))root;")
ex21_net_k1_set1c = readTopology("(BHV5,(((K22,(B589)#H1),(SP1777,#H1)),(216_II,((Cooper,((C46,Titanium))#H2),(C33,#H2)))))root;")
ex21_net_k2_set1c = readTopology("(BHV5,(((K22,(B589)#H1),(SP1777,#H1)),(216_II,((Cooper,((C46,Titanium))#H2),(C33,#H2)))))root;")
ex21_net_k4_set1c = readTopology("(BHV5,(((#H1,((Cooper,((C46,Titanium))#H3),(C33,#H3))),((K22,(B589)#H2),(SP1777,#H2))),(216_II)#H1))root;")
ex21_net_k5_set1c = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H3))),(Cooper)#H3)),((K22,(B589)#H2),(SP1777,#H2))),(216_II)#H1))root;")
ex21_net_k6_5_set1c = readTopology("(BHV5,(((#H1,((Cooper,(C33,(C46)#H3)),(Titanium,#H3))),((K22,(B589)#H2),(SP1777,#H2))),(216_II)#H1))root;")
ex21_net_k8_set1c = readTopology("(BHV5,(((#H1,((Cooper,(C33,(C46)#H3)),(Titanium,#H3))),((K22,(B589,#H2)),(SP1777)#H2)),(216_II)#H1))root;")
ex21_net_k10_set1c = readTopology("(BHV5,((B589,(SP1777,(K22,(#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))),(216_II)#H1))root;")
ex21_net_k20_set1c = readTopology("(BHV5,((B589,(SP1777,(K22,(#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))),(216_II)#H1))root;")

# Todo: rotate and re-root the ex21_networks to make them pretty

# Save lists of the ex21_networks
ex21_nets_set1b = [ex21_net_k0_5_set1b, ex21_net_k1_set1b, ex21_net_k2_set1b, ex21_net_k4_set1b, ex21_net_k5_set1b, ex21_net_k6_5_set1b, ex21_net_k8_set1b, ex21_net_k10_set1b, ex21_net_k20_set1b]
ex21_nets_set1c = [ex21_net_k0_5_set1c, ex21_net_k1_set1c, ex21_net_k2_set1c, ex21_net_k4_set1c, ex21_net_k5_set1c, ex21_net_k6_5_set1c, ex21_net_k8_set1c, ex21_net_k10_set1c, ex21_net_k20_set1c]

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
kappa_values = ["0.5","1","2","4","5","6.5","8","10","20"]

# Rotate networks to make them pretty
for i in [1,2,3] rotate!(ex21_nets_set1b[i],-10); end
for i in [1,2,3] rotate!(ex21_nets_set1b[i],-3); end
for k in [-4,-5,-10] rotate!(ex21_nets_set1b[4],k); end
for k in [-4,-5,-6,-7,-8,-9] rotate!(ex21_nets_set1b[5],k); end
for k in [-4,-5,-6,-7,-8,-9] rotate!(ex21_nets_set1b[6],k); end
for k in [-4,-5,-6,-7,-8,-9] rotate!(ex21_nets_set1b[7],k); end
for k in [-4,-5,-6,-7,-8,-9] rotate!(ex21_nets_set1b[8],k); end
for k in [-4,-5,-6,-7,-8,-9] rotate!(ex21_nets_set1b[9],k); end
for i in [1,2,3] rotate!(ex21_nets_set1c[i],-13); end
for i in [1,2,3] rotate!(ex21_nets_set1c[i],-7); end
for k in [-4,-5,-10,-14] rotate!(ex21_nets_set1c[4],k); end
for k in [-4,-5,-14,-6,-7,-8,-9] rotate!(ex21_nets_set1c[5],k); end
for k in [-4,-5,-10,-14] rotate!(ex21_nets_set1c[6],k); end
for k in [-4,-5,-10] rotate!(ex21_nets_set1c[7],k); end
for k in [-7,-8,-9,-10,-11] rotate!(ex21_nets_set1c[8],k); end
for k in [-7,-8,-9,-10,-11] rotate!(ex21_nets_set1c[9],k); end
# good
for k in [-3,-4,-5] rotate!(ex21_nets_set1c[4],k); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1c[5],k); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1c[6],k); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1c[7],k); end
for k in [-3,-4,-5,-6,-7] rotate!(ex21_nets_set1c[8],k); end
for k in [-3,-4,-5,-6,-7] rotate!(ex21_nets_set1c[9],k); end
for i in [1,2,3] rotate!(ex21_nets_set1b[i],-3); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1b[4],k); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1b[5],k); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1b[6],k); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1b[7],k); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1b[8],k); end
for k in [-3,-4,-5] rotate!(ex21_nets_set1b[9],k); end
# better


# Plot the networks in two 9x9 panels
R"cairo_pdf"("trilonet-major-breakpoint-only.pdf", height=5, width=12);
R"layout"([1 2 3 0 10 11 12;
           4 5 6 0 13 14 15;
           7 8 9 0 16 17 18],
    widths=[2 2 2 1 2 2 2]);
R"par"(mar=[0,.5,0,.5]);
for (i,net) in enumerate(ex21_nets_set1b)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1b-xpad,xmax_set1b[i]+xpad], ylim=[ymin_set1b,ymax_set1b[i]],shownodenumber=false,tipcex=.8) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
for (i,net) in enumerate(ex21_nets_set1c)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1c-xpad,xmax_set1c[i]+xpad], ylim=[ymin_set1c,ymax_set1c[i]],shownodenumber=false,tipcex=.8) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"();


# EXPERIMENT 22 NETWORKS (all breakpoints)

# we run the following, which gives the code for the networks listed below
get_networks("set1b","experiment-22-with-all-breakpoints")
get_networks("set1c","experiment-22-with-all-breakpoints")

# Networks for set1b from directory experiment-22-with-all-breakpoints
ex22_net_k0_5_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex22_net_k1_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex22_net_k2_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex22_net_k4_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex22_net_k5_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1))))root;")
ex22_net_k6_5_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1))))root;")
ex22_net_k8_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1))))root;")
ex22_net_k10_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1))))root;")
ex22_net_k20_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")

# Networks for set1c from directory experiment-22-with-all-breakpoints
ex22_net_k0_5_set1c = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex22_net_k1_set1c = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex22_net_k2_set1c = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex22_net_k4_set1c = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((Cooper,((C46,Titanium))#H2),(C33,#H2)))))root;")
ex22_net_k5_set1c = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))root;")
ex22_net_k6_5_set1c = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))root;")
ex22_net_k8_set1c = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))root;")
ex22_net_k10_set1c = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((C33,(Titanium,(C46,#H2))),(Cooper)#H2))))root;")
ex22_net_k20_set1c = readTopology("(BHV5,(((B589,(SP1777)#H1),(K22,#H1)),(216_II,((Cooper,(C33,(C46)#H2)),(Titanium,#H2)))))root;")

# Todo: rotate and re-root the networks to make them pretty

# Save lists of the networks
ex22_nets_set1b = [ex22_net_k0_5_set1b, ex22_net_k1_set1b, ex22_net_k2_set1b, ex22_net_k4_set1b, ex22_net_k5_set1b, ex22_net_k6_5_set1b, ex22_net_k8_set1b, ex22_net_k10_set1b, ex22_net_k20_set1b]
ex22_nets_set1c = [ex22_net_k0_5_set1c, ex22_net_k1_set1c, ex22_net_k2_set1c, ex22_net_k4_set1c, ex22_net_k5_set1c, ex22_net_k6_5_set1c, ex22_net_k8_set1c, ex22_net_k10_set1c, ex22_net_k20_set1c]

# rotate the networks
rotate!(ex22_nets_set1b[2],-10)
rotate!(ex22_nets_set1b[3],-10)
rotate!(ex22_nets_set1b[4],-10)
rotate!(ex22_nets_set1b[4],-7)
rotate!(ex22_nets_set1b[9],-10)
rotate!(ex22_nets_set1c[2],-11)
rotate!(ex22_nets_set1c[3],-11)
for k in [-9,-10] rotate!(ex22_nets_set1c[4],k); end
rotate!(ex22_nets_set1c[5],-7)
rotate!(ex22_nets_set1c[6],-7)
rotate!(ex22_nets_set1c[7],-7)
rotate!(ex22_nets_set1c[8],-7)
rotate!(ex22_nets_set1c[9],-7)
rotate!(ex22_nets_set1c[9],-13)




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
kappa_values = ["0.5","1","2","4","5","6.5","8","10","20"]
                    
# Plot the networks in two 9x9 panels
R"cairo_pdf"("trilonet-with-all-breakpoints.pdf", height=5, width=12);
R"layout"([1 2 3 0 10 11 12;
           4 5 6 0 13 14 15;
           7 8 9 0 16 17 18],
    widths=[2 2 2 1 2 2 2]);
R"par"(mar=[0,.5,0,.5]);
for (i,net) in enumerate(ex22_nets_set1b)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1b-xpad,xmax_set1b[i]+xpad], ylim=[ymin_set1b,ymax_set1b[i]],shownodenumber=false,tipcex=.8) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
for (i,net) in enumerate(ex22_nets_set1c)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1c-xpad,xmax_set1c[i]+xpad], ylim=[ymin_set1c,ymax_set1c[i]],shownodenumber=false,tipcex=.8) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"();


# EXPERIMENT 23 NETWORKS (no breakpoints)

# we run the following, which gives the code for the networks listed below
get_networks("set1b","experiment-23-with-no-breakpoints")
get_networks("set1c","experiment-23-with-no-breakpoints")


# Ex23_Networks for set1b from directory experiment-23-with-no-breakpoints
ex23_net_k0_5_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex23_net_k1_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex23_net_k2_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex23_net_k4_set1b = readTopology("(BHV5,((B589,SP1777),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
ex23_net_k5_set1b = readTopology("(BHV5,((216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1))),(SP1777,B589)))root;")
ex23_net_k6_5_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
ex23_net_k8_set1b = readTopology("(BHV5,(((#H1,((C33,(Titanium,(C46,#H2))),(Cooper)#H2)),(B589,SP1777)),(216_II)#H1))root;")
ex23_net_k10_set1b = readTopology("(BHV5,(((Cooper,(Titanium,(C33,(216_II)#H1))),(B589,SP1777)),(C46,#H1)))root;")
ex23_net_k20_set1b = readTopology("(BHV5,((SP1777,(B589,(216_II,(C33,(Titanium,(C46,#H1)))))),(Cooper)#H1))root;")

# Ex23_Networks for set1c from directory experiment-23-with-no-breakpoints
ex23_net_k0_5_set1c = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex23_net_k1_set1c = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,(Cooper,(C33,(C46,Titanium))))))root;")
ex23_net_k2_set1c = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,((C46,Titanium))#H1),(C33,#H1)))))root;")
ex23_net_k4_set1c = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
ex23_net_k5_set1c = readTopology("(BHV5,((B589,(SP1777,K22)),(216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1)))))root;")
ex23_net_k6_5_set1c = readTopology("(BHV5,((216_II,((C33,(Titanium,(C46,#H1))),(Cooper)#H1)),(B589,(SP1777,K22))))root;")
ex23_net_k8_set1c = readTopology("(BHV5,((216_II,((Cooper,(C33,(C46)#H1)),(Titanium,#H1))),(B589,(SP1777,K22))))root;")
ex23_net_k10_set1c = readTopology("(BHV5,((B589,((Cooper,(Titanium,(C33,(216_II)#H1))),(SP1777,K22))),(C46,#H1)))root;")
ex23_net_k20_set1c = readTopology("(BHV5,((K22,(B589,(216_II,(SP1777,(C33,(Titanium,(C46,#H1))))))),(Cooper)#H1))root;")

# Save lists of the networks
ex23_nets_set1b = [ex23_net_k0_5_set1b, ex23_net_k1_set1b, ex23_net_k2_set1b, ex23_net_k4_set1b, ex23_net_k5_set1b, ex23_net_k6_5_set1b, ex23_net_k8_set1b, ex23_net_k10_set1b, ex23_net_k20_set1b]
ex23_nets_set1c = [ex23_net_k0_5_set1c, ex23_net_k1_set1c, ex23_net_k2_set1c, ex23_net_k4_set1c, ex23_net_k5_set1c, ex23_net_k6_5_set1c, ex23_net_k8_set1c, ex23_net_k10_set1c, ex23_net_k20_set1c]

# Rotate the networks
rotate!(ex23_nets_set1b[3],-10)
rotate!(ex23_nets_set1b[4],-10)
rotate!(ex23_nets_set1b[5],-9)
rotate!(ex23_nets_set1b[5],-3)
rotate!(ex23_nets_set1b[6],-3)
rotate!(ex23_nets_set1b[7],-3)
for k in [-3,-5,-6,-7] rotate!(ex23_nets_set1b[8],k); end
for i in [3,4,5] rotate!(ex23_nets_set1c[i],-11); end
rotate!(ex23_nets_set1c[6],-3)
rotate!(ex23_nets_set1c[7],-3)
rotate!(ex23_nets_set1c[7],-9)
rotate!(ex23_nets_set1c[8],-11)
rotate!(ex23_nets_set1c[8],-5)

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
R"cairo_pdf"("trilonet-with-no-breakpoints.pdf", height=5, width=12);
R"layout"([1 2 3 0 10 11 12;
           4 5 6 0 13 14 15;
           7 8 9 0 16 17 18],
    widths=[2 2 2 1 2 2 2]);
R"par"(mar=[0,.5,0,.5]);
for (i,net) in enumerate(ex23_nets_set1b)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1b-xpad,xmax_set1b[i]+xpad], ylim=[ymin_set1b,ymax_set1b[i]],shownodenumber=false,tipcex=.8) # tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
for (i,net) in enumerate(ex23_nets_set1c)
    k = kappa_values[i]
    plot(net, xlim=[xmin_set1c-xpad,xmax_set1c[i]+xpad], ylim=[ymin_set1c,ymax_set1c[i]],shownodenumber=false,tipcex=.8)# tipcex=1.3
    R"mtext"("κ=$k", side=3, line=-1.5, adj=0.1, cex=0.8)
end
R"dev.off"();

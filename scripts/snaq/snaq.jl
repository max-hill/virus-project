#= to run snaq, either from gene trees, or from SNPs.

assumed structure of the working directory from which to run the script:
├── gitrepo
│   ├── analysis
│   │   ├── snaq
│   │   │   ├── fromIQtrees
│   │   │   │   └── set1b
│   │   │   └── fromSNPs
│   │   └── ...
│   ├── ...
│   └── scripts
│       ├── ...
│       ├── snaq
│       └── ...
├── osf
│   ├── genetrees
│   │   ├── 15genes_blksize10000_set1b.treefile
│   │   ├── ...
│   └── msa
│       ├── ...
│       ├── set1b.fasta
│       ├── set1c.fasta
│       └── ...

how to run: copy-paste in julia session
input files: in `osf/{genetrees,msa}`
output files: in `gitrepo/analysis/snaq/`

versions used:
Julia v1.7.2
[336ed68f] CSV v0.10.4
[33ad39ac] PhyloNetworks v0.15.1
[c0d5b6db] PhyloPlots v1.0.0
[6f49c342] RCall v0.13.13

for using SNPs directly: R functions from Olave M. & Meyer A (2020):
at https://github.com/melisaolave/SNPs2CF
downloaded SNPs2CF "functions_v1.5.R" on 2022-08-15 with:
curl https://raw.githubusercontent.com/melisaolave/SNPs2CF/master/functions_v1.5.R > gitrepo/scripts/SNPs2CF.R

from sessionInfo():
R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
install.packages("foreach", repos="http://R-Forge.R-project.org") # foreach_1.5.1
install.packages("doMC", repos="http://R-Forge.R-project.org") # doMC_1.3.5
install.packages("pegas") # pegas_1.1

To get MSAs in phylip format instead of fasta format, downloaded python script:
curl https://raw.githubusercontent.com/tkchafin/scripts/master/fasta2phylip.py > gitrepo/scripts/fasta2phylip.py
Python 3.9.8
=#

using PhyloNetworks, CSV, Random, RCall, PhyloPlots
gitrepo = "/home/zergling/.emacs.d/virus-project" # change this to the absolute path of your git repo
outdir_root = "$gitrepo/analysis/snaq"
indir_root = "$gitrepo/../osf/"

#---- from 15 (unrooted) gene trees (10_000 sites/block), taxon sets 1b and 1c
# run on 2022-08-15 for set1b & set1c

# taxonset = "set1b"
taxonset = "set1b" # then re-run interactively (copy-paste) what's below.
indir = joinpath(indir_root, "genetrees")
outdir = joinpath(outdir_root, "fromIQtrees", taxonset)
isdir(outdir) || mkpath(outdir)
iqtrees = readMultiTopology("$indir/15genes_blksize10000/15genes_blksize10000_$taxonset.treefile");
iqtrees
q,t = redirect_stdio(stdout=devnull) do
    countquartetsintrees(iqtrees)
end;
df = writeTableCF(q,t) # 70 4-taxon sets with set1b, 126 4-taxon sets with set1c
cffile = joinpath(outdir, "tableCF.csv")
CSV.write(cffile, df);
iqtreeCF = readTableCF(cffile)
# starting with RF-Net tree: first topology in the file, with h=0
starttree = readTopology(joinpath(indir_root, "rfnet/networks/15genes_blksize10000",
                "15genes_blksize10000_$(taxonset)_rooted_r10.newick"))
writeTopology(starttree)
# set1b: "((((((C46,Cooper),Titanium_IBR_MLV_vaccine),C33),216_II),(B589,SP1777)),BHV5);"
# set1c: "((((((C46,Cooper),Titanium_IBR_MLV_vaccine),C33),216_II),((SP1777,B589),K22)),BHV5);"

# To plot the starting trees, run the above and then run the following commented out code:
# starttree1b = starttree
# R"pdf"("$outdir_root/snaq-starting-tree-1b.pdf", height=3.8, width=4);
# R"par"(mar=[.1,.1,.1,.1]); R"layout"([1])
# PhyloPlots.plot(starttree1b,tipcex=.7, xlim=[0.2,9])
# R"dev.off"()
# OR:
# starttree1c = starttree
# R"pdf"("$outdir_root/snaq-starting-tree-1c.pdf", height=3.8, width=4);
# R"par"(mar=[.1,.1,.1,.1]); R"layout"([1])
# PhyloPlots.plot(starttree1c,tipcex=.7, xlim=[0.2,11])
# R"dev.off"()

net = Vector{HybridNetwork}(undef,4)
if taxonset == "set1b"
    Random.seed!(765331)
else
    Random.seed!(873215)
end
seeds = round.(Int, 1_000_000 .* rand(4))
@time for i in 1:4
    hmax = i-1
    fileroot = joinpath(outdir, "net$hmax")
    isfile("$fileroot.out") && continue # to skip this hmax if completed already. but seed not con
    startT = (hmax==0 ? starttree : net[i-1])
    net[i] = snaq!(startT, iqtreeCF, hmax=hmax, filename=fileroot, runs=20, seed=seeds[i])
end

# timing, using single-processor, cecile-mac-2018 used for many other processes:
# set1b: 2.30 hours (from `grep "time" .../*.out`)
# set1c: 10.66h (from @time)
# see how the scores improve from 0 to 3 reticulations:
print([n.loglik for n in net])
# set1b: [365.1587301188808, 153.94391099485654, 77.30057399899863, 77.3005739986875]
# set1c: [513.0663771591117, 312.7311634534902, 187.5802997215012, 187.5802997214814]

# To plot these likelihoods, uncomment and run the following code:
# b=[365.1587301188808, 153.94391099485654, 77.30057399899863, 77.3005739986875]
# c=[513.0663771591117, 312.7311634534902, 187.5802997215012, 187.5802997214814]
# x=[0,1,2,3]
# import Plots
# Plots.plot(x,[b,c], label=["Set 1b" "Set 1c"], title="",linewidth=3,marker=[:hex :d])
# Plots.xlabel!("h: maximum number of reticulations")
# Plots.ylabel!("-log likelihood")
# Plots.savefig("$outdir_root/snaq-likelihoods.png")

# combine h=0:3 into 1 file:
bestnet_firstline = [readline(joinpath(outdir, "net$h.out")) for h in 0:3]
bestnet_file = joinpath(outdir,"snaq_net_0123.out")
write(bestnet_file, join(bestnet_firstline,"\n"))
# attempt to root these networks with BHV5
net = readMultiTopology(bestnet_file)
if taxonset == "set1b"
  for i in 1:2 # fails for with h=2 (i=3)
    @info "attempt to re-root at BHV5, h=$(i-1)"
    rootatnode!(net[i], "BHV5")
  end
  rootonedge!(net[3],9) # major hybrid edge above hybrid node just above BHV5
  rootonedge!(net[4],10)
  hardwiredClusterDistance(net[3],net[4],true) # 0: net3 = net4
  for i in [9,-6,-5] rotate!(net[1],i); end
  for i in [10,-7,-6,-5] rotate!(net[2],i); end
  for i in [-6,-5,11,-7] rotate!(net[3],i); end
else
  for i in 1:4
    @info "attempt to re-root at BHV5, h=$(i-1)"
    rootatnode!(net[i], "BHV5")
  end
  hardwiredClusterDistance(net[3],net[4],true) # 0: net3=net4
  for i in [10,-6,-5] rotate!(net[1],i); end
  for i in [11,-8,-7,-6] rotate!(net[2],i); end
  for i in [12,-3,-2,-12,-5,-7] rotate!(net[3],i); end
end
# plot the results
R"pdf"("$outdir_root/snaq_fromIQtree_$taxonset.pdf", height=3.8, width=12);
R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2 3])
plot(net[1], showedgelength=true, showgamma=true, xlim=[0.2, 9])
plot(net[2], showedgelength=true, showgamma=true, xlim=[0.2,11])
plot(net[3], showedgelength=true, showgamma=true, xlim=[0.2,13])
R"dev.off"()

#---- from SNPs, taxon sets 1b and 1c, using SNPs2CF
# run on on 2022-08-16 for set1b, 2022-08-17 for set1c
###
taxonset = "set1b"
#taxonset = "set1c" # then re-run interactively (copy-paste) what's below.
indir = joinpath(indir_root, "msa")
msafile_root = joinpath(indir, "$taxonset")
msafile = msafile_root * ".phylip"
if !isfile(msafile)
  run(`python3 $gitrepo/scripts/fasta2phylip.py -f $msafile_root.fasta`)
  # WARNING: had tabs separating taxon names from sequences.
  # manually changed to a single space. Tabs caused SNPs2CF to fail.
end

outdir = joinpath(outdir_root, "fromSNPs", taxonset)
isdir(outdir) || mkpath(outdir)
cffile = joinpath(outdir, "tableCF.csv")
if !isfile(cffile)
  R"source"("$gitrepo/scripts/SNPs2CF.R")
  @time R"SNPs2CF"(seqMatrix=msafile, outputName=cffile)
  # 7.8 (set1b) and 12.5 (set1c) minutes. 1 single core
end

snpCF = readTableCF(cffile)
# starting with RF-Net tree: first topology in the file, with h=0
starttree = readTopology(joinpath(indir_root, "rfnet/networks/15genes_blksize10000",
                "15genes_blksize10000_$(taxonset)_rooted_r10.newick"))
# same starting tree as when using gene trees as input
net = Vector{HybridNetwork}(undef,4)
Random.seed!( (taxonset == "set1b" ? 651001 : 985132) )
seeds = round.(Int, 1_000_000 .* rand(4))
@time for i in 1:4
    hmax = i-1
    fileroot = joinpath(outdir, "net$hmax")
    isfile("$fileroot.out") && continue # to skip this hmax if completed already. but seed not con
    startT = (hmax==0 ? starttree : net[i-1])
    net[i] = snaq!(startT, snpCF, hmax=hmax, filename=fileroot, runs=20, seed=seeds[i])
end
# timing, using single-processor, on cecile-mac-2018:
# 4.0 (set1b) and 8.7 (set1c) hours
# see how the scores improve from 0 to 3 reticulations:
print([n.loglik for n in net])
# set1b: [501.7797654686652, 295.145876646058, 233.87695612845857, 233.87695612806655]
# set1b: [660.7009851999568, 400.83826191436407, 315.80346871955686, 312.09545236068]
# combine h=0:3 into 1 file:
bestnet_firstline = [readline(joinpath(outdir, "net$h.out")) for h in 0:3]
bestnet_file = joinpath(outdir,"snaq_net_0123.out")
write(bestnet_file, join(bestnet_firstline,"\n"))
# attempt to root these networks with BHV5
net = readMultiTopology(bestnet_file)
if taxonset == "set1b"
  for i in 1:2 # fails for with h=2 (i=3)
    @info "attempt to re-root at BHV5, h=$(i-1)"
    rootatnode!(net[i], "BHV5")
  end
  rootonedge!(net[3],7) # major hybrid edge above hybrid node just above BHV5
  rootonedge!(net[4],12)
  hardwiredClusterDistance(net[3],net[4],true) # 0: net3 = net4
  for i in [9,-6,-5] rotate!(net[1],i); end
  for i in [10,-6,-5,-4] rotate!(net[2],i); end
  for i in [12,-6,-5,-4,-2] rotate!(net[3],i); end
else
  for i in 1:2 # fails for with h=2 (i=3)
    @info "attempt to re-root at BHV5, h=$(i-1)"
    rootatnode!(net[i], "BHV5")
  end
  rootonedge!(net[3],11)
  # h=3: no way to have BHV5 as outgroup with best network loglik=312.09545236068
  #      so going with second, third of 4th best, loglik in [313.08734576453514, 313.8895930203436]
  net[4] = readMultiTopology(joinpath(outdir, "net3.networks"))[2]
  rootonedge!(net[4],6) # possible to root above hybrid edge above BHV6
  hardwiredClusterDistance(net[3],net[4],true) # 13: net3 != net4
  for i in [10,-7,-5] rotate!(net[1],i); end
  for i in [11,-6,-5] rotate!(net[2],i); end
  for i in [12,-6,-5] rotate!(net[3],i); end
  for i in [13,-11,-10,-6,-5] rotate!(net[4],i); end
end
# plot the results
R"pdf"("$outdir_root/snaq_fromSNPs_$taxonset.pdf", height=3.8,
       width=(taxonset == "set1b" ? 12 : 16));
R"par"(mar=[.1,.1,.1,.1]);
R"layout"((taxonset == "set1b" ? [1 2 3] : [1 2 3 4]))
plot(net[1], showedgelength=true, showgamma=true, xlim=[0.2,11])
plot(net[2], showedgelength=true, showgamma=true, xlim=[0.2,11])
plot(net[3], showedgelength=true, showgamma=true, xlim=[0.2,12])
if taxonset == "set1c"
  plot(net[4], showedgelength=true, showgamma=true, xlim=[0.2,13])
end
R"dev.off"()

###
#= -------- bootstrap with h=1 ------------------
- not h=2 because of rooting issues, and for simplicity of interpretation
- on franklin03 to paralellize the runs, using tmux & stashticket to log out:

ssh franklin03.stat.wisc.edu
top # checking for no other users / big processes running
stashticket
tmux new-session -s bootsnaq
ssh franklin03
cd private/concordance/herpesvirus


from SNPs: done 2022-08-26/27.
failed when giving multiple processors to julia with `julia --project -p 30`
because ...wtf?... first bootstrap rep: all good.
second bootstrap rep: error saying "cannot serialize a running Task"
instead: paralellize the bootstrap reps, each with a single processor: see
`snaq_bootstrap.jl` to do a single bootstrap rep, and
`snaq_bootstrap_submit.sh` to run 8 reps serially

below: code to concatenate the 100 bootstrap reps and summarize them.
=#

using PhyloNetworks, CSV, DataFrames
using PhyloPlots, RCall
outdir_root = "$gitrepo/analysis/snaq"

function rootaboveoutgroup!(net::HybridNetwork, outgroup)
  good = true
  # direct rooting at the outgroup is possible, or
  # the direct parent of the outgroup is a hybrid
  try
    rootatnode!(net, outgroup)
  catch e1
    isa(e1, PhyloNetworks.RootMismatch) || rethrow(e1)
    directEdges!(net)
    oi = findfirst(n -> n.name == outgroup, net.node)
    pa = PhyloNetworks.getParent(net.node[oi].edge[1])
    for ntrials in eachindex(net.node)
      pa.hybrid && break
      good = false
      oi = findfirst(n -> n === pa, net.node)
      pa = PhyloNetworks.getMajorParent(net.node[oi])
    end
    pae = PhyloNetworks.getMajorParentEdge(pa)
    try
      rootonedge!(net, pae)
    catch e2
      isa(e2, PhyloNetworks.RootMismatch) || rethrow(e2)
      directEdges!(net)
      msg = "can't root above $outgroup's 1st ancestral hybrid: the major edge is still below a reticulation"
      throw(PhyloNetworks.RootMismatch(msg)) # e2.msg * msg
    end
  end
  msg = (good ? "" : "can't root above $outgroup directly: its parent isn't a hybrid")
  return(msg)
end

R"pdf"("$outdir_root/snaq_bootstrap.pdf", height=7.6, width=9);
R"par"(mar=[.1,.1,.1,.1], oma=[0,0,1,0]); R"layout"([1 2; 3 4])
for datainput in ("SNPs", "genetrees")
  for taxonset in ("set1b","set1c")
    # concatenate the 100 bootstrap networks into a single file
    outdir = joinpath(outdir_root,
                  (datainput == "genetrees" ? "fromIQtrees" : "fromSNPs"),
                  taxonset)
    infiles = joinpath.(outdir, "bootsnaq_" .* lpad.(0:99, 2, "0") .* ".out")
    bootfile = joinpath(outdir, "bootsnaq_all.out")
    if !isfile(bootfile) || filesize(bootfile)==0
      @info "write $bootfile"
      open(bootfile, "w") do bfh
        for infile in infiles
          bnet = readlines(infile, keep=true)[1] # 1 line only: 1 network
          write(bfh, bnet)
        end
      end
    end
    # read best network with h=1
    net1 = readMultiTopology(joinpath(outdir,"snaq_net_0123.out"))[2] # 2nd line: h=1
    rootatnode!(net1, "BHV5")
    rotateindices = (datainput == "SNPs" ?
                      (taxonset == "set1b" ? [-4,-6,-5] : [-7,-5]) :
                      (taxonset == "set1b" ? [-6,-5,-7] : [-6,-7,-8]))
    for i in rotateindices rotate!(net1,i); end
    # read & root bootstrap networks
    bootnet = readMultiTopology(bootfile)
    problematicreps = Int[]
    for (i,n) in enumerate(bootnet)
      try
        msg = rootaboveoutgroup!(n, "BHV5")
        msg == "" || push!(problematicreps, i)
      catch e
        isa(e, PhyloNetworks.RootMismatch) || rethrow(e)
        push!(problematicreps, i)
        @error "from $datainput, $taxonset: can't reroot bootstrap rep $i\n(e.msg)"
      end
    end
    isempty(problematicreps) ||
      @error "from $datainput, $taxonset: can't reroot bootstrap net for $(length(problematicreps)) reps: $(join(problematicreps,','))"
    # SNPs set1b: rep 19. set1c: none
    # genetrees set1b: 10 reps (6,9,16,17,23,31,36,59,74,83). set1c: 9 reps (3,6,28,31,46,82,90,97,98)
    # calculate bootstrap support for major tree in best network
    BSe_tree, tree1 = treeEdgesBootstrap(bootnet,net1);
    plot(net1, edgelabel=BSe_tree)
    R"mtext"("support for major tree edges", side=1, line=-1.5, cex=0.8)
    BSn, BSe, BSc, BSgam, BSedgenum = hybridBootstrapSupport(bootnet, net1);
    plot(net1, edgelabel=BSe[:,[:edge,:BS_hybrid_edge]],
          nodelabel=BSn[:,[:hybridnode,:BS_hybrid_samesisters]],
          nodelabelcolor="red4");
    R"mtext"("support for hybrid edges & full reticulation", side=1, line=-1.5, cex=0.8)
    R"mtext"("from $datainput, sets 1b (top) and 1c (bottom)", side=3, line=-1, outer=true)
    if datainput == "genetrees" && taxonset == "set1c"
      for i in [3,6,29,32]
        plot(bootnet[i])
        R"mtext"("bootstrap replicate $i", side=1, line=-1, cex=0.8)
      end
      R"mtext"("""example bootstrap: often (~10%) the direction of gene flow is reversed
                  (from gene trees, $taxonset)""", side=3, outer=true, line=-2)
    end
  end
end
R"dev.off"()

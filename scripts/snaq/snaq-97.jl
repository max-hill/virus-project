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
indir = joinpath(indir_root, "genetrees")



#---- from 97 (unrooted) gene trees (1500 sites/block) on taxon set set1c
# direction="L2R"

direction="L2R" # run interactively the folloing code block with direction="L2R" and direction="R2L"

### Begin code block
if direction == "L2R"
    outdir = joinpath(outdir_root, "fromIQtrees", "set1c-97genes-L2R")
    iqtrees = readMultiTopology("$indir/97genes_blksize1500/97genes_blksize1500_set1c.treefile");
elseif direction == "R2L"
    outdir = joinpath(outdir_root, "fromIQtrees", "set1c-97genes-R2L")
    iqtrees = readMultiTopology("$indir/97genes_blksize1500/97genes_blksize1500_rev_set1c.treefile");
end
isdir(outdir) || mkpath(outdir)

q,t = redirect_stdio(stdout=devnull) do
    countquartetsintrees(iqtrees)
end;
df = writeTableCF(q,t) #  126 4-taxon sets
cffile = joinpath(outdir, "tableCF.csv")
CSV.write(cffile, df);
iqtreeCF = readTableCF(cffile)
# starting with RF-Net tree: first topology in the file, with h=0
if direction == "L2R"
    starttree = readTopology(joinpath(indir_root, "rfnet/networks/97genes_blksize1500",
                                      "97genes_blksize1500_set1c_rooted_r20.newick"))
elseif direction == "R2L"
    starttree = readTopology(joinpath(indir_root, "rfnet/networks/97genes_blksize1500",
                                      "97genes_blksize1500_rev_set1c_rooted_r20.newick"))
end
writeTopology(starttree)
# L2R: ((((((Titanium_IBR_MLV_vaccine,C46),Cooper),C33),216_II),(B589,(K22,SP1777))),BHV5);
# R2L: ((((((C46,Cooper),Titanium_IBR_MLV_vaccine),C33),216_II),((K22,SP1777),B589)),BHV5);

net = Vector{HybridNetwork}(undef,4)
Random.seed!(873215)
seeds = round.(Int, 1_000_000 .* rand(4))
@time for i in 1:4
    hmax = i-1
    fileroot = joinpath(outdir, "net$hmax")
    isfile("$fileroot.out") && continue # to skip this hmax if completed already. but seed not con
    startT = (hmax==0 ? starttree : net[i-1])
    net[i] = snaq!(startT, iqtreeCF, hmax=hmax, filename=fileroot, runs=20, seed=seeds[i])
end
# timing, using zergling: (CPU(s): 12, Intel(R) Core(TM) i5-10400 CPU @ 2.90GHz)
# set1c, 97 genes:  (from `grep "time" .../*.out`)
# L2R:
# R2L: 20321 seconds ≈ 5.6 hours

# see how the scores improve from 0 to 3 reticulations:
# print([n.loglik for n in net])
# L2R: [91.3453760347144, 64.73722495478751, 57.007570441752435, 57.00757044168647]
# R2L: [120.18754736724615, 88.80499919683237, 75.308988155416, 74.21115492808421]

# combine h=0:3 into 1 file:
bestnet_firstline = [readline(joinpath(outdir, "net$h.out")) for h in 0:3]
bestnet_file = joinpath(outdir,"snaq_net_0123.out")
write(bestnet_file, join(bestnet_firstline,"\n"))
net = readMultiTopology(bestnet_file)

# Reroot as close to BHV5 as possible, untangle the networks, and plot side-by-side.
if direction == "R2L"
    rootatnode!(net[1], "BHV5")
    rootatnode!(net[2], -4)
    rootatnode!(net[3],-5)
    rootatnode!(net[4],-5)
    for i in [-3,-2,-10] rotate!(net[2],i); end
    for i in [-3,-4] rotate!(net[3],i); end
    for i in [-4] rotate!(net[4],i); end
    R"pdf"("$outdir_root/snaq_fromIQtree_set1c97genes_$(direction).pdf", height=3.8, width=16);
    R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2 3 4])
    plot(net[1], showedgelength=true, showgamma=true, xlim=[0.2, 9])
    plot(net[2], showedgelength=true, showgamma=true, xlim=[0.2,11])
    plot(net[3], showedgelength=true, showgamma=true, xlim=[0.2,13])
    plot(net[4], showedgelength=true, showgamma=true, xlim=[0.2,14])
    R"dev.off"()
    
elseif direction== "L2R"
    rootatnode!(net[1], "BHV5")
    rootatnode!(net[2], -9)
    rootatnode!(net[3], -11)
    hardwiredClusterDistance(net[3],net[4],true) # 0: net3 = net4
    for i in [-8,-7] rotate!(net[2], i); end
    for i in [-9,-10,-3,-5] rotate!(net[3], i); end
    R"pdf"("$outdir_root/snaq_fromIQtree_set1c97genes_$(direction).pdf", height=3.8, width=16);
    R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2 3])
    plot(net[1], showedgelength=true, showgamma=true, xlim=[0.2, 11])
    plot(net[2], showedgelength=true, showgamma=true, xlim=[0.2,13])
    plot(net[3], showedgelength=true, showgamma=true, xlim=[0.2,15])
    R"dev.off"()
end
### End code block










#---- from SNPs, taxon set 1c, using SNPs2CF
# run on 2023-04-16

direction = "L2R" # then re-run interactivley (copy-paste) using "R2L"

### Begin code block
indir = joinpath(indir_root, "msa")
msafile_root = joinpath(indir, "set1c")
msafile = msafile_root * ".phylip"
if !isfile(msafile)
  run(`python3 $gitrepo/scripts/fasta2phylip.py -f $msafile_root.fasta`)
  # WARNING: had tabs separating taxon names from sequences.
  # manually changed to a single space. Tabs caused SNPs2CF to fail.
end
outdir = joinpath(outdir_root, "fromSNPs", "set1c-97genes-$(direction)")
isdir(outdir) || mkpath(outdir)
cffile = joinpath(outdir, "tableCF.csv")
if !isfile(cffile)
  R"source"("$gitrepo/scripts/SNPs2CF.R")
  @time R"SNPs2CF"(seqMatrix=msafile, outputName=cffile)
end

snpCF = readTableCF(cffile)
# starting with RF-Net tree: first topology in the file, with h=0
if direction == "L2R"
    starttree = readTopology(joinpath(indir_root,
                                      "rfnet/networks/97genes_blksize1500",
                                      "97genes_blksize1500_set1c_rooted_r20.newick"))
    #((((((Titanium_IBR_MLV_vaccine,C46),Cooper),C33),216_II),(B589,(K22,SP1777))),BHV5);
elseif direction == "R2L"
    starttree = readTopology(joinpath(indir_root,
                                      "rfnet/networks/97genes_blksize1500",
                                      "97genes_blksize1500_rev_set1c_rooted_r20.newick"))
end

# same starting tree as when using gene trees as input
net = Vector{HybridNetwork}(undef,4)
Random.seed!(985132)
seeds = round.(Int, 1_000_000 .* rand(4))
@time for i in 1:4
    hmax = i-1
    fileroot = joinpath(outdir, "net$hmax")
    isfile("$fileroot.out") && continue # to skip this hmax if completed already. but seed not con
    startT = (hmax==0 ? starttree : net[i-1])
    net[i] = snaq!(startT, snpCF, hmax=hmax, filename=fileroot, runs=20, seed=seeds[i])
end
# timing, using zergling (12 procesors):
# L2R:
# R2L: 
# see how the scores improve from 0 to 3 reticulations:
print([n.loglik for n in net])


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
R"pdf"("$outdir_root/snaq_fromSNPs_set1c_97_genes_$(direction).pdf", height=3.8,
       width= 16);
R"par"(mar=[.1,.1,.1,.1]);
R"layout"( [1 2 3 4])
plot(net[1], showedgelength=true, showgamma=true, xlim=[0.2,11])
plot(net[2], showedgelength=true, showgamma=true, xlim=[0.2,11])
plot(net[3], showedgelength=true, showgamma=true, xlim=[0.2,12])
plot(net[4], showedgelength=true, showgamma=true, xlim=[0.2,13])
R"dev.off"()

### end code block

# To plot set1c in a 2x2 matrix, uncomment and run the following:
# R"pdf"("$outdir_root/snaq_fromSNPs_$(taxonset)_2by2.pdf", height=7.6,
#        width=16);
# R"par"(mar=[.1,.1,.1,.1]);
# R"layout(matrix(1:4,nrow=2,ncol=2,byrow=TRUE))"
# plot(net[1], showedgelength=true, showgamma=true, xlim=[0.2,11])
# plot(net[2], showedgelength=true, showgamma=true, xlim=[0.2,13])
# plot(net[3], showedgelength=true, showgamma=true, xlim=[0.2,12])
# plot(net[4], showedgelength=true, showgamma=true, xlim=[0.2,13])
# R"dev.off"()



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

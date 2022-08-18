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
│       ├── snaq.jl
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
outdir_root = "gitrepo/analysis/snaq"
indir_root = "osf/"

#---- from 15 (unrooted) gene trees (10_000 sites/block), taxon sets 1b and 1c
# run on 2022-08-15 for set1b & set1c

# taxonset = "set1b"
taxonset = "set1c" # then re-run interactively (copy-paste) what's below.
indir = joinpath(indir_root, "genetrees")
outdir = joinpath(outdir_root, "fromIQtrees", taxonset)
isdir(outdir) || mkpath(outdir)
iqtrees = readMultiTopology("$indir/15genes_blksize10000_$taxonset.treefile");

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
# set1b: [513.0663771591117, 312.7311634534902, 187.5802997215012, 187.5802997214814]
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
taxonset = "set1b"
taxonset = "set1c" # then re-run interactively (copy-paste) what's below.
indir = joinpath(indir_root, "msa")
msafile_root = joinpath(indir, "$taxonset")
msafile = msafile_root * ".phylip"
if !isfile(msafile)
  run(`python3 gitrepo/scripts/fasta2phylip.py -f $msafile_root.fasta`)
  # WARNING: had tabs separating taxon names from sequences.
  # manually changed to a single space. Tabs caused SNPs2CF to fail.
end

outdir = joinpath(outdir_root, "fromSNPs", taxonset)
isdir(outdir) || mkpath(outdir)
cffile = joinpath(outdir, "tableCF.csv")
if !isfile(cffile)
  R"source"("gitrepo/scripts/SNPs2CF.R")
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

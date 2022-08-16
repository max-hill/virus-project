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
install.packages("foreach", repos="http://R-Forge.R-project.org") foreach_1.5.1
install.packages("doMC", repos="http://R-Forge.R-project.org") doMC_1.3.5
install.packages("pegas", repos="http://R-Forge.R-project.org") pegas_1.1

=#

using PhyloNetworks, CSV, Random, RCall, PhyloPlots
outdir_root = "gitrepo/analysis/snaq"
indir_root = "osf/"

#---- from 15 (unrooted) gene trees (10_000 sites/block), taxon sets 1b and 1c
# run on 2022-08-15 for set 1b

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
  for i in 1:4 # fails for with h=2 (i=3)
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

#=---- from SNPs, taxon sets 1b and 1c, using SNPs2CF
source("gitrepo/scripts/SNPs2CF.R")
snp2cf <- SNPs2CF(wd=getwd(),seqMatrix, ImapName=NULL,
      rm.outgroup=FALSE, outgroupSp="outgroup",
      indels.as.fifth.state=FALSE, bootstrap=TRUE,
      boots.rep=100, outputName="SNPs2CF.csv",
      n.quartets="all", between.sp.only=FALSE,
      starting.sp.quartet=1, max.SNPs=NULL,
      max.quartets=100000, save.progress=TRUE, cores=1);
seqMatrix: example in examples/5taxa-30K_SNPs.phy
outputName: by default "SNPs2CF.csv"
later:
CF = readTableCF("SNPs2CF.csv")
=#



#= does 1 bootstrap replicate
use like this, with "SNPs" or "genetrees" for the data input and
"b" or "c" for set1b or set1c, from ~/private/concordance/herpesvirus :
last argument: go from 0 to 99 -- see associated submit file for this

julia --project --color=yes gitrepo/scripts/snaq/snaq_bootstrap.jl SNPs b 0

=#

using PhyloNetworks, CSV, DataFrames
outdir_root = "gitrepo/analysis/snaq"
indir_root = "osf/"

nboot = 1 # do 100 separately, then concatenate together: see snaq.jl for that
nruns = 30 # more than 20 before, because starting net has h=0 instead of 1
datainput = ARGS[1] # "SNPs" or "genetrees"
taxonset = "set1" * ARGS[2]
# to replace the easier option (that failed with multiple processors:)
# for datainput in ("genetrees","SNPs")
#   for taxonset in ("set1b","set1c")
irep = parse(Int, ARGS[3])

# starting network: from RF-Net, h=0
starttree = readTopology(joinpath(indir_root,
        "rfnet/networks/15genes_blksize10000",
        "15genes_blksize10000_$(taxonset)_rooted_r10.newick"))
outdir = joinpath(outdir_root,
                  (datainput == "genetrees" ? "fromIQtrees" : "fromSNPs"),
                  taxonset)
outfile = joinpath(outdir, "bootsnaq_" * lpad(irep, 2, "0"))
logfile = outfile * ".screen"
@info """bootstrap SNaQ for $taxonset, with $datainput as input.
        bootstrap replicate $irep, $nruns runs. output will go $outfile.
        check log at $logfile."""

# read in the data: bootstrap gene trees, or table of CFs with CIs
dat = nothing
if datainput == "genetrees"
    boottreefile = joinpath(indir_root, "genetrees", "15genes_blksize10000_$taxonset.ufboot")
    # don't use readBootstrapTrees because all bootstrap trees in a single file
    # WARNING: assuming the following order: all trees from 1st gene, then all trees from 2nd gene, etc.
    ngenes = 15 # ?could parse boottreefile name to get this?
    nboottrees = 1000 # IQ-Tree was run for ultra-fast bootstrap with 1000 reps
    dat = Array{Vector{HybridNetwork}}(undef, ngenes)
    open(boottreefile) do bootfh
      for igene in 1:ngenes
        dat[igene] = Vector{HybridNetwork}(undef, nboottrees)
        for j in 1:nboottrees
          trestr = readline(bootfh)
          dat[igene][j] = readTopology(trestr)
        end
      end
    end
    @info "bootstrap from $(length(dat)) loci. #trees/locus in $(extrema(length(d) for d in dat))"
else
    cffile = joinpath(outdir, "tableCF.csv")
    dat = DataFrame(CSV.File(cffile); copycols = false)
    @info "data file at $cffile"
end
# run SNaQ on bootstrap dataset
open(logfile, "w") do logio
    redirect_stdio(stdout=logio, stderr=logio)
    bootnet = bootsnaq(starttree, dat, hmax=1, nrep=nboot, runs=nruns, filename=outfile)
end

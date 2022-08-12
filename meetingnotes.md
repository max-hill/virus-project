# 2022-08-11

to do next:

- find a regular time that works for all during the Fall semester!

Write bullet-point lists for:
1. Summarize what can we conclude so far about BHV? and
   what are we still unclear about?
2. Describe the nuts and bolts of analyzing these data sets.
   What is the gold standard?
   bootscan: says that some portion of the genome from strain X came
   from strain B instead of strain A, but doesn't say when these events occured.
   Something like this was used for SARS-CoV-2.
   List the caveats of using these programs.
   What is causing the various output networks? biases? stopping criterion?

Also:
3. Shuqi: continue the `bacter` analysis
   If doing another analysis of this type: use the reassortment model,
   not the single-recombination-breakpoint model.
4. Revise the RF-Net analyses? summarize their results similarly to how
   Aaron summarized his HSV-1 results? Basically:
   * from multiple output networks (from his 3 runs), extract the major tree,
     then find the recombinations found in 2/3 of the network,
   * make a table of block x recombination, with 1 color for each separate
     recombination
   * map the blocks as consecutive linear segments of the genome, and estimate
     the density of the # of recombination.
   Possibly:
   * for each data set, run multiple independent runs, and only keep the network
     from the "best" run (best likelihood or embedding cost score)
   * create 100 bootstrap data sets (resample blocks?), re-do the same on each
   * summarize the 100 bootstrap networks like Aaron summarized his 3 runs.
     functions can help do this automatically in PhyloNetworks, to do this
     on 100 networks instead of 3: see this section to summarize a
     [bootstrap set](http://crsl4.github.io/PhyloNetworks.jl/latest/man/bootstrap/#support-for-tree-edges)

presentation & discussion: Ben & Max

blocks of 10,000 bp -> IQ-Tree -> timetree -> 15 rooted gene trees
(7 distinct unrooted topologies).
each one used as a starting tree for Net-RAX on the whole alignment,
using the default "linked" mode.
The "unlinked" option doesn't work: Max filed an issue about this.
Also, it would require forbidding triangles.
Q: is Net-RAX sensitive to the starting network?

Very different results even with similar starting trees.
8 taxa "titanium set"
over half (8) of the 15 gene trees have the same topology.
in the 15 output networks: the major tree is always one of 3 (A,B or C).
All have: (the 1.1 clade ((C46,Cooper),Ti) + C33).
Uncertainty: the last 4 taxa.
Among most output network, it's always unrooted BHV5,216II | B589 SP1717.
The (1.1,C33) node attaches to the quartet onto: 216II, or
on the internal edge (best likelihood)

- How to summarize the results?
- How many reticulations are too many?
  * use Net-RAX's BIC?
  * stopping criterion to stop the search: too lenient

70% GC-rich, silimar to HSV.

Aaron: data from Bascom-Palmer, pac-bio (Illumina has trouble with high GC areas)
mostly HSV-1, one HSV-2 outgroup.
- splitstree with gaps deleted, GTR+G+I
  32 isolates, from conjunctiva (C), cornea (K), eyelid (L) + a few more:
  ~50 strains total.
- bootscan analysis
- RF-Net analysis: 178 kb total divided into 10_000 bp blocks ->
  IQ-Tree with bootstrap for each block -> gene trees intput to RF-Net2
  3 independent runs. 4 or 5 reticulations present in all 3 runs.
  matrix block x reticulation event, color-coded
  by the way: RF-Net assumes rooted gene trees: need to root them using an
  outgroup (or some other method using branch lengths, e.g. )

Amazingly: similar recombination rate function across the genome as inferred
~8 years ago with Marc Craven, after creating recombinants in the lab, then ran
bootscan before estimating the density of recombination "hotspots".

When the HSV genome is packaged: it has many single-strand breaks => could initiate
multiple recombination events at the same time and same multi-infection.
There can also be multiple events through a single round of infection.

Coronaviruses: very different mechanism for recombination, which supposedly
leads to a single breakpoint after a recombination event.

# 2022-07-07

2 datasets: "Titanium" and "K22"
- Each has 8 taxa - BHV5, 216_II, and taxa from BHV1.1 and BHV1.2. 
- "Titanium": BHV5, 216_II, {C33, C46, Ti, Cooper}, {B589, SP1777} 
- "K22": BHV5, 216_II, {Cooper, MN3, C36}, {K22, MN2, SM023}

Partition: 14 consecutive blocks of 10,000 bp, the 15th block has 4551 bp.

Methods tried out:

1. NetRAX
- Ran NetRAX 2 times for "Titanium"
  - Both runs completed in about 12 mins and detected 7-8 reticulations. 
  - The output networks were visually very similar.
- Ran NetRAX 3 times for "K22"
  - Run-1 completed in less than a minute
  - Run-2 was terminated before completion after running for close to 8 hrs
  - Run-3 was initialized using the final state of Run-2, and completed in 10
  minutes.
  - Run-1 detected 3 reticulations, Run-2/3 detected 7 reticulations.
2. TriloNet
- TriloNet recovered (essentially) the same network as NetRAX for "Titanium".
- For "K22" the networks recovered by TriloNet and NetRAX are less easy to
compare.
3. SnappNet + PhyloNet
- SnappNet runs much slower when BHV5 (outgroup) is included, and an
"outOfMemory" error occurs after a large number of MCMC iterations.
- PhyloNet offers multiple ways to summarize networks (displayed tree,
backbone tree, maximum networks).
- For "Titanium", there are 2 max-displayed trees and 2 backbone networks with
51% support, and there is 1 max network with 25% support.
4. Bacter
- Fast to run. Its outputs have many more reticulations than the other methods.
- The authors recommend using the ACGAnnotator package in BEAST to summarize the
MCMC output. PhyloNet can't be used directly (some amount of string formatting
has to be done).
5. RF-Net + TreeTime
- TreeTime was used to root the gene trees. The method used to root the gene
trees assumes a strict molecular clock. The method failed to root some gene
trees.
- The average (across gene trees) embedding cost is much lower when 15 genes
are used than when 58 genes are used. 
- For "Titanium", if 58 genes (2,500 bp per gene) are used, then a reticulation
from BHV5 to 216_II is detected.
- For "K22", if 15 genes (10,000 bp per gene) are used, then a reticulation
from 216_II to K22 is detected.
- If we look at which genes "support" the above 2 reticulations, then site
positions match with the bootscan figure provided by Aaron. 

Next steps:
- Try NetRAX with restricted number of reticulations.
- Look more closely at the IQ-tree fitted gene trees given the partitions used.
- Add K22 to "Titanium" since we have a good idea of the topology of "Titanium".
Try adding more taxa to "Titanium" to build on the current result for that set.
- Figure out "outOfMemory" issue for SnappNet as well as the effect on runtime
of adding BHV5. For every set of taxa, consider running SnappNet with and
without BHV5 (for runtime considerations).
- Try running methods on HSV-1, HSV2 dataset with more known recombinants to
test current understanding of how these methods work?

# 2022-06-02

NetRAX
- 14 taxa, partition size 1500: didn't finish in a week.
  best inferred until then: 11 reticulations, many ancestral to BHV5
- 5 taxa, partition size 1500: best has 6 reticulations.
  best with h=1 recovers the 126_II as a hybrid (also it looks like it's BHV5)
- set1 taxa, 2500 bp: 4 taxa with Cooper, C33, C46, Titanium.
  but 3-blob: both displayed trees have same topology. explains rate variation.
- increasing to 8 taxa, 2500 bp: 30 seconds.
  again: one reticulation explains rate variation
- same 8 taxa, 500 bp: 11 minutes

SnappNet
data set 3, 7 taxa: 2 days
need tools to summarize the posterior distribution programmatically (not manually)
2 networks sampled from the posterior: don't detect 216 as a hybrid between
BHV5 and BHV1. Some from BHV1.1 (e.g. MN2) and some from BHV1.2 (e.g. Cooper, MN3)

data set on 4 taxa: C46 C33 Cooper Titanium
(C33,(Titanium,(Cooper,C46))): PP ~ 50/800
network with h=1: PP ~ 400/800 (unrooted version identical to that from RF-Net)
networks with h=2: multiple.

RF-Net
2500 bp: 58 genes with IQ-Tree < 0.5h
set3 on ~12 taxa.
reticulation with 216_II appears with later reticulations, but later than 8.
Perhaps because the block transferred from BHV5 to 216 was a small and single block.

1500 bp: 97 genes (1 of them matches the transferred segment almost perfectly)
Would need to re-root at BHV5.

set with 4 taxa: got 2 different networks from 1500 bp and 2500 bp,
but those are identical if we unroot them.
2500 bp: (C33,((Cooper,C46),Ti)) + edge from C33 to Cooper
1500 bp: (C33,((Ti,C46),Cooper)) + edge from C33 to Ti

C46Cooper | TiC33 and C33Cooper | TiC46
C46Ti | CooperC33 and C33Ti | CooperC46

gene trees:
C46Ti 16
C46Cooper 29
TiCooper 13

To get uncertainty: bootstrap the sites, or bootstrap the loci

messages for a publication:
- instability between methods
- impact of assumptions, e.g.:
  * NetRAX ignores rate variation, and can infer extra reticulations to explain
    rate variation (of the form of 3-cycles: all displayed trees have the same
    topology for differing edge lengths)
  * RF-Net ignores ILS, and may infer reticulations to explain gene tree
    discordance due to ILS
- technical difficulties, e.g.:
  * some methods infer "unrooted" networks, or rather, undirected, but does not
    say so, or does not recall this fact prominently. Multiple runs may seem
    to recover different networks (because they must be written with a root somewhere)
    but in fact equal semi-directed networks.
- model selection is difficult:
  * easiest with Bayesian methods: SnappNet
  * *what model selection method do they recommend for NetRAX? for RF-Net?*
  * slope heuristic: still just a heuristic
  * RF-Net: fast enough that we can bootstrap sites (via bootstrap trees) and/or
    loci to quantify stability of the estimate. Could bootstrapping help select
    the number of reticulations?
  * cross-validation: what criterion to use with large gene trees, to compare
    the network obtained from the training data (say 90% of the gene trees)
    with the remaining (10%) gene trees?
- lack of SNPs could affect accuracy, and choice of partition size
- computational burden and running time: for method comparisons
  *warning*: analyses may need to be re-run, so that all analyses were run on the
  same machine to make running times comparable.
- lack of identifiability: not taken into account by some methods (e.g. RF-Net and 4-cycles)
  echo the "reproducibility crisis" from Maier et al. 2022 preprint

next steps
* throw in another taxon: 5 or 6 taxa
* from BHV1.1 only, from a mix of BHV1.1 and BHV1.2
  to compare data sets with lower / higher level of variation overall,
  e.g. lower / higher number of SNPs.
* add find_graphs

# 2022-04-28

analysis of 3 data sets: 0,A,B, using SnappNet & NetRAX
* inferred 2 reticulations from 0 and A, none from B.
* SnappNet cannot handle as many taxa as NetRax.
* NetRax ran somewhat with 50 taxa, but not to completion: 50 is too many.

next steps:

- partition into loci:
  * try partition with chunks of 2500-bp?
  * or chunks of 100 SNPs? For delimiting chunks with a fixed
    number of SNPs, see the MDL part of the
    [TICR](https://github.com/nstenz/TICR#mdl) pipeline.
  * Aaron can send coordinate the genes that are known to have
    been transferred from BHV5 into 216.

- taxon sampling: sample taxa around possible reticulations:
  * around known hybrid 216 (between BHV5 and BHV1):
    to use as a positive control, if we can find that some methods
    perform well enough to detect it.
  * around MN10?
  * around the vaccine-derived strain: that was our original question.

- try RF-Net-2?
  * run either IQTree or RAxML on each gene
    (I think IQTree can do all gene alignments within a folder)
  * then use this as input to RF-Net-2

# 2022-03-31

IQ-Tree on all BHV1.1 strains + MN2
on the 7-taxon:
- full alignment: 140,000 sites. Only 231 are biallelic, on the 7-taxon subset.
- NetRax: was given blocks of 500 sites. Ran very fast.
- SnappNet: MCMC for 8 million generations, took ~8h, 1 core. Max of h=2.
- both SnappNet and NetRax returned a tree.
  SnappNet: there may be an influence from the prior on the # of reticulations.
  The set of 7 strains may have been chosen as well "separated" on the
  splitstree graph. If so, the true phylogeny for these 7 may well be a tree.

Next steps:
- use larger blocks of 5kb, to increase the average # of SNPs / block
- get gene trees from these blocks: to be used as input to NetRax, or to SNaQ
- increase number of strains with NetRax: progressively?
- try [RF-Net-2](https://github.com/flu-crew/rf-net-2)?
  [They](https://doi.org/10.1093/bioinformatics/btac075)
  analyzed 200 flu strains, which has 8 genes. It takes gene trees as input.
- Curtis & Aaron: look for data with more variability.
  ex: HCMV (Human cytomegalovirus, a beta-herpesvirus).
  Large variation, even within a host: higher substitution & recombination rates.

# 2022-01-27

- focus on BHV1.1 + MN2 as outgroup (from BHV1.2)
- run IQ-Tree. Aaron will generate the data in the desired format
- next, subset the alignment and run NetRAX then SnappNet
  * subset choice: based on disease phenotype, vaccine vs cow, and
    groups found by splitstree (phenotypes are available in box).
    From Aaron & Curtis:

        C14_CSU_034_10640
        C33
        MN5
        MN12
        MN3
        C46
        BoviShield_MLV

  * to subset the alignment: re-do it from scratch using MAFFT (Aaron)
    or use shell commands (or some external program) instead
    to keep the alignment constant
    (for later comparisons between methods, taxon sampling etc.)

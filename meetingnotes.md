# 2022-06-02

NexRAX
- 14 taxa, partition size 1500: didn't finish in a week.
  best inferred until then: 11 reticulations, many ancestral to BHV5
- 5 taxa, partition size 1500: best has 6 reticulations.
  best with h=1 recovers the 126_II as a hybrid (also it looks like it's BHV5)
- set1 taxa, 2500 bp: 4 taxa with Cooper, C33, C46, Titanium.
  but 3-blob: both displayed trees have same topology. explains rate variation.
- increasing to 8 taxa, 2500 bp: 30 seconds.
  again: one reticulation explains rate variation
- same 8 taxa, 500 bp: 11 minutes

SbappNet
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

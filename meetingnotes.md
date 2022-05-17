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

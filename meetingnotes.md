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

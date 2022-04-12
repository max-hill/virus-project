# Interpreting *snappnet* output
If we run *snappnet* on the example `JDD1.xml` file and save `stdout` to
`run_snappnet_bg.log`:
```bash
#!/bin/bash
...
#SBATCH -o run_snappnet_bg.log
...

java -jar SnappNetProjectToRun.jar JDD1.xml 2> error.txt
```
then the first few lines of `run_snappnet_bg.log` will read:
```
Seq: seq_W1943_Or3
Seq: seq_W2036_Or3
Seq: seq_B204_Jap
Seq: seq_IRIS_313-11924_Jap
Seq: seq_IRIS_313-11058_Aus
Seq: seq_IRIS_313-11737_Aus
Seq: seq_IRIS_313-11819_Ind
Seq: seq_IRIS_313-11796_Ind
Seq: seq_IRIS_313-11062_Aro
Seq: seq_IRIS_313-11825_Aro
Seq: seq_W3105_Or1
Seq: seq_W1559_Or1
6 taxa
5682 sites with weight 12000
309 patterns
```
1. The first 12 lines list the sequence ids of the 12 lineages in our dataset.
2. "6 taxa" tells us that these 12 lineages are grouped into 6 species. In fact,
the `taxonset` elements within `JDD1.xml` indicate that we have 2 lineages per
species.
```xml
<taxonset id="Or3" spec="TaxonSet">
    <taxon id="W1943_Or3" spec="Taxon"/>
    <taxon id="W2036_Or3" spec="Taxon"/>
</taxonset>
<taxonset id="Jap" spec="TaxonSet">
    <taxon id="B204_Jap" spec="Taxon"/>
    <taxon id="IRIS_313-11924_Jap" spec="Taxon"/>
</taxonset>
<taxonset id="Aus" spec="TaxonSet">
    <taxon id="IRIS_313-11058_Aus" spec="Taxon"/>
    <taxon id="IRIS_313-11737_Aus" spec="Taxon"/>
</taxonset>
<taxonset id="Ind" spec="TaxonSet">
    <taxon id="IRIS_313-11819_Ind" spec="Taxon"/>
    <taxon id="IRIS_313-11796_Ind" spec="Taxon"/>
</taxonset>
<taxonset id="Aro" spec="TaxonSet">
    <taxon id="IRIS_313-11062_Aro" spec="Taxon"/>
    <taxon id="IRIS_313-11825_Aro" spec="Taxon"/>
</taxonset>
<taxonset id="Or1" spec="TaxonSet">
    <taxon id="W3105_Or1" spec="Taxon"/>
    <taxon id="W1559_Or1" spec="Taxon"/>
</taxonset>
```
3. "5682 sites with weight 12000" tells us that:  
i. the MSA has length 12000 (i.e. each sequence is 12000 sites long)  
ii. there are only 5682 "informative" sites (i.e. non-monoallelic sites)
4. "309 patterns" tells us that there are 309-2=307 unique site patterns if we
exclude monoallelic sites (e.g. 11...1 or 00...0).

## Unique site patterns
### Single lineage per species
*Unique site pattern* has a straightforward definition when there is a single
lineage per species. Suppose we have 4 species (each with 1 lineage) and an
4-long MSA:
```
s1: 1110
s2: 1011
s3: 1001
s4: 1010
```
Site 1 is monoallelic (since all 1s). Sites 2-4 are "informative" and correspond
to 3 unique site patterns (1000, 1101, 0110). If we were to run this through
*snappnet*, we would see the output:
```
4 taxa
3 sites with weight 4
5 patterns
```
### Multiple lineages per species
When there are multiple lineages per species, *unique site pattern* is defined
so that any two site patterns that match each others on allele frequencies
within each species are considered equivalent. Suppose we extend the above
example so that species 1 now has 2 lineages and the MSA has 1 more site:
```
s1, l1: 11100
s1, l2: 11001
    s2: 10111
    s3: 10010
    s4: 10101
```
Sites 2-5 are "informative" but sites 3 and 5 are considered to have the same
pattern (11000, {10101,01101}, 00110), so the *snappnet* output would be:
```
4 taxa
4 sites with weight 5
5 patterns
```

## Default SNP-calling algorithm
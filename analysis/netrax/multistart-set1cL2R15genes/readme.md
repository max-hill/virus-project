First we make the directories

```
/analysis/netrax/multistart-set1cL2R15genes/
/analysis/netrax/multistart-set1cR2L15genes/
/analysis/netrax/multistart-set1cL2R97genes/
```

Then we make and copy the partition files to the analysis folder directories.

Then we go to the `/data/` directory and run the command

```
grep -A1 -E '>C33$|>C46$|>Titanium_IBR_MLV_vaccine$|>Cooper$|>SP1777$|>B589$|>BHV5$|>216_II$|>K22$' BHV1-plus-BHV5-outgroup-alignment.fasta | grep -v -- "^--$" > set1c.fasta
```

Then copy set1c.fasta to the appropriate experiment directories with the following, which we run from the directory `/data/`

```
cp set1c.fasta ../analysis/netrax/multistart-set1cL2R15genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cR2L15genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cL2R97genes/set1c.fasta
```

Then navigate to `/data/bhv/genetrees/ogrooted/` and run

```
cp 97genes_blksize1500_set1c.treefile ../../../../analysis/netrax/multistart-set1cL2R97genes/97genes_blksize1500_set1c.treefile

cp 15genes_blksize10000_set1c.treefile ../../../../analysis/netrax/multistart-set1cL2R15genes/15genes_blksize10000_set1c.treefile

cp 15genes_blksize10000_rev_set1c.treefile ../../../../analysis/netrax/multistart-set1cR2L15genes/15genes_blksize10000_rev_set1c.treefile

```

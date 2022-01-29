# analysis of herpesvirus DNA data

with "explicit" network methods to look for evidence of reticulation.  
data source: files from Aaron in [box](https://uwmadison.app.box.com/folder/147319895420)

first focus: BHV1.1 because
1. fairly closely related: avoid heterogeneity of GC content (base frequencies)
2. question about recombination between vaccines and natural infections.


# reproducible script
## Using SplitsTree to make implicit tree
Here we show how to use SplitsTree to create a splits tree for the bovine herpes
virus. (Source: D. H. Huson and D. Bryant, Application of Phylogenetic Networks
in Evolutionary Studies, Mol. Biol. Evol., 23(2):254-267, 2006.)

Instructions for Debian. First install SplitsTree 5. To open splitstree 5,
navigate to the directory `scripts/splitstree5/` and then run the commmand

```
./SplitsTree5
```

Then using the file menu, open the datafile `BHV1_plus BHV5
outgroup_alignment.txt` in the `data` directory.

## Using IQ-TREE to construct explicit tree
Here we use IQ-TREE to construct an explicit tree with the BHV1 data. 

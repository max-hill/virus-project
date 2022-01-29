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
Here we use IQ-TREE to construct an explicit tree with the BHV1 data. The following installation and instruction sections were copied from [this repo](https://github.com/UWMadison-computingtools-2020/fp-group-6/blob/master/readme.md)


### Installing IQ-TREE

[Download](http://www.iqtree.org) the latest version for your operating system,
move the package in a convenient place, and move the executable in your
`~/bin/` directory to easily run the program from anywhere.
Alternatively, you could make a link in your `~/bin/` folder to point
to wherever the executable is, e.g.

```shell
ln -s ~/apps/iqtree-2.1.1-MacOSX/bin/iqtree2 ~/bin/iqtree
iqtree --version
```

If your Mac refuses to run iqtree because it was download from the web, then:
choose Apple menu > System Preferences, click Security & Privacy, then click
General. You should see a note saying that iqtree was blocked. Click "allow
anyway". Then re-try running iqtree, like `iqtree --version`.

For Linux and Windows users, see [here](http://www.iqtree.org/doc/Quickstart)
for installation, using packages managers like `apt-get`, `conda`, `brew` or `pkg`.

### how to run IQ-TREE

As an example, let's say that one of your alignment block is like below,
in a file named `fake.phy` in directory `alignement/`:

    4 10
    Aa_0  ATTTGGTTAC
    Abd_0 AATTGATTCT
    Ag_0  AATTGATTAT
    Ak_1  ATTTGGTTAT

Run this command from directory `iqtree/`, say:

```shell
iqtree -s ../alignments/fake.phy -m HKY+G -T 2 -pre fake
```

IQ-TREE will create output files named `fake.*` in the current directory
because of the `-pre fake` option.
The main output file containing the tree will be in `fake.treefile`.
`fake.iqtree` would contain a quick visualization of the tree, if you are curious.
In this example the tree looks like this,
showing that Abd_0 and Ag_0 are sister:

    % cat fake.treefile
    (Aa_0:0.1072067619,(Abd_0:0.1329274929,Ag_0:0.0000010000):0.2869958964,Ak_1:0.0000010000);

In this "parenthetical description" of the tree, the numbers give branch lengths
measured in average number of DNA changes per site (i.e. per column in the alignment).


#### options to adapt

- the `-s` option is for the sequence file: adapt it to your input file name.
- adapt the `-pre` option value to change the name of output files.
- `-T 2` is to use 2 threads: adapt this to your machine.
  Do *not* use more threads than your machine has! IQ-TREE also has an option
  `-T AUTO` to automatically choose the number of threads appropriate for your machine.
  (This option was `-nt` in versions 1.x.)

Keep the option `-m HKY+G`, which describes the statistical model for sequence evolution:
the choice here is both computationally tractable and biologically flexible
(e.g. allowing for fast-evolving and slow-evolving sites).

You may also consider these other useful options:
- `--no-log` to suppress the creation of the `.log` file,  
  `--no-iqtree` to suppress the creation of the `.iqtree` file
  <!-- `--no-uniqueseq` to suppress the `.uniqueseq.phy` file --no longer created -->
- `-djc` to further suppress the output files `.bionj` and `.mldist`
- `-quiet` to suppress printing information to the screen
- `-redo` to rerun an analysis that had already been started (or finished)
- `-fast` to make the search much faster, but far less accurate --
  definitely document this option if you need to use it.

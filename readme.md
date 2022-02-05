# analysis of herpesvirus DNA data

with "explicit" network methods to look for evidence of reticulation.  
data source: files from Aaron in [box](https://uwmadison.app.box.com/folder/147319895420)

first focus: BHV1.1 because
1. fairly closely related: avoid heterogeneity of GC content (base frequencies)
2. question about recombination between vaccines and natural infections.


# Reproducible Script (Using Debian 10)
## 1. SplitsTree
Here we show how to use SplitsTree to create an implicit tree (a splits tree)
for the bovine herpes virus. Our goal in this section is to reproduce Aaron's
tree.

These instructions are for Debian 10. First install SplitsTree5, available
[here](https://software-ab.informatik.uni-tuebingen.de/download/splitstree5/welcome.html).
To open splitstree 5, which has a GUI, navigate to the directory
`scripts/splitstree5/` and then run the commmand

```
./SplitsTree5
```

Then using the file menu, open the datafile `BHV1_plus BHV5
outgroup_alignment.txt` in the `data` directory.



## 2. IQ-TREE
Here we use IQ-TREE to construct an explicit tree with the BHV1 data.

### Install IQ-TREE
Download IQ-TREE [here](http://www.iqtree.org/) and save it to the `scripts/`
directory. We used the COVID-19 release 2.1.3 (April 21, 2021) for linux, and
are using Debian 10 Buster stable. For other operating systems, see
[here](http://www.iqtree.org/doc/Quickstart). Step-by-step instruction for Mac
can also be found
[here](https://github.com/UWMadison-computingtools-2020/fp-group-6/blob/master/stepsinstructions.md#install-iq-tree).

Navigate to the scripts directory. To extract IQ-TREE from the archive, run the
command 

`tar -xf iqtree-2.1.3-Linux.tar.gz`

Run the following command to give the iqtree2 file executable permission:

`sudo chmod a+rx iqtree-2.1.3-Linux/bin/iqtree2`

Then copy the file iqtree2 to `/usr/local/bin` with the command:

`sudo cp iqtree-2.1.3-Linux/bin/iqtree2 /usr/local/bin`

We are now able to run IQ-TREE by running the command `iqtree2` in any directory.

### Download Data
We used the datafile [BHV1_plus BHV5
outgroup_alignment.txt](https://uwmadison.box.com/s/phga81fq7jacu80cc42m0e4866wfnzbg)
(7mb fasta file). Download this file to the `data/` directory. Rename the file to 
`BHV1-plus-BHV5-outgroup-alignment.fasta`.

### Running IQ-TREE

From the `data/` directory, run the command

`iqtree2 -nt AUTO -s BHV1-plus-BHV5-outgroup-alignment.fasta`

This should take a few minutes to run; when it is done, the output files will be indicated:

	```
	Analysis results written to: 
	  IQ-TREE report:                BHV1-plus-BHV5-outgroup-alignment.fasta.iqtree
	  Maximum-likelihood tree:       BHV1-plus-BHV5-outgroup-alignment.fasta.treefile
	  Likelihood distances:          BHV1-plus-BHV5-outgroup-alignment.fasta.mldist
	  Screen log file:               BHV1-plus-BHV5-outgroup-alignment.fasta.log

	```

Then move all output files to the `data` directory.

### Visualizing the Tree
Here we attempt to visualize the tree made by IQ-TREE, which is located in the file `BHV1-plus-BHV5-outgroup-alignment.fasta.treefile`

#### Attempt 1: FigTree
For our first attempt at data visualization, we will use a program called
[FigTree](https://github.com/rambaut/figtree/).

##### Installing FigTree
Download
[FigTree_v1.4.4.tgz](https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree_v1.4.4.tgz)
to the `scripts/` directory. Next, navigate to the `scripts/` directory and
extract the archive by running the command

`tar -xf FigTree_v1.4.4.tgz`

There is an error with this version of FigTree, so we will follow the steps
found
[here](http://labsergen.langebio.cinvestav.mx/bioinformatics/jacob/?p=1200) to
patch the figtree shell script in order to get it to run. The script we need to
edit is a shell script called `figtree` located in `scripts/FigTree_v1.4.4/bin`.
Navigate to the directory `scripts/FigTree_v1.4.4/bin` and run the following to
edit the file automatically:

    ```
	echo '#!/bin/sh' > FigTree_v1.4.4/bin/figtree
    echo 'FIGTREE_HOME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."' >> FigTree_v1.4.4/bin/figtree
    echo 'java -Xms64m -Xmx512m -jar $FIGTREE_HOME/lib/figtree.jar $*' >> FigTree_v1.4.4/bin/figtree
    ```

##### Running FigTree
We can now run Figtree by running the following command from the `scripts/` directory:

`bash FigTree_v1.4.4/bin/figtree`

FigTree has a GUI. Using the File menu, open the treefile 
`BHV1-plus-BHV5-outgroup-alignment.fasta.treefile`
from the `data/` directory.

And there you have a tree!

One problem with the tree here is that the evolutionary distance of our outgroup
BHV5 is so great that the other parts of the tree are too small. Indeed, going
in to the file `BHV1-plus-BHV5-outgroup-alignment.fasta.treefile` we see the
edge length for BHV5 is 0.3570603406 whereas we have many other edges on the
order of 0.00001. In order to clearly show these we will edite the treefile to
change the edge length of BHV5 from 0.3570603406 to .01. To do this run the
following command from the `data/` directory:

`sed --regexp-extended 's/BHV5:0.[0-9]+/BHV5:0.01/' BHV1-plus-BHV5-outgroup-alignment.fasta.treefile > BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile`

Alternatively, this edit can be done manually, but make sure to save the edited treefile as 

`BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile` 


Next, we open this edited file with FigTree. We export the file in svg and pdf
format to the `analysis/figtree/` directory with the following names:

	```
	BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile.pdf
	BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile.svg
	BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile-unrooted.pdf
	BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile-unrooted.svg
	```


#### Attempt 2: IcyTree
Our second attempt at data visualization used the online tool
[IcyTree](https://icytree.org/). This tool has a GUI, so we just load the file
`BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile` and then can play
around with the options. Trees produced using this tool are found in the
directory `analysis/icytree/`.

# SnappNet and NetRax
Next we will try to obtain a tree for a few select virus strains chosen from the
BHV1 dataset. The following are the strains chosen for the initial SnappNet and
NetRax analysis:

    ```
	C14_CSU_034_10640
	C33
	MN5
	MN12
	MN3
	C46
	BoviShield_MLV
	```
## Data Processing
Here we describe how to restrict the large dataset `BHV1-plus-BHV5-outgroup-alignment.fasta` to just the strains that we are intersted in. 


# Appendix: Help with IQ-TREE
This material is copied from  [here](https://github.com/UWMadison-computingtools-2020/fp-group-6/blob/master/stepsinstructions.md#install-iq-tree).

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

IQ-TREE will create output files named `fake.*` in the current directory because
of the `-pre fake` option. The main output file containing the tree will be in
`fake.treefile`. `fake.iqtree` would contain a quick visualization of the tree,
if you are curious. In this example the tree looks like this, showing that Abd_0
and Ag_0 are sister:

    % cat fake.treefile
    (Aa_0:0.1072067619,(Abd_0:0.1329274929,Ag_0:0.0000010000):0.2869958964,Ak_1:0.0000010000);

In this "parenthetical description" of the tree, the numbers give branch lengths
measured in average number of DNA changes per site (i.e. per column in the
alignment).


#### options to adapt

- the `-s` option is for the sequence file: adapt it to your input file name.
- adapt the `-pre` option value to change the name of output files.
- `-T 2` is to use 2 threads: adapt this to your machine. Do *not* use more
  threads than your machine has! IQ-TREE also has an option `-T AUTO` to
  automatically choose the number of threads appropriate for your machine. (This
  option was `-nt` in versions 1.x.)

Keep the option `-m HKY+G`, which describes the statistical model for sequence
evolution: the choice here is both computationally tractable and biologically
flexible (e.g. allowing for fast-evolving and slow-evolving sites).

You may also consider these other useful options:
- `--no-log` to suppress the creation of the `.log` file,  
  `--no-iqtree` to suppress the creation of the `.iqtree` file
  <!-- `--no-uniqueseq` to suppress the `.uniqueseq.phy` file --no longer created -->
- `-djc` to further suppress the output files `.bionj` and `.mldist`
- `-quiet` to suppress printing information to the screen
- `-redo` to rerun an analysis that had already been started (or finished)
- `-fast` to make the search much faster, but far less accurate -- definitely
  document this option if you need to use it.

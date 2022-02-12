# Analysis of Herpesvirus DNA Data

with "explicit" network methods to look for evidence of reticulation.  
data source: files from Aaron in [box](https://uwmadison.app.box.com/folder/147319895420)

first focus: BHV1.1 because
1. fairly closely related: avoid heterogeneity of GC content (base frequencies)
2. question about recombination between vaccines and natural infections.






# Reproducible Script
Unless otherwise noted, the following instructions are for Debian 10 Buster
(Stable). The project (so far) is divided into four parts:

- In Part 0, we describe how to obtain the main dataset and how to construct
  subset datasets from it.
- In Part 1, we use SplitsTree to obtain an implicit tree for the BHV1 dataset. 
- In Part 2, we use IQ-Tree in combination with data visualization software to
  obtain an explicit tree for the BHV1 dataset.
- In Part 3, we will use SnappNet or NetRax to obtain a phylogeny for a subset
  of six taxa taken from the BHV1 dataset.

## Part 0. Obtaining the Data

### Downloading the Full BHV1 Dataset
We used the datafile [BHV1\_plus\ BHV5\
outgroup\_alignment.txt](https://uwmadison.box.com/s/phga81fq7jacu80cc42m0e4866wfnzbg)
(7mb fasta file). Download this file to the `data/` directory. Rename the file to 

`BHV1-plus-BHV5-outgroup-alignment.fasta`

Throughout this document, this file will be referred to as the *BHV1 Dataset*.
This is a fasta file without any linebreaks within the sequences. We will work
with this full BHV1 dataset directly as well as with subsets of it.

### Subsetting the BHV1 Dataset
In general, to construct a dataset using only a subset of the taxa from the BHV1
dataset, there are two options:

First, use the R code in [this file](scripts/Rutilies.Rmd).

Alternatively, one can use the following shell command. To build a dataset with
three specified sequences with names "SEQ1" "SEQ2" and "SEQ3", run the following
shell command from the `data/` directory:

`grep -A1 -E '>SEQ1|>SEQ2|>SEQ3' BHV1-plus-BHV5-outgroup-alignment.fasta | grep -v -- "^--$" > OUTPUT.fasta`

This will output the dataset as `OUTPUT.fasta` to the `data/` directory.

## Part 1. SplitsTree
Here we show how to use SplitsTree to create an implicit tree (a splits tree)
for the bovine herpes virus. Our goal in this section is to reproduce Aaron's
tree.

First install SplitsTree5, available
[here](https://software-ab.informatik.uni-tuebingen.de/download/splitstree5/welcome.html).
To open splitstree 5, which has a GUI, navigate to the directory
`scripts/splitstree5/` and then run the commmand

```
./SplitsTree5
```

Then using the file menu, open the datafile `BHV1_plus BHV5
outgroup_alignment.txt` in the `data` directory. Then fiddle around with the GUI
until you are satisfied.


## Part 2. IQ-TREE
Here we use IQ-TREE to construct an explicit tree with the BHV1 data. In
addition to the detailed instructions below, additional material for help and options
when using IQ-TREE can be found
[here](https://github.com/UWMadison-computingtools-2020/fp-group-6/blob/master/stepsinstructions.md#install-iq-tree).

### Install IQ-TREE
Download IQ-TREE [here](http://www.iqtree.org/) and save it to the `scripts/`
directory. We used the COVID-19 release 2.1.3 (April 21, 2021) for Linux.
Instructions for downloading and using IQ-TREE on other operatins systems can be
found [here](http://www.iqtree.org/doc/Quickstart) (general instructions) and
[here](https://github.com/UWMadison-computingtools-2020/fp-group-6/blob/master/stepsinstructions.md#install-iq-tree)
(for Mac OS).

Once you have downloaded IQ-TREE, it needs to be extract. To do this, run the
following command from the `scripts/` directory:

`tar -xf iqtree-2.1.3-Linux.tar.gz`

Run the following command to give the iqtree2 file executable permission:

`sudo chmod a+rx iqtree-2.1.3-Linux/bin/iqtree2`

Then copy the file iqtree2 to `/usr/local/bin` with the command:

`sudo cp iqtree-2.1.3-Linux/bin/iqtree2 /usr/local/bin`

We are now able to run IQ-TREE by running the command `iqtree2` in any directory.



### Running IQ-TREE
To run IQ-TREE with the full BHV1 dataset, run the following command from the `data/`
directory:

`iqtree2 -nt AUTO -s BHV1-plus-BHV5-outgroup-alignment.fasta`

This should take a few minutes to run; when it is done, the output files will be indicated:

	```
	Analysis results written to: 
	  IQ-TREE report:                BHV1-plus-BHV5-outgroup-alignment.fasta.iqtree
	  Maximum-likelihood tree:       BHV1-plus-BHV5-outgroup-alignment.fasta.treefile
	  Likelihood distances:          BHV1-plus-BHV5-outgroup-alignment.fasta.mldist
	  Screen log file:               BHV1-plus-BHV5-outgroup-alignment.fasta.log

	```

Then move all of the output files to the `data/` directory.

### Visualizing the IQ-Tree Output
Here we show how to visualize the tree made by IQ-TREE, which is located in the
file `BHV1-plus-BHV5-outgroup-alignment.fasta.treefile`. So far we have
considered two methods: FigTree and IcyTree ("I see tree", get it? Heh..).

#### Method 1: FigTree
For our first attempt at data visualization, we will use a program called
[FigTree](https://github.com/rambaut/figtree/).

##### Installing and Patching FigTree
Download
[FigTree_v1.4.4.tgz](https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree_v1.4.4.tgz)
to the `scripts/` directory. Next, navigate to the `scripts/` directory and
extract the archive by running the command

`tar -xf FigTree_v1.4.4.tgz`

Unfortunately there is an error with this version of FigTree (possibly only on
Linux?), so we will follow the steps found
[here](http://labsergen.langebio.cinvestav.mx/bioinformatics/jacob/?p=1200) to
get it to work. We will need to edit the shell script
`scripts/FigTree_v1.4.4/bin/figtree`. To do this, navigate to the directory
`scripts/FigTree_v1.4.4/bin` and run the following to edit the file
automatically:

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
order of 0.00001. In order to clearly show these we will edit the treefile to
change the edge length of BHV5 from 0.3570603406 to .01. To do this run the
following command from the `data/` directory:

`sed --regexp-extended 's/BHV5:0.[0-9]+/BHV5:0.01/' BHV1-plus-BHV5-outgroup-alignment.fasta.treefile > BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile`

This will produce a new file, called

`BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile`

that will be sufficient for now. Alternatively, this edit can be done manually,
but make sure to save the edited treefile as

`BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile` 

so you don't overwrite the original treefile.

Next, we open this edited file with FigTree. We export the file in svg and pdf
format to the `analysis/figtree/` directory with the following names:

	```
	BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile.pdf
	BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile.svg
	BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile-unrooted.pdf
	BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile-unrooted.svg
	```

#### Method 2: IcyTree
Our second attempt at data visualization used the online tool
[IcyTree](https://icytree.org/). This browser-based tool has a GUI, so we just
load the file `BHV1-plus-BHV5-outgroup-alignment-EDITED.fasta.treefile` and then
can play around with the options. Trees produced using this tool are found in
the directory `analysis/icytree/`.

# Part 3. SnappNet and NetRax

## Obtain the Dataset
The following are the strains chosen for the initial SnappNet and NetRax
analysis:

    ```
	C14_CSU_034_10640
	C33
	MN5
	MN12
	MN3
	C46
	BoviShield_MLV
	```

To restrict the large dataset `BHV1-plus-BHV5-outgroup-alignment.fasta` to just
these six strains, do the following:

First, check that you have the file
`BHV1-plus-BHV5-outgroup-alignment.fasta` saved in the `data/` directory.


Second, running the following command from the `data/` directory:

`grep -A1 -E '>C14_CSU_034_10640|>C33|>MN5|>MN12|>MN3|>C46|>BoviShield_MLV' BHV1-plus-BHV5-outgroup-alignment.fasta | grep -v -- "^--$" > BHV1-6-clinical-isolates.fasta`

This command saves the restricted dataset to the `data/` directory as
`BHV1-6-clinical-isolates.fasta`.


## Installing SnappNet

1. Download Beast2 for Linux from [this website](https://www.beast2.org/) to the `scripts/` directory. 

2. Extract the downloaded archive by running the following command from the `scripts/` directory:

`tar -xf BEAST_with_JRE.v2.6.6.Linux.tgz`

3. Beast requires Java, so check that you have java by running

`java -version`

If you don't have java, you will need to install it. Our results indicate that we do have it:

```
openjdk version "11.0.14" 2022-01-18
OpenJDK Runtime Environment (build 11.0.14+9-post-Debian-1deb10u1)
OpenJDK 64-Bit Server VM (build 11.0.14+9-post-Debian-1deb10u1, mixed mode, sharing)
```

4. Next we will follow the instructions for installing SnappNet from [this tutorial file](scripts/MySnappNet/workspace-Package-Beast/SnappNet/doc/TutoSnappNet). Navigate to your home directory `~/` by running the command `cd`. 

5. From the `~/` directory, run 

`mkdir .beast/2.6/SnappNet`

You should check your version of Beast2. We used version 2.6.6, so that in the above command we put the subdirectory 2.6. If you have a different version of beast, use 2.X where X is the major version of Beast. 

6. Next we need to copy the file `SnappNet.addon.zip` to the directory we just created in the previous step. To do this, navigate back to `scripts/` and run the command

`cp MySnappNet/workspace-Package-Beast/SnappNet/tmp/SnappNet.addon.zip ~/.beast/2.6/SnappNet/`


7. Unzip the file by running

`unzip ~/.beast/2.6/SnappNet/SnappNet.addon.zip`

8. Copy the snappnet version file to the beast directory by running the following command from `scripts/`

`cp MySnappNet/workspace-Package-Beast/SnappNet/release/package/version.xml ~/.beast/2.6/SnappNet/`


9. Navigate to `scripts/beast/` and run the command `./bin/beauti`. This will run Beauti.

10. In the GUI the pops up,

```
   File ---> Clear Class Paths
   File ---> Exit
```


10. Run Beauti again (by repeating the previous step). The in the GUI check the following:

```
File ----> Template 
You should now see SnappNetTemplate

File ----> Manage Packages
You should now see SnappNet
```


11. Read [A Rough Guide to SnappNet](scripts/MySnappNet/workspace-Package-Beast/SnappNet/doc/SnappNet.pdf). (This link will work if you've downloaded SnappNet to the correct directory in the previous steps.


## Running SnappNet with Our Data


### Convert to Nexus Format

We'll need to convert our file `BHV1-6-clinical-isolates.fasta` Nexus format. We will use seqmagick.
1. Install seqmagick by running the following command from `scripts/`

`git clone https://github.com/fhcrc/seqmagick.git`

2. Run

`pip install seqmagick`

3. Do the converstion by running the following command from `data/`,=

`seqmagick convert --output-format nexus --alphabet dna BHV1-6-clinical-isolates.fasta BHV1-6-clinical-isolates.nex`


### Prepare xml file with Beauti

1. Run Beauti. To do this, navigate to `scripts/beast/` and run the command `./bin/beauti`. This will run Beauti.

2. File --> Template -- choose snappnet

3. File --> Add Alignment -- choose the file `data/BHV1-6-clinical-isolates`

4.  Simply go to File --> Save as -- and save the file as `BHV1-6-clinical-isolates.xml` in the `data/` directory.

### Run the SnapNet thing with beast

1. Navigate to `scripts/beast/bin/` and run `./beast` (you should pipe the standard output to a file because this produces millions of lines of output)

2. Select the file `data/BHV1-6--clinical-isolates`

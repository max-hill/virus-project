# Analysis of Herpesvirus DNA Data

with "explicit" network methods to look for evidence of reticulation.  
data source: files from Aaron in [box](https://uwmadison.app.box.com/folder/147319895420)

cd scfirst focus: BHV1.1 because
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
- In Part 3, we will use SnappNet or NetRax to obtain a phylogenetic network for
  a subset of six taxa taken from the BHV1 dataset.
- In Part 4, we attempt to use TriLoNet to obtain a phylogenetic network.

## Part 0. Obtaining the Data

### Downloading the Full BHV1 Dataset
We used the datafile [BHV1\_plus\ BHV5\

(7mb fasta file). Download this file to the `data/` directory. Rename the file to 

`BHV1-plus-BHV5-outgroup-alignment.fasta`

Throughout this document, this file will be referred to as the *BHV1 Dataset*.
This is a fasta file without any linebreaks within the sequences. We will work
with this full BHV1 dataset directly as well as with subsets of it.

### Subsetting the Dataset (Required for SnappNet and NetRAX)
The following are the strains chosen for SnappNet and NetRax analysis:

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

The output is the file `data/BHV1-6-clinical-isolates.fasta` which consists of
the six clinical isolates and one vaccine strain listed above.

If this does not work, an alternative method is to use the R code in [this
file](scripts/Rutilies.Rmd).



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
1. Download IQ-TREE [here](http://www.iqtree.org/) and save it to the `scripts/`
directory. 

Note: we used the COVID-19 release 2.1.3 (April 21, 2021) for Linux.
Instructions for downloading and using IQ-TREE on other operating systems can be
found [here](http://www.iqtree.org/doc/Quickstart) (general instructions) and
[here](https://github.com/UWMadison-computingtools-2020/fp-group-6/blob/master/stepsinstructions.md#install-iq-tree)
(for Mac OS).

2. Extract IQ-TREE: from `scripts/`, run

`tar -xf iqtree-2.1.3-Linux.tar.gz`

3. Give the file `iqtree2` executable permission: from `scripts/`, run

`sudo chmod a+rx iqtree-2.1.3-Linux/bin/iqtree2`

4. Copy the file `iqtree2` to `/usr/local/bin` with the command:

`sudo cp iqtree-2.1.3-Linux/bin/iqtree2 /usr/local/bin`

We are now able to run IQ-TREE by running the command `iqtree2` in any
directory.

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

## Part 3A. SnappNet

### Installing SnappNet

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

10. In the GUI that pops up,

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


### Running SnappNet with Our Data


#### Convert to Nexus Format

We'll need to convert our file `BHV1-6-clinical-isolates.fasta` Nexus format. We will use seqmagick.
1. Install seqmagick by running the following command from `scripts/`

`git clone https://github.com/fhcrc/seqmagick.git`

2. Run

`pip install seqmagick`

3. Do the converstion by running the following command from `data/`,=

`seqmagick convert --output-format nexus --alphabet dna BHV1-6-clinical-isolates.fasta bhv1-6-clinical-isolates.nex`


#### Prepare xml file with Beauti
Beauti creates an .xml file that will be used as input by BEAST. Beauti provides
a GUI interface which allows the user to set parameters for the BEAST run.

1. To open Beauti, navigate to `scripts/beast/` and run the command
   `./bin/beauti`. This will run Beauti.

2. File --> Template -- choose SnappNetTemplate

3. File --> Add Alignment -- choose the file
   `data/BHV1-6-clinical-isolates.nex`. (Update: Or you could choose the file
   `data/BHV1-6-clinical-isolates.fasta`, it loads just fine.)

(It will take a moment to load the nexus file.)

4. Click the 'Model Parameters' tab, and then press the button 'Calc mutation
   rates'. Click the MCMC tab. Set "Chain Length" to 400000. Set "Store Every"
   to 10000. Under 'screenlog', 'tracelog', and 'specieslog', set "Log Every"
   to 10000. Actually, set 'screenlog' to 10000. Set the file name of tracelog
   to 'bhv1-6-clinical-isolates.trace.log' and set the file name of specieslog
   to 'bhv1-6-clinical-isolates.species.trees'

5. File --> Save as -- and save the file as `bhv1-6-clinical-isolates.xml` in
   the `data/` directory.

6. Close Beauti.

#### Run the SnapNet with Beast
This does not work. Produces likelihood calculation errors. Don't know how to
fix this issue.

1. Navigate to `scripts/beast/bin/` and run `./beast >screen-output.txt` (or
   `./beast -overwrite >screen-output.txt`). This will run BEAST, keep track of
   time, and pipe screen output to the file 'screen-output.txt'.

2. You will be promted to select an XML input file. Select the file we just
   created in Beauti: `data/bhv1-6-clinical-isolates.xml`.

3. Output: Unfortunately, we get output

```
At sample 10000
Likelihood incorrectly calculated: -38522.29225573676 != -43958.14044651774(5435.848190780984) Operator: snappNetProject.operators.ChangeGamma(ChangeGamma)
At sample 20000
Likelihood incorrectly calculated: -41492.85106857587 != -43973.661745134435(2480.810676558569) Operator: snappNetProject.operators.ChangeAllGamma(ChangeAllGamma)
At sample 30000
Likelihood incorrectly calculated: -41346.937857411365 != -43914.57251333519(2567.6346559238227) Operator: snappNetProject.operators.ChangeGamma(ChangeGamma)
At sample 40000
Likelihood incorrectly calculated: -39491.49321965366 != -43845.896530480604(4354.403310826943) Operator: snappNetProject.operators.ChangeAllGamma(ChangeAllGamma)
At sample 50000
Likelihood incorrectly calculated: -40549.58602114914 != -43918.23198724979(3368.645966100652) Operator: snappNetProject.operators.NodeUniform(NodeUniform:species)
At sample 60000
Likelihood incorrectly calculated: -40581.94689045801 != -43925.31812409661(3343.3712336385943) Operator: snappNetProject.operators.ChangeGamma(ChangeGamma)
At sample 70000
Likelihood incorrectly calculated: -36914.30121048426 != -43994.82679086487(7080.525580380614) Operator: snappNetProject.operators.ChangeGamma(ChangeGamma)
At sample 80000
Likelihood incorrectly calculated: -43723.29599458741 != -43781.59053304736(58.29453845995158) Operator: snappNetProject.operators.ChangeGamma(ChangeGamma)
```

## Part 3B. NetRAX
### Installing NetRAX
To install NetRAX, we will mostly follow the instruction in [the NetRAX
github](https://github.com/lutteropp/NetRAX) but with some modifications.

1. Install dependencies (same as github instructions)

`sudo apt-get install flex bison libgmp3-dev cmake doxygen libmpfrc++-dev libopenmpi-dev`

2. Build Instructions (this part is different from the github instructions). Run
   the following commands from the `scripts` directory:

```
git clone --recurse-submodules https://github.com/lutteropp/NetRAX.git
cd NetRAX
sed -i 's/master/main/' CMakeLists.txt.in
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON ..
make
```

(The modification was that we needed to change the word "master" to "main" in
the file `CMakeLists.txt.in` to avoid an error, and for some reason we also had
to run the command `make`).

3. Set binary path: Open the file `/scripts/NetRAX/netrax.py` and edit the line

`NETRAX_CORE_PATH = "/home/luttersh/NetRAX/bin/netrax"`

so that it gives the correct path for your NetRAX installation. Since we
installed NetRAX in the `scripts` directory, you should change the path to
something like

`NETRAX_CORE_PATH = "/.../virus-project/scripts/NetRAX/bin/netrax"`

where `...` is replaced by the appropriate absolute path on your computer. (I
used `NETRAX_CORE_PATH =
"/home/mutalisk/virus-project/scripts/NetRAX/bin/netrax"`)

4. Check that you can run Netrax: from `NetRAX/`, run

`python3 netrax.py --help`

or run

`./netrax --help`

from the directory `NetRAX/bin`

### Running NetRAX

The input of NetRAX consists of:
- fasta file
- newick format initial network(s)
- some other stuff which I don't understand yet

We aim to build a network consisting of the six clinical isolates plus the BHV5
outgroup. Before running NetRAX, we need to obtain the necessary inputs files:

1. Generate alignment: from the `data/` directory, run the command
 
`grep -A1 -E '>C14_CSU_034_10640|>C33|>MN5|>MN12|>MN3|>C46|>BoviShield_MLV' BHV1-plus-BHV5-outgroup-alignment.fasta | grep -v -- "^--$" > bhv1-6-clinical-isolates.fasta`

2. Install IQ-TREE using the instructions in Part 2.

3. Generate an initial tree with iqtree: from the `data/` directory, run the command

`iqtree2 -nt AUTO -s  bhv1-6-clinical-isolates.fasta`

4. Cleanup: create a folder for this analysis by running the following command from `data/`

`mkdir -p ../analysis/iqtree-output/6-clinical-isolates`

and then move the iqtree output to that directory by running

`mv bhv1-6-clinical-isolates.fasta.* ../analysis/iqtree-output/6-clinical-isolates/`

5. To test NetRAX, run it with the input files we just generated: from `/scripts/NetRAX` run

`mpiexec /home/mutalisk/virus-project/scripts/NetRAX/bin/netrax --start_network ~/virus-project/analysis/iqtree-output/6-clinical-isolates/bhv1-6-clinical-isolates.fasta.treefile --msa ~/virus-project/data/bhv1-6-clinical-isolates.fasta  --output ./bhv1-6-clinical-isolates-netrax-network.txt`

and it works, the get a network with no reticulations. 

6. Cleanup: delete the output files from our test run. From `/scripts/NetRAX/` run

`rm bhv1-6-clinical-isolates-netrax-network.*`


### Running NetRAX on partitioned data

Netrax requires partitioned MSA data in order to generate a network. That means
we need a parition file in addition to our full MSA. We will partition our MSA
into uniform blocks of length 500 base pairs. [Todo: we can also try using
partitionfinder.]

1. From `data`, run

`awk '{print length}' bhv1-6-clinical-isolates.fasta`

and the output tells us that the length of the gene sequences are 144552 each
(144551 base pairs + 1 newline character). We will make a dummy file

`dummy-partition-file-for-bhv1-6.txt`

```
DNA, gene_59=1-10000
DNA, gene_281=10001-20000
DNA, gene_281=20001-30000
DNA, gene_281=30001-40000
DNA, gene_281=40001-50000
DNA, gene_281=50001-60000
DNA, gene_281=60001-70000
DNA, gene_281=70001-80000
DNA, gene_281=80001-90000
DNA, gene_281=90001-100000
DNA, gene_281=100001-110000
DNA, gene_281=110001-120000
DNA, gene_281=120001-130000
DNA, gene_281=130001-140000
DNA, gene_281=140001-144551
```
Now we try runinng the following command from `scripts/NetRAX/`

`mpiexec /home/mutalisk/research/virus-project/scripts/NetRAX/bin/netrax --start_network ~/research/virus-project/analysis/iqtree-output/6-clinical-isolates-plus-bhv5-outgroup/bhv1-6-clinical-isolates-plus-bhv5-outgroup.fasta.treefile --msa ~/research/virus-project/data/bhv1-6-clinical-isolates-plus-bhv5-outgroup.fasta  --output ./bhv1-6-clinical-isolates-plus-bhv5-outgroup-netrax-network.txt --model ~/research/virus-project/data/dummy-partition-file-for-bhv1-6.txt`

and we get a network! But the outgroup is making visualization difficut. We should try running this on the six clinical isolates without the outgroup, We should also try using a partition with uniform lengths 500 or 1000. Also we can try using partitionfinder to get a network. 


## Part 4. Using TriLoNet

https://www.uea.ac.uk/groups-and-centres/computational-biology/software/trilonet

1. Download trilonet from
https://www.uea.ac.uk/documents/96135/3673668/TriLoNet3.zip/be80ee73-f11b-b2e0-565f-773efeaabdfc?t=1605188054795

and save it to the `scripts/` directory.

2. Extract trilonet into the `scripts/` directory

3. To run TriLoNet, we next want to follow instructons in the [TriLoNet
   Manual](scripts/TriLoNet3/TriLoNet/manual.pdf). Unfortunately, this returns a
   "file not found" error if we use an input file with any capital letters. To
   avoid this, navigate to the directory `scripts/TriLoNet3/TriLoNet/TriLoNet`
   and run

```
cp ../../../../data/BHV1-6-clinical-isolates.nex bhv1-6-clinical-isolates.nex
```

which will copy the data to the TriLoNet directory, using a name that has only
lowecase letters. Then we run the command


```
java -jar TriLoNet.jar bhv1-6-clinical-isolates.nex
```

Unfortunately, this return an error

```
Exception in thread "main" java.lang.NullPointerException
        at trilonet.ReadInput.readNEXUS(ReadInput.java:88)
        at trilonet.Main.main(Main.java:176)
```


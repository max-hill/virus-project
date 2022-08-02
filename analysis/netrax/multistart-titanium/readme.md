# README

Here we run netrax using the python wrapper with multiple start trees. The dataset is the titanium dataset with 8 taxa

The start trees are in `15genes_blksize10000_set1b_rooted.treefile` were generated using IQ-tree 


## Setting up
First we copied the dataset and partition files `titanium-dataset.fasta` and
`partition.txt` to this folder. Then we copied the file
`15genes_blksize10000_set1b_rooted.treefile`, which contains the start trees, to
this folder.

## Running Netrax

After completing setup, navigate to the directory `/scripts/NetRAX/` and run the following commands:

```
path="../../analysis/netrax/multistart-titanium"
start_networks="${path}/15genes_blksize10000_set1b_rooted.treefile"
partition_file="${path}/partition.txt"
msa_file="${path}/titanium-dataset.fasta"

python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 42 --likelihood_type average  --name multistart-titanium

```

The above command will expand to

```
python3 netrax.py --msa_path ../../analysis/netrax/multistart-titanium/titanium-dataset.fasta --partitions_path ../../analysis/netrax/multistart-titanium/partition.txt --start_networks ../../analysis/netrax/multistart-titanium/15genes_blksize10000_set1b_rooted.treefile --seed 42 --likelihood_type average --name multistart-titanium

```

## About the start trees
First, unrooted gene-trees were generated using IQ-Tree with the following options:

```

$iqtree2 -s $ALN_FILE -S $PARTITION_FILE -B $breps --prefix $prefix -wbt -T AUTO

```

These options are:

```
-S FILE|DIR: separate tree inference (not edge-linked)
-B: resample sites within partitions for bootstrapping (I didn't upload the bootstrap trees though)
-T AUTO: detect the best no. of CPU cores
-wbt: "turn on" writing bootstrap trees

```

Then the gene-trees which were rooted using TreeTime.

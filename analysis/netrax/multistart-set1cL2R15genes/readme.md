# Instructions for setting up this experiment
## Set up experiment directories
Here we set up the following four experiments simultaneously:

```
multistart-set1cL2R15genes
multistart-set1cR2L15genes
multistart-set1cL2R97genes
multistart-set1cR2L97genes
```

To do this, first make the directories

```
/analysis/netrax/multistart-set1cL2R15genes/
/analysis/netrax/multistart-set1cR2L15genes/
/analysis/netrax/multistart-set1cL2R97genes/
/analysis/netrax/multistart-set1cR2L97genes/
```

Then go to the `/data/` directory and run the command

```
grep -A1 -E '>C33$|>C46$|>Titanium_IBR_MLV_vaccine$|>Cooper$|>SP1777$|>B589$|>BHV5$|>216_II$|>K22$' BHV1-plus-BHV5-outgroup-alignment.fasta | grep -v -- "^--$" > set1c.fasta
```

Then copy the MSA file `set1c.fasta` and the partition files to the appropriate experiment directories. To do this, navigate to `/data/` and run

```
cp set1c.fasta ../analysis/netrax/multistart-set1cL2R15genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cR2L15genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cL2R97genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cR2L97genes/set1c.fasta

cp bhv/part/15genes_blksize10000.txt ../analysis/netrax/multistart-set1cL2R15genes/partition.txt
cp bhv/part/15genes_blksize10000_rev.txt ../analysis/netrax/multistart-set1cR2L15genes/partition.txt
cp bhv/part/97genes_blksize1500_rev.txt ../analysis/netrax/multistart-set1cR2L97genes/partition.txt
cp bhv/part/97genes_blksize1500.txt ../analysis/netrax/multistart-set1cL2R97genes/partition.txt
```


Then we copy the starting networks to the appropriate directories by navigating
to `/data/bhv/genetrees/ogrooted/` and running:

```
cp 97genes_blksize1500_set1c.treefile ../../../../analysis/netrax/multistart-set1cL2R97genes/starting-networks.treefile

cp 97genes_blksize1500_rev_set1c.treefile ../../../../analysis/netrax/multistart-set1cR2L97genes/starting-networks.treefile

cp 15genes_blksize10000_set1c.treefile ../../../../analysis/netrax/multistart-set1cL2R15genes/starting-networks.treefile

cp 15genes_blksize10000_rev_set1c.treefile ../../../../analysis/netrax/multistart-set1cR2L15genes/starting-networks.treefile

```

## Install NetRAX correctly

Navigate to `/scripts/` and run the following

```
sudo apt-get install flex bison libgmp3-dev cmake doxygen libmpfrc++-dev libopenmpi-dev
git clone --recurse-submodules https://github.com/lutteropp/NetRAX.git
cd NetRAX
git reset --hard 542d4c12
sed -i 's/master/main/' CMakeLists.txt.in
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON ..
make
```


Then in the file `NetRAX/netrax.py`, observe that the variable `NETRAX_CORE_PATH` is set to `"/home/luttersh/NetRAX/bin/netrax"`. Change this to the correct path on your machine. To do this, I ran the following command from `/scripts/NetRAX`

```
sed -i 's/\/home\/luttersh\/NetRAX\/bin\/netrax/\/home\/zergling\/.emacs.d\/virus-project\/scripts\/NetRAX\/bin\/netrax/' netrax.py
```

# Running the experiments
After completing the setup, navigate to the directory `/scripts/NetRAX/` and do the following: 

First, run exactly one of the following lines (corresponding to the desired experiment to run):

```
experiment_name="multistart-set1cL2R15genes"
experiment_name="multistart-set1cR2L15genes"
experiment_name="multistart-set1cL2R97genes"
experiment_name="multistart-set1cR2L97genes"
```

Then run all of the following:

```
path="../../analysis/netrax/${experiment_name}"
start_networks="${path}/starting-networks.treefile"
partition_file="${path}/partition.txt"
msa_file="${path}/set1c.fasta"

python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 42 --likelihood_type average  --name ${experiment_name}

```


When the run has finished, we copy the output from `scripts/NetRAX/` to the
experiment folder manually.


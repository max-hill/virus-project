# NOTE

This experiment was never performed. I killed this process early on in the run
when it became clear that the 97 gene case was going to take a long time. For 97
genese I only kept the L2R case running (but will probably kill that run early
as well). 2023-06-04 max hill


# Instructions for setting up this experiment

## Install NetRAX correctly
In order to get NetRAX to work correctly, it is important to install NetRAX exactly in the manner described in this section.

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


## Set up the four experiment directories
In this section explain how we set up the following four simultaneously:

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

Then, to generate the correct data files (a multiple species alignment with the appropriate taxa in fasta format), go to the `/data/` directory and run the command

```
grep -A1 -E '>C33$|>C46$|>Titanium_IBR_MLV_vaccine$|>Cooper$|>SP1777$|>B589$|>BHV5$|>216_II$|>K22$' BHV1-plus-BHV5-outgroup-alignment.fasta | grep -v -- "^--$" > set1c.fasta
```

Then copy the MSA file `set1c.fasta` to the experiment directories by running the following commands from the `/data/` directory:

```
cp set1c.fasta ../analysis/netrax/multistart-set1cL2R15genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cR2L15genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cL2R97genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cR2L97genes/set1c.fasta

```

Next, the partition files have already been created and are located in the directory `/data/bhv/part/`. We will copy these to the experiment directories as well, by running the following commands again from the `/data/` directory:

```
cp bhv/part/15genes_blksize10000.txt ../analysis/netrax/multistart-set1cL2R15genes/partition.txt
cp bhv/part/15genes_blksize10000_rev.txt ../analysis/netrax/multistart-set1cR2L15genes/partition.txt
cp bhv/part/97genes_blksize1500_rev.txt ../analysis/netrax/multistart-set1cR2L97genes/partition.txt
cp bhv/part/97genes_blksize1500.txt ../analysis/netrax/multistart-set1cL2R97genes/partition.txt
```


The starting networks have also already been created. We copy the starting networks to the appropriate directories by navigating
to `/data/bhv/genetrees/ogrooted/` and running:

```
cp 97genes_blksize1500_set1c.treefile ../../../../analysis/netrax/multistart-set1cL2R97genes/starting-networks.treefile

cp 97genes_blksize1500_rev_set1c.treefile ../../../../analysis/netrax/multistart-set1cR2L97genes/starting-networks.treefile

cp 15genes_blksize10000_set1c.treefile ../../../../analysis/netrax/multistart-set1cL2R15genes/starting-networks.treefile

cp 15genes_blksize10000_rev_set1c.treefile ../../../../analysis/netrax/multistart-set1cR2L15genes/starting-networks.treefile

```

# Instructions for running the experiment
After completing the setup and installing NetRAX as described in the previous sections, we can now run the experiments by following the instructions in the two subsections below. The experiments for 15 genes (`multistart-set1cR2L15genes` and `multistart-set1cR2L15genes`) and 97 genes (`multistart-set1cR2L97genes` and `multistart-set1cR2L97genes`) are run slightly differently. Since we expected the experiments with 97 genes to take a long time to finish, we wanted to run them on a remote server that we could log off of, so we wrote scripts to handle that.

## Instructions for running the Experiments with 97 genes 
This subsection describes how to run the experiements `multistart-set1cR2L97genes` and `multistart-set1cR2L97genes`. To do this, first copy the script `runL2R97.sh` (or `runR2L97.sh`) to the directory `/scripts/NetRAX`. Then run the script from that directory. Detailed instructions for running those scripts (especially on a remote server) can be found as comments in the scripts themselves. 

When the run has finished, we copy the output from `scripts/NetRAX/` to the experiment folder manually.

## Instructions for running the Experiments with 15 genes
This subsection describes how to run the experiements `multistart-set1cR2L15genes` and `multistart-set1cR2L15genes`.

After completing the setup and installing NetRAX, navigate to the directory `/scripts/NetRAX/` and do the following: 

First, run exactly one of the following lines, depending on whether you wish to run the L2R or the R2L experiment:

```
experiment_name="multistart-set1cL2R15genes"
experiment_name="multistart-set1cR2L15genes"
```

Then run all of the following:

```
path="../../analysis/netrax/${experiment_name}"
start_networks="${path}/starting-networks.treefile"
partition_file="${path}/partition.txt"
msa_file="${path}/set1c.fasta"

python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 42 --likelihood_type average  --name ${experiment_name}

```

When the run has finished, we copy the output from `scripts/NetRAX/` to the experiment folder manually.


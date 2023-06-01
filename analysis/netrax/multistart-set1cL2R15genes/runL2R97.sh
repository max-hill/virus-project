#!bin/bash

# This is for running the experiment named below. Run this script from the directory scripts/NetRAX/

experiment_name="multistart-set1cL2R97genes"
path="../../analysis/netrax/${experiment_name}"
start_networks="${path}/starting-networks.treefile"
partition_file="${path}/partition.txt"
msa_file="${path}/set1c.fasta"

echo "Command Used:">multistart-set1cL2R97genes-script-output.txt
echo "">>multistart-set1cL2R97genes-script-output.txt
echo "python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 42 --likelihood_type average  --name ${experiment_name}" >>multistart-set1cL2R97genes-script-output.txt
echo "">>multistart-set1cL2R97genes-script-output.txt
echo "">>multistart-set1cL2R97genes-script-output.txt

time python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 42 --likelihood_type average  --name ${experiment_name} >> multistart-set1cL2R97genes-script-output.txt



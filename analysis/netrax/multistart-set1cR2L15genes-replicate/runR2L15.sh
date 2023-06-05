#!bin/bash

# 2023-06-04 max hill

# This is for running the experiment named below. Run this script from the directory scripts/NetRAX/

# To run this on a server, run the commands (from scripts/NetRAX):
# chmod u+rwx run R2L15.sh
# nohup bash runR2L15.sh &
# echo $! > PID_for_R2L15.txt

# To kill this process run the command 'kill PID' where PID is the number in PID_for_R2L15.txt.

experiment_name="multistart-set1cR2L15genes-replicate"
path="../../analysis/netrax/${experiment_name}"
start_networks="${path}/starting-networks.treefile"
partition_file="${path}/partition.txt"
msa_file="${path}/set1c.fasta"

echo "Command Used:">multistart-set1cR2L15genes-replicate-script-output.txt
echo "">>multistart-set1cR2L15genes-replicate-script-output.txt
echo "python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 42 --likelihood_type average  --name ${experiment_name}" >>multistart-set1cR2L15genes-replicate-script-output.txt
echo "">>multistart-set1cR2L15genes-replicate-script-output.txt
echo "">>multistart-set1cR2L15genes-replicate-script-output.txt

time python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 42 --likelihood_type average  --name ${experiment_name} >> multistart-set1cR2L15genes-replicate-script-output.txt



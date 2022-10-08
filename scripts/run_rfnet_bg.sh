#!/bin/bash
#SBATCH -o run_rfnet_bg.log
#SBATCH -J rfnet
#SBATCH -t 2:00:00
#SBATCH -p debug
#SBATCH -n 1 # allocate 1 CPU
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1 # no. of threads per process

cd ../../RF-Net-2.0.4

input=../virusanalysis/genetrees/15genes_blksize10000_set1c_rooted/15genes_blksize10000_set1c_rooted.treefile
output=../virusanalysis/rfnet_analysis/15genes_blksize10000_set1c_rooted_r10.newick

# redirect stderr to `error.txt`
java -jar RF-Net-2.0.4.jar -i $input -o $output -r 10 -e 2> rfnet_error.txt

mv rfnet_error.txt ../virusanalysis/rfnet_analysis
#mv full.newick ../virusanalysis/rfnet_analysis
#mv full.newick.embeddings.tre ../virusanalysis/rfnet_analysis

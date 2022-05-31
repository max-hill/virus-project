#!/bin/bash
#_______________________________________________________________________________
#
# run-netrax-experiment --- Command line script to generate explicit networks
#                           with netrax on subsets of BHV1 herpesvirus dataset.
#_______________________________________________________________________________

# Authors: Max Hill, Ben Teo, Shuqi Yu
# (Last updated 2022-05-30)
#
# INSTRUCTIONS: This script is intended to be run from the directory
# virus-project/scripts/ using the command line. To run this script, use a
# command of the form 'bash run-netrax-experiment.sh experiment-label
# partition-block-size taxa-1 taxa-2 ... taxa-2' where experiment-label is a
# string representing the name of your experiment (e.g. "X"),
# partition-block-size is a number (e.g. 500), and taxa-1,...,taxa-n are names
# of virus taxa in the BHV1 dataset (if you misspell a name, the program will
# prompt you with a list of acceptable names.
#
# Example usage: "bash run-netrax-experiment X 500 MN1 MN13 BHV5"


### ---Begin script---

################################################################################

# Check if script is run from the correct directory
dir=$(pwd | awk -F'/' '{print $NF}')
if [[ $dir != scripts ]]
then
    echo "ERROR: this script must be run from the directory 'virus-project/scripts/'."
    exit
fi

# Name the input variable. 
experiment_label="$1"
partition_block_size="$2"

# Check that all taxon names are entered correctly
echo "Analyzing the following set of taxa:"

for ((i=3; i<=$#; i++))
do
    taxon="${!i}"
    if [[ $taxon != +(BHV5|MN1|MN2|MN3|MN4|MN5|MN6|MN7|MN8|MN9|MN10|MN11|MN12|MN13|MN14|MN15|PA1|PA2|PA3|C14_CSU_034_10640|C18|C26|C28_55771|C29|C33|C35_1839_9847|C36_876_459|C42|C43|C44|C45|C46|C47|Nasalgen_IP_MLV_vaccine|TSV-2_Nasal_MLV_vaccine|BoviShield_Gold_FP5_MLV_vaccine|BovSh_IBR_MLV_vaccine|Vista_IBR_MLV_vaccine|Pyramid_IBR_MLV_vaccine|Express1_IBR_MLV_vaccine|Titanium_IBR_MLV_vaccine|Arsenal_IBR_MLV_vaccine|VR188_Los_Angeles|NVSL_challenge_97_11|216_II|SP1777|SM023|K22|B589|Cooper) ]]
    then
        echo "ERROR: TAXON ENTERED INCORRECTLY: $taxon"
        echo -e "Taxon names must be entered exactly as they are written in the following list:
\n216_II\nArsenal_IBR_MLV_vaccine\nB589\nBHV5\nBoviShield_Gold_FP5_MLV_vaccine\nBovSh_IBR_MLV_vaccine\nC14_CSU_034_10640\nC18\nC26\nC28_55771\nC29\nC33\nC35_1839_9847\nC36_876_459\nC42\nC43\nC44\nC45\nC46\nC47\nCooper\nExpress1_IBR_MLV_vaccine\nK22\nMN1\nMN2\nMN3\nMN4\nMN5\nMN6\nMN7\nMN8\nMN9\nMN10\nMN11\nMN12\nMN13\nMN14\nMN15\nNasalgen_IP_MLV_vaccine\nNVSL_challenge_97_11\nPA1\nPA2\nPA3\nPyramid_IBR_MLV_vaccine\nSP1777\nSM023\nTitanium_IBR_MLV_vaccine\nTSV-2_Nasal_MLV_vaccine\nVista_IBR_MLV_vaccine\nVR188_Los_Angeles"
        exit
    else
        echo "$taxon"
    fi
done

# Make the directory for experiment
if [ -d "../analysis/netrax/experiment-$experiment_label" ]
then
    echo "ERROR: the directory ../analysis/netrax/experiment-$experiment_label"
    echo "already exists. Choose a different experiment name."
    exit
else
    echo "Making directory for the experiment:"
    echo "~/virus-project/analysis/netrax/experiment-$experiment_label"
    mkdir ../analysis/netrax/experiment-$experiment_label
fi

# Make partition file
bash generate-partition-file.sh 144551 $partition_block_size >../analysis/netrax/experiment-$experiment_label/partition.txt
echo "Using a partition consisting of blocks of length $partition_block_size".
echo "Partition file written to the directory"
echo "../analysis/netrax/experiment-$experiment_label/partition.txt"


# Make alignment file
for ((i=3; i<=$#; i++))
do
    taxon="${!i}"
    grep -A1 -E ">$taxon$" ../data/BHV1-plus-BHV5-outgroup-alignment.fasta | grep -v -- "^--$" >> ../analysis/netrax/experiment-$experiment_label/experiment-$experiment_label-dataset.fasta
done

# Use iqtree to generate an ML tree
cd ../analysis/netrax/experiment-$experiment_label/
iqtree2 -nt AUTO -s experiment-$experiment_label-dataset.fasta -pre experiment-$experiment_label


# Run netrax using the partition, alignment, and ML tree as input
cd ../../../scripts/NetRAX/bin/
experiment_path="../../../analysis/netrax/experiment-$experiment_label"


START=$(date +%s)
mpiexec ./netrax --name experiment$experiment_label --msa $experiment_path/experiment-$experiment_label-dataset.fasta --model $experiment_path/partition.txt --average_displayed_tree_variant --start_network $experiment_path/experiment-$experiment_label.treefile --output $experiment_path/experiment-$experiment_label-netrax-output.txt --seed 42
END=$(date +%s)
netrax_runtime=$(( $END - $START ))

# Create a readme with documentation and copy this script to the experiment folder
cp "../../run-netrax-experiment.sh" "$experiment_path/run-netrax-experiment.sh"

echo "# Readme" >$experiment_path/readme.md

echo -e "\nThis folder and its contents were produced with the shell script
\`run-netrax-experiment.txt\`. In particular, the script was run from the
directory \`scripts/\` and the command used was \`bash run-netrax-experiment.sh
$@\`" >>$experiment_path/readme.md

echo -e "\nThe best inferred network by netrax can be found in the file
\`experiment-$experiment_label-netrax-output.txt\`" >>$experiment_path/readme.md

echo -e "\nThe best inferred tree by iqtree can be found in the file
\`experiment-$experiment_label.treefile\`" >>$experiment_path/readme.md

echo -e "\nNetrax runtime: $netrax_runtime seconds" >>$experiment_path/readme.md

echo -e "\nLast updated: $(date)"  >>$experiment_path/readme.md

# Tell the user where to find the output
echo -e "\nAll output written to the folder"
echo -e "\n../analysis/netrax/experiment-$experiment_label/"
echo -e "\nSee readme.md in that folder for further details."

# Return to scripts/ directory
cd ../../

################################################################################

### ---Script ends here---

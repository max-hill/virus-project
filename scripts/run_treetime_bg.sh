#!/bin/bash
#SBATCH -o run_treetime_bg.log
#SBATCH -J treetime
#SBATCH -t 2:00:00
#SBATCH -p debug
#SBATCH -n 1 # allocate 1 CPU
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1


treetime=../treetime/bin/treetime
trees=15genes_blksize10000_set1c/15genes_blksize10000_set1c.treefile
dates=../sampling_times/sampling_times.csv

# <Rooting gene/locus trees>
genenum=1
numgenes=15
while read -r line
do
    echo "$line" > tree.newick # unrooted gene tree
    if [ $genenum -ne $numgenes ]; then
        seqlen=10000
    else
        seqlen=44551 # the last gene is shorter
    fi
    # root gene tree
    $treetime clock --tree tree.newick --dates $dates --sequence-len $seqlen --reroot least-squares --outdir clock_results
    # append rooted gene tree to file
    cat clock_results/rerooted.newick >> 15genes_blksize10000_set1c_rooted.treefile
    rm tree.newick
    let genenum++
done < "$trees" 2> treetime_error.txt

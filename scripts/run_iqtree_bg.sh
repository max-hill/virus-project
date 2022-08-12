#!/bin/bash
#SBATCH -o run_iqtree_bg.log
#SBATCH -J iqtree
#SBATCH -t 2:00:00
#SBATCH -p debug
#SBATCH -n 1 # allocate 1 CPU
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1

#ALN_FILE=../msa/BHV1-plus-BHV5-outgroup-alignment.fasta
ALN_FILE=../msa/set1c.fasta

PARTITION_FILE=../part/15genes_blksize10000.txt

#iqtree=../../iqtree-1.6.12-Linux/bin/iqtree
iqtree2=../../iqtree-2.2.0-Linux/bin/iqtree2

prefix=15genes_blksize10000_set1c
#seed=722261
breps=1000

# <Inferring gene/locus trees>
# (http://www.iqtree.org/doc/Concordance-Factor#inferring-genelocus-trees)
# --prefix loci: all output files will be loci.*
# -T AUTO: detect the best no. of CPU cores
# -S FILE|DIR: separate tree inference (not edge-linked)
# -B: resample sites within partitions
# -wbt: turn on writing bootstrap trees to .ufboot
$iqtree2 -s $ALN_FILE -S $PARTITION_FILE -B $breps -pre $prefix -wbt -T AUTO -redo 2> iqtree_error.txt

# <Inferring species tree>
#$iqtree2 -s $ALN_FILE -p $PARTITION_FILE --prefix $prefix -T AUTO -redo 2> iqtree_error.txt

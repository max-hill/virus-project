#!/bin/bash
#_______________________________________________________________________________
#
# generate-partition-file.sh --- script to make a partition file for an MSA
#_______________________________________________________________________________
# Author: Max Hill
# (Last updated: 2022-03-12)
#
# INSTRUCTIONS: This file is intended to be run from the command line. To run
# this script, navigate to the scripts/ directory and run the command 'bash
# generate-partition-file.sh A B > output-file.txt' where A and B are positive
# integers and output-file.txt is the desired name of the partition file. The
# inputs A and B are:
#   A = total alignment length (i.e. the length of the DNA sequences in your MSA
#     that you seek to partition).
#   B = desired length of partition blocks. Must have B<A.
#
# OUTPUT: An MSA partition file is piped to standard output. If the block length
# B does not evenly divide the alignment length A, the last partition will
# consist of fewer than B base pairs.
#
# DEPENDENCIES: This script requires Steel Bank Common Lisp (SBCL)
# (http://www.sbcl.org/)
#
# EXAMPLE USAGE:
# 'bash generate-partition-file.sh 144551 500 >../data/partition-A144551-B500.txt'
# Needs to be run from the directory virus-analysis/scripts/

### ---Begin script---

################################################################################

# Name the input variables
A="$1"
B="$2"

# Call the lisp functions
sbcl --noinform --eval '
(progn 
       (load "generate-partition-file-aux.lisp")
       (generate-partition-file *A* *B*)
       (quit)) 
' $A $B 2>/dev/null

# Note that SBCL ('steel bank common lisp') is an implmentation of common lisp.
# When loading files, SBCL outputs style warnings by default. Since this might
# confuse users not familiar with lisp, such behavior is surpressed by piping
# standard error into the black hole of /dev/null.

################################################################################

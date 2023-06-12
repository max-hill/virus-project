#!/bin/bash
#SBATCH -o run_mc3_test.log
#SBATCH -t 8-00:00:00 # run job for <= 8 days
#SBATCH -p long # use long partition
#SBATCH -n 1 # allocate 1 CPU
#SBATCH --mem-per-cpu=2600M # RAM per CPU/processor
#SBATCH --cpus-per-task=9 # no. of threads
export BEAST_PACKAGE_PATH="/workspace/bteo/.beast/2.7"

# -overwrite: overwrite log files
# -resume: append to log file
# -validate: parse XML, but do not run
# -threads: no. of threads to use (default: 1, no. of cores: -1)
# -instances: divide site patterns amongst no. of threads (use with -threads option)
# -packagedir: set user package directory instead of using default
bash /workspace/bteo/beast/bin/beast -seed 333 -overwrite -threads -1 mc3_test.xml 2> mc3_test.txt

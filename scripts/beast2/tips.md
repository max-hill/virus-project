## General *beast*

### Installation (local machine)
- [BEAST-2.6.6 for Windows (w/ Java)](http://www.beast2.org/download-windows-with-jre/)


### Installing new packages (e.g. SNAPP)
- Do `package manager -add SNAPP`. The `package manager` executable is located
within the BEAST installation (e.g. `beast/bin/packagemanager`). Reading it is
helpful for understanding the BEAST/JAVA environment variables and how they are
used.
- After a new package is installed, this change should be reflected (i.e.
`<new package name>.addon.jar` should be appended to `package.path`) in the
`beauti.properties` file located in the package directory.

### Multithreading

### Adjusting transition kernels

### Debugging
- [Beast-users google-group](https://groups.google.com/g/beast-users)
---

## Running *snapp* in the foreground
```
#!/bin/bash

# You may need to set this environment variable. The `beast` executable will
# look here for additional packages you installed.
export BEAST_PACKAGE_PATH="/workspace/bteo/.beast/2.6"

# <path to executable> <options> <path to BEAUti-generated .xml script>
# `>` at the end is for piping standard output to a log file
/workspace/bteo/beast/bin/beast -seed 123 -overwrite snapp_default.xml >
snapp_default_fg.log
```

Suppose you saved the above script as `run_snapp_default.sh`, then you would
do
```
$ ./run_snapp_default.sh
```

Note that `-overwrite` means that the entire analysis is restarted and the
log/state files generated from a previous analysis will be overwritten. If
`-resume` is used instead
```
/workspace/bteo/beast/bin/beast -seed 123 -resume snapp_default.xml >
snapp_default_fg.log
```
then the analysis will start from the end state of the previous analysis,
which is read from a `.xml.state` file. In this case, you might want to append
to `snapp_default_fg.log ` using `>>` rather than overwriting it using `>`.

---

## Submitting a *snapp* job on the cluster 
Suppose you want to submit a job to a HPC cluster that schedules jobs using
SLURM, then you might write a script like this
```
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bteo@wisc.edu
#SBATCH -o run_snapp_bg.log
#SBATCH -J snapp
#SBATCH -t 02:30:00
#SBATCH -p short # send this job to the debug partition
#SBATCH -n 1 # allocate 1 CPU
#SBATCH --mem-per-cpu=150M
#SBATCH --cpus-per-task=1 # no. of threads per process

export BEAST_PACKAGE_PATH="/workspace/bteo/.beast/2.6"
bash /workspace/bteo/beast/bin/beast -seed 29 -resume snapp_default.xml
```
To submit the job and leave it to run in the background, you would do
```
$ sbatch run_snapp_default.sh
```

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
#SBATCH -p short # send this job to the `short` partition
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

---

## Submitting a *snappnet* job on the cluster
Instead of running the `beast` executable, we run the
[`SnappNetProjectToRun.jar`](https://github.com/rabier/MySnappNet/tree/master/workspace-Package-Beast/SnappNet/deliverable)
executable.
```
#!/bin/bash
#SBATCH -o run_snappnet_bg.log
#SBATCH -J snappnet
#SBATCH -t 02:30:00
#SBATCH -p short
#SBATCH -n 1 # allocate 1 CPU
#SBATCH --mem-per-cpu=150M
#SBATCH --cpus-per-task=1 # no. of threads per process

# run MCMC analysis with SnappNet
java -jar SnappNetProjectToRun.jar JDD1.xml 2> error.txt
```

The following line in the header redirects `stdout` to `run_snappnet_bg.log`. 
```
#SBATCH -o run_snappnet_bg.log
```

`stderr` may contain useful information (e.g. likelihood calculation errors), so
we pipe it (e.g. using `2>`) to `error.txt`. Displayed below is what `error.txt`
contained after 35000 samples.
```
File: JDD1.xml seed: 127 threads: 1
Or3: 5682 3
Jap: 5682 3
Aus: 5682 3
Ind: 5682 3
Aro: 5682 3
Or1: 5682 3
Operator                                                                         Tuning    #accept    #reject      Pr(m)  Pr(acc|m)
ScaleOperator(divrRateScale:species)                                            0.23854        393        436    0.02326    0.47407 Try setting scaleFactor to about 0.057
ScaleOperator(turnOverScale:species)                                            0.18339        339        411    0.02326    0.45200 Try setting scaleFactor to about 0.038
snappNetProject.operators.InheritanceProbUniform(gammaProbUniform:species)            -          3        796    0.02326    0.00375
snappNetProject.operators.InheritanceProbRndWalk(gammaProbRndWalk:species)            -         22        800    0.02326    0.02676
snappNetProject.operators.OriginMultiplier(originMultiplier:species)            1.00000        263        161    0.01163    0.62028
snappNetProject.operators.AddReticulation(addReticulation:species)                    -          2        782    0.02326    0.00255
snappNetProject.operators.DeleteReticulation(deleteReticulation:species)              -          1        803    0.02326    0.00124
snappNetProject.operators.NetworkMultiplier(networkMultiplier:species)                -        121        286    0.01163    0.29730
snappNetProject.operators.FlipReticulation(flipReticulation:species)                  -          0        827    0.02326    0.00000
snappNetProject.operators.RelocateBranch(relocateBranch:species)                      -         16        816    0.02326    0.01923
snappNetProject.operators.NodeSlider(nodeSlider:species)                              -        378        440    0.02326    0.46210
snappNetProject.operators.NodeUniform(NodeUniform:species)                            -         58        776    0.02326    0.06954
snappNetProject.operators.RelocateBranchNarrow(relocateBranchNarrow:species)          -         14        777    0.02326    0.01770
snappNetProject.operators.ChangeUAndV(ChangeUAndV)                                    -         28        766    0.02326    0.03526
snappNetProject.operators.ChangeGamma(ChangeGamma)                              0.79253       4281       7844    0.34884    0.35307
snappNetProject.operators.ChangeAllGamma(ChangeAllGamma)                              -       1368      10866    0.34884    0.11182

     Tuning: The value of the operator's tuning parameter, or '-' if the operator can't be optimized.
    #accept: The total number of times a proposal by this operator has been accepted.
    #reject: The total number of times a proposal by this operator has been rejected.
      Pr(m): The probability this operator is chosen in a step of the MCMC (i.e. the normalized weight).
  Pr(acc|m): The acceptance probability (#accept as a fraction of the total proposals for this operator).

```

For example, we see
- initial seed (127), no. of threads (1)
- identifiers for the 6 taxa
- no. of sites (5682), [state count](https://www.beast2.org/xml/beast.evolution.alignment.AscertainedAlignment.html) (3)
- [operator report](http://www.beast2.org/2021/07/20/operator-tuning.html)

To resume an analysis, we run a different executable.
```
java -jar SnappNetProjectToRunResume.jar -resume JDD1.xml
```
If the `-resume` flag is used, the tree log (i.e. `.species.trees`) and
trace log (i.e. `.trace.log`) are appended to. For example, if you had (1) run
the analysis for 35000 samples previously and resumed for another 1000 samples,
and (2) logged information every 1000 samples, then both log files would
contain information for each of the 36 time points (i.e. 0, 1000, ..., 36000). 

However, `run_snappnet_bg.log` gets overwritten
(e.g. `tail -n 2 run_snappnet_bg.log`).
```
Sample      posterior ESS(posterior)     likelihood          prior
 35000    -33552.8232              N    -33384.4882      -168.3350 --
 36000    -33552.1381         2.0       -33381.0335      -171.1045 --
```
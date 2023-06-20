# Instructions for setting up this experiment

## Install NetRAX correctly
In order to get NetRAX to work correctly, it is important to install NetRAX exactly in the manner described in this section.

Navigate to `/scripts/` and run the following

```
sudo apt-get install flex bison libgmp3-dev cmake doxygen libmpfrc++-dev libopenmpi-dev
git clone --recurse-submodules https://github.com/lutteropp/NetRAX.git
cd NetRAX
git reset --hard 542d4c12
sed -i 's/master/main/' CMakeLists.txt.in
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON ..
make
```

Then in the file `NetRAX/netrax.py`, observe that the variable `NETRAX_CORE_PATH` is set to `"/home/luttersh/NetRAX/bin/netrax"`. Change this to the correct path on your machine. To do this, I ran the following command from `/scripts/NetRAX`

```
sed -i 's/\/home\/luttersh\/NetRAX\/bin\/netrax/\/home\/zergling\/.emacs.d\/virus-project\/scripts\/NetRAX\/bin\/netrax/' netrax.py
```


## Set up the four experiment directories
In this section explain how we set up the following four simultaneously:

```
multistart-set1cL2R15genes
multistart-set1cR2L15genes
multistart-set1cL2R97genes
multistart-set1cR2L97genes
```

To do this, first make the directories

```
/analysis/netrax/multistart-set1cL2R15genes/
/analysis/netrax/multistart-set1cR2L15genes/
/analysis/netrax/multistart-set1cL2R97genes/
/analysis/netrax/multistart-set1cR2L97genes/
```

Then, to generate the correct data files (a multiple species alignment with the appropriate taxa in fasta format), go to the `/data/` directory and run the command

```
grep -A1 -E '>C33$|>C46$|>Titanium_IBR_MLV_vaccine$|>Cooper$|>SP1777$|>B589$|>BHV5$|>216_II$|>K22$' BHV1-plus-BHV5-outgroup-alignment.fasta | grep -v -- "^--$" > set1c.fasta
```

Then copy the MSA file `set1c.fasta` to the experiment directories by running the following commands from the `/data/` directory:

```
cp set1c.fasta ../analysis/netrax/multistart-set1cL2R15genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cR2L15genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cL2R97genes/set1c.fasta
cp set1c.fasta ../analysis/netrax/multistart-set1cR2L97genes/set1c.fasta

```

Next, the partition files have already been created and are located in the directory `/data/bhv/part/`. We will copy these to the experiment directories as well, by running the following commands again from the `/data/` directory:

```
cp bhv/part/15genes_blksize10000.txt ../analysis/netrax/multistart-set1cL2R15genes/partition.txt
cp bhv/part/15genes_blksize10000_rev.txt ../analysis/netrax/multistart-set1cR2L15genes/partition.txt
cp bhv/part/97genes_blksize1500_rev.txt ../analysis/netrax/multistart-set1cR2L97genes/partition.txt
cp bhv/part/97genes_blksize1500.txt ../analysis/netrax/multistart-set1cL2R97genes/partition.txt
```


The starting networks have also already been created. We copy the starting networks to the appropriate directories by navigating
to `/data/bhv/genetrees/ogrooted/` and running:

```
cp 97genes_blksize1500_set1c.treefile ../../../../analysis/netrax/multistart-set1cL2R97genes/starting-networks.treefile

cp 97genes_blksize1500_rev_set1c.treefile ../../../../analysis/netrax/multistart-set1cR2L97genes/starting-networks.treefile

cp 15genes_blksize10000_set1c.treefile ../../../../analysis/netrax/multistart-set1cL2R15genes/starting-networks.treefile

cp 15genes_blksize10000_rev_set1c.treefile ../../../../analysis/netrax/multistart-set1cR2L15genes/starting-networks.treefile

```

# Instructions for running the experiment
After completing the setup and installing NetRAX as described in the previous sections, we can now run the experiments by following the instructions in the two subsections below. The experiments for 15 genes (`multistart-set1cR2L15genes` and `multistart-set1cR2L15genes`) and 97 genes (`multistart-set1cR2L97genes` and `multistart-set1cR2L97genes`) are run slightly differently. Since we expected the experiments with 97 genes to take a long time to finish, we wanted to run them on a remote server that we could log off of, so we wrote scripts to handle that.

## Instructions for running the Experiments with 97 genes 
This subsection describes how to run the experiements `multistart-set1cR2L97genes` and `multistart-set1cR2L97genes`. To do this, first copy the script `runL2R97.sh` (or `runR2L97.sh`) to the directory `/scripts/NetRAX`. Then run the script from that directory. Detailed instructions for running those scripts (especially on a remote server) can be found as comments in the scripts themselves. 

When the run has finished, we copy the output from `scripts/NetRAX/` to the experiment folder manually.

## Instructions for running the Experiments with 15 genes
This subsection describes how to run the experiements `multistart-set1cR2L15genes` and `multistart-set1cR2L15genes`.

After completing the setup and installing NetRAX, navigate to the directory `/scripts/NetRAX/` and do the following: 

First, run exactly one of the following lines, depending on whether you wish to run the L2R or the R2L experiment:

```
experiment_name="multistart-set1cL2R15genes"
experiment_name="multistart-set1cR2L15genes"
```

Then run all of the following:

```
path="../../analysis/netrax/${experiment_name}"
start_networks="${path}/starting-networks.treefile"
partition_file="${path}/partition.txt"
msa_file="${path}/set1c.fasta"

python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 42 --likelihood_type average  --name ${experiment_name}

```

When the run has finished, we copy the output from `scripts/NetRAX/` to the
experiment folder manually.

# Output and Plots
This analysis was run for 71 days, after which it was terminated prior to
completion of all 97 runs. Only 20 runs had been attempted at that point. Most
of the runs 0-18 took only a few days each, but the software stalled in run19
for 31 days.

Since the 97 genes run was terminated before completion, we have only 19
complete runs (out of 97), so I had to determine the best network manually going
through the logfiles and comparing the best networks from each of the 19
complete runs plus the 20th run which was terminated early. The best network and
its logl and bic score is at the end of the log file for each of the runs
0,1,...,18. Some of the runs appear to have not terminated correctly, but the
best network from that run is still recorded as the last entry in the logfile
(also it is located in the files with names like `run_n_inferred_network.nw`)
The very last run (the 20th run, called run19) was terminated prior to
completion; the best network for that last run is at the very end of the 54mb
file `multistart-set1cL2R97genes-script-output.txt`. It turned out that run19
was actually the run with the best bic score! 

The networks from each of the 20 runs can be found in plots.jl, which also
contains the code for plotting the best network.

This analysis was run on a really fast computer (magma2), on which NetRAX maxed
out 35 cores for the duration of the run. The lscpu results for magma2 are:

```

Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         46 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  64
  On-line CPU(s) list:   0-63
Vendor ID:               GenuineIntel
  Model name:            Intel(R) Xeon(R) CPU E5-2698 v3 @ 2.30GHz
    CPU family:          6
    Model:               63
    Thread(s) per core:  2
    Core(s) per socket:  16
    Socket(s):           2
    Stepping:            2
    CPU max MHz:         3600.0000
    CPU min MHz:         1200.0000
    BogoMIPS:            4594.74
    Flags:               fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush
                          dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_
                         tsc arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc cpuid aperfmperf 
                         pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pd
                         cm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx 
                         f16c rdrand lahf_lm abm cpuid_fault epb invpcid_single pti ssbd ibrs ibpb stibp
                          tpr_shadow vnmi flexpriority ept vpid ept_ad fsgsbase tsc_adjust bmi1 avx2 sme
                         p bmi2 erms invpcid cqm xsaveopt cqm_llc cqm_occup_llc dtherm ida arat pln pts 
                         md_clear flush_l1d
Virtualization features: 
  Virtualization:        VT-x
Caches (sum of all):     
  L1d:                   1 MiB (32 instances)
  L1i:                   1 MiB (32 instances)
  L2:                    8 MiB (32 instances)
  L3:                    80 MiB (2 instances)
NUMA:                    
  NUMA node(s):          2
  NUMA node0 CPU(s):     0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,
                         56,58,60,62
  NUMA node1 CPU(s):     1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,
                         57,59,61,63
Vulnerabilities:         
  Itlb multihit:         KVM: Mitigation: VMX disabled
  L1tf:                  Mitigation; PTE Inversion; VMX conditional cache flushes, SMT vulnerable
  Mds:                   Mitigation; Clear CPU buffers; SMT vulnerable
  Meltdown:              Mitigation; PTI
  Mmio stale data:       Mitigation; Clear CPU buffers; SMT vulnerable
  Retbleed:              Not affected
  Spec store bypass:     Mitigation; Speculative Store Bypass disabled via prctl and seccomp
  Spectre v1:            Mitigation; usercopy/swapgs barriers and __user pointer sanitization
  Spectre v2:            Mitigation; Retpolines, IBPB conditional, IBRS_FW, STIBP conditional, RSB filli
                         ng, PBRSB-eIBRS Not affected
  Srbds:                 Not affected
  Tsx async abort:       Not affected

```

# Instructions for setting up this experiment

This experiment was run in a similar manner to the experiment in the folder
`analysis/multistart-set1cL2R15genes/` and the readme in that folder may be
useful as it contains more details. We are just replicating that experiment here
with a different random seed, so to install NetRAX and set up the experimental
directories, follow the instructions theres.

# Instructions for running the experiment
The experiments are run similarly to the experiment for 97 genes in the L2R
case, i.e. using the script `runL2R15.sh`. Detailed instructions for the script
(especially on a remote server) can be found as comments in the script itself.



After completing the setup and installing NetRAX, navigate to the directory
`/scripts/NetRAX/` and do the following:

First, run
```
experiment_name="multistart-set1cL2R15genes-replicate"
```

Then run all of the following:

```
path="../../analysis/netrax/${experiment_name}"
start_networks="${path}/starting-networks.treefile"
partition_file="${path}/partition.txt"
msa_file="${path}/set1c.fasta"

time python3 netrax.py --msa_path ${msa_file} --partitions_path ${partition_file} --start_networks ${start_networks} --seed 41 --likelihood_type average  --name ${experiment_name}

```

Note the we have used a different seed from the original run `multistart-set1cL2R15genes`

When the run has finished, we copy the output from `scripts/NetRAX/` to the
experiment folder (i.e. this folder) manually.

# Plotting the network
The code to plot the best network is found in the file
`/analysis/netrax/multistart-set1cL2R15genes-replicate/plot.jl`


# Note about specs
This experiment was run on magma3

```
Last login: Sun Jun  4 23:02:54 2023 from 144.92.166.43
bacharach@rossby:~$ htop
bacharach@rossby:~$ lspcu
-bash: lspcu: command not found
bacharach@rossby:~$ lscpu
Architecture:            x86_64
  CPU op-mode(s):        32-bit, 64-bit
  Address sizes:         46 bits physical, 48 bits virtual
  Byte Order:            Little Endian
CPU(s):                  40
  On-line CPU(s) list:   0-39
Vendor ID:               GenuineIntel
  Model name:            Intel(R) Xeon(R) CPU E5-2660 v3 @ 2.60GHz
    CPU family:          6
    Model:               63
    Thread(s) per core:  2
    Core(s) per socket:  10
    Socket(s):           2
    Stepping:            2
    CPU max MHz:         3300.0000
    CPU min MHz:         1200.0000
    BogoMIPS:            5193.43
    Flags:               fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc arch_perfmon pebs bts re
                         p_good nopl xtopology nonstop_tsc cpuid aperfmperf pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_deadli
                         ne_timer aes xsave avx f16c rdrand lahf_lm abm cpuid_fault epb invpcid_single pti ssbd ibrs ibpb stibp tpr_shadow vnmi flexpriority ept vpid ept_ad fsgsbase tsc_adjust bmi1 avx2 smep bmi
                         2 erms invpcid cqm xsaveopt cqm_llc cqm_occup_llc dtherm ida arat pln pts md_clear flush_l1d
Virtualization features: 
  Virtualization:        VT-x
Caches (sum of all):     
  L1d:                   640 KiB (20 instances)
  L1i:                   640 KiB (20 instances)
  L2:                    5 MiB (20 instances)
  L3:                    50 MiB (2 instances)
NUMA:                    
  NUMA node(s):          2
  NUMA node0 CPU(s):     0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38
  NUMA node1 CPU(s):     1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39
Vulnerabilities:         
  Itlb multihit:         KVM: Mitigation: VMX disabled
  L1tf:                  Mitigation; PTE Inversion; VMX conditional cache flushes, SMT vulnerable
  Mds:                   Mitigation; Clear CPU buffers; SMT vulnerable
  Meltdown:              Mitigation; PTI
  Mmio stale data:       Mitigation; Clear CPU buffers; SMT vulnerable
  Retbleed:              Not affected
  Spec store bypass:     Mitigation; Speculative Store Bypass disabled via prctl and seccomp
  Spectre v1:            Mitigation; usercopy/swapgs barriers and __user pointer sanitization
  Spectre v2:            Mitigation; Retpolines, IBPB conditional, IBRS_FW, STIBP conditional, RSB filling, PBRSB-eIBRS Not affected
  Srbds:                 Not affected
  Tsx async abort:       Not affected
```

# Trilonet experiments
These experiments were run using TriLoNet Version 1.2 on Debian Buster (stable).

This folder contains the results of trilonet experiments. The paper for trilonet is https://academic.oup.com/mbe/article/33/8/2151/2578738

To setup trilonet, first get the software from https://www.uea.ac.uk/groups-and-centres/computational-biology/software/trilonet#2

Unzip TriLoNet3.zip into the folder `virus-project/scripts/`. If there are any
questions, follow the instructions in the file manual.pdf. The documentation is
good and it's easy to run.

For the following experiments, the output files are manually copied to the respective directories.
# Experiment 1 - set1c
Navigate to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`
then run

`java -jar TriLoNet.jar ../../../../analysis/trilonet/set1c.fasta set1c-output.dot set1c-output.txt`

The output files are `set1c-output.dot` and `set1c-output.txt`. 

# Experiment 2 - set1b
Same instructions as experiment 1, but instead this time run

`java -jar TriLoNet.jar ../../../../analysis/trilonet/set1b.fasta set1b-output.dot set1b-output.txt`

# Experiment 3 - set1c with kappa=4
Default value of kappa is 6.5. In this experiment, we repeat experiment 1 but
with kappa=4. I do not know what kappa is and am unable to find any
documentation about it, other than that it affects the inference of trinets from
the sequence data.
From  `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet` run
`java -jar TriLoNet.jar ../../../../analysis/trilonet/set1c.fasta --k4.0 set1c-kappa4-output.dot set1c-kappa4-output.txt`

# Experiment 4 - set1b with kappa=4
Repeat experiment 2 but with kappa=4.

From  `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet` run
`java -jar TriLoNet.jar ../../../../analysis/trilonet/set1b.fasta --k4.0 set1b-kappa4-output.dot set1b-kappa4-output.txt`

# Experiment 5 - set1c with breakpoints
Navigate to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`
then run
`taxonset="set1c"`
and 
`java -jar TriLoNet.jar ../../../../analysis/trilonet/${taxonset}.fasta --b60800,61000,81000,82700,109900,110300 ${taxonset}-with-breakpoints-output.dot ${taxonset}-with-breakpoints-output.txt`

The output files are `set1c-with-breakpoints-output.dot` and `set1c-with-breakpoints-output.txt`. 


# Experiment 6 - set1c with kappa=1
Starting with this experiment, we fully automate the process of making the experiment directory and outputting files to that director. 

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="6" # Experiment number
K="1.0" # Kappa value
taxonset="set1c" # or "set1b"
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

The output files are automatically saved to `outdir`.

    
# Experiment 7 - set1c with kappa=10

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="7" # Experiment number
K="10" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

The output files are automatically saved to `outdir`

    
    
# Experiment 8 - set1c with kappa=20

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="8" # Experiment number
K="20" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```
    
# Experiment 9 - set1c with kappa=4 and breakpoints

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="9" # Experiment number
K="4" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}-with-breakpoints"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --b60800,61000,81000,82700,109900,110300 --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

The output files are automatically saved to `outdir`
    
# Experiment 10 - set1b with kappa=4 and breakpoints

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="10" # Experiment number
K="4" # Kappa value
taxonset="set1b"
output_name="experiment-${N}--${taxonset}-kappa${K}-with-breakpoints"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --b60800,61000,81000,82700,109900,110300 --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

The output files are automatically saved to `outdir`

# Experiment 11 - set1c with kappa=4 and major breakpoints
Only breakpoints included here are those corresponding to the 1600bp BHV5 -> 216_II reticulation.

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="11" # Experiment number
K="4" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}-with-major-breakpoints"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --81000,82700 --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

The output files are automatically saved to `outdir`

# Experiment 12 - set1c with kappa=6.5 and major breakpoints
Only breakpoints included here are those corresponding to the 1600bp BHV5 -> 216_II reticulation.

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="12" # Experiment number
K="6.5" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}-with-major-breakpoints"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --81000,82700 --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

The output files are automatically saved to `outdir`

    
# Experiment 13 - set1c with kappa=2.0

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="13" # Experiment number
K="2" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```
# Experiment 14 - set1c with kappa=8.0

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="14" # Experiment number
K="8" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

# Experiment 15 - set1c with kappa=0.5

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="15" # Experiment number
K="0.5" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

The output files are automatically saved to `outdir`

# Experiment 16 - set1c with kappa=5

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```

N="16" # Experiment number
K="5" # Kappa value
taxonset="set1c"
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/${output_name}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt

```

The output files are automatically saved to `outdir`

# Experiment 17 - set1b with kappa=0.5,1,2,4,5,6.5,8,10,20

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```
taxonset="set1b"
N="17" # Experiment number
for K in 0.5 1 2 4 5 6.5 8 10 20
do
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/experiment-${N}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt
done

```
    

# Experiment 18 - set1c with kappa=0.5,1,2,4,5,6.5,8,10,20

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```
taxonset="set1c"
N="18" # Experiment number
for K in 0.5 1 2 4 5 6.5 8 10 20
do
output_name="experiment-${N}--${taxonset}-kappa${K}"
outdir="../../../../analysis/trilonet/experiment-${N}"
mkdir ${outdir}
indir="../../../../analysis/trilonet"
java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt
done

```



# Experiment 19 - set1c and set1b with kappa=4,6.5,8 and with breakpoints

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```
N="19" # Experiment number
indir="../../../../analysis/trilonet"
outdir="../../../../analysis/trilonet/experiment-${N}-with-breakpoints"
mkdir ${outdir}

for taxonset in "set1c" "set1b"
do
    for K in 4 6.5 8
    do
        output_name="experiment-${N}--${taxonset}-k${K}-with-breakpoints"
        java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --b60800,61000,81000,82700,109900,110300 --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt
    done
done
```


# Experiment 20 - set1c and set1b with kappa=4,6.5,8 with only major breakpoints

First go to `virus-project/scripts/TriLoNet3/TriLoNet/TriLoNet`

Then run the following code:

```
N="20" # Experiment number
indir="../../../../analysis/trilonet"
outdir="../../../../analysis/trilonet/experiment-${N}-with-major-breakpoints"
mkdir ${outdir}

for taxonset in "set1c" "set1b"
do
    for K in 4 6.5 8
    do
        output_name="experiment-${N}--${taxonset}-k${K}-with-major-breakpoints"
        java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --b81000,82700 --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt
    done
done
```




# Plots
Code for creating the plots is found in plot.jl

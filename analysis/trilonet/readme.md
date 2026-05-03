# Trilonet experiments
These experiments were run using TriLoNet Version 1.2 on Debian Buster (stable).

This folder contains the results of trilonet experiments and code for making the plots in the paper. The paper for trilonet is https://academic.oup.com/mbe/article/33/8/2151/2578738

To setup trilonet, first get the software from https://www.uea.ac.uk/groups-and-centres/computational-biology/software/trilonet#2

Unzip TriLoNet3.zip into the folder `virus-project/scripts/`. If there are any
questions, follow the instructions in the file manual.pdf. The documentation is
good and it's easy to run.

Experiments 21, 22, and 23 are used to create the plots in the paper.
Experiments 1-20 were exploratory; for these experiments, all output files
were manually copied to their respective directories.

# Experiment 1 - set1c
Navigate to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`
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
From  `virus-project/scripts/trilonet3/TriLoNet/TriLoNet` run
`java -jar TriLoNet.jar ../../../../analysis/trilonet/set1c.fasta --k4.0 set1c-kappa4-output.dot set1c-kappa4-output.txt`

# Experiment 4 - set1b with kappa=4
Repeat experiment 2 but with kappa=4.

From  `virus-project/scripts/trilonet3/TriLoNet/TriLoNet` run
`java -jar TriLoNet.jar ../../../../analysis/trilonet/set1b.fasta --k4.0 set1b-kappa4-output.dot set1b-kappa4-output.txt`

# Experiment 5 - set1c with breakpoints
Navigate to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`
then run
`taxonset="set1c"`
and 
`java -jar TriLoNet.jar ../../../../analysis/trilonet/${taxonset}.fasta --b60800,61000,81000,82700,109900,110300 ${taxonset}-with-breakpoints-output.dot ${taxonset}-with-breakpoints-output.txt`

The output files are `set1c-with-breakpoints-output.dot` and `set1c-with-breakpoints-output.txt`. 


# Experiment 6 - set1c with kappa=1
Starting with this experiment, we fully automate the process of making the experiment directory and outputting files to that director. 

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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

First go to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

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


<!-- # Experiment 21 - major breakpoint only -->

<!-- First download trilonet and extract to the `scripts/` directory. Then, in shell, -->
<!-- navigate to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet` -->

<!-- Then run the following code: -->

<!-- ``` -->
<!-- N="21" # Experiment number -->
<!-- suffix="with-major-breakpoint-only" -->
<!-- indir="../../../../analysis/trilonet" -->
<!-- outdir="../../../../analysis/trilonet/experiment-${N}-${suffix}" -->
<!-- mkdir ${outdir} -->

<!-- time for taxonset in "set1c" "set1b" -->
<!-- do -->
<!--     for K in {0.5,1,2,4,5,6.5,8,10,20} -->
<!--     do -->
<!--         output_name="experiment-${N}--${taxonset}-k${K}-${suffix}" -->
<!--             java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --b81000,82700 --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt -->
<!--     done -->
<!-- done -->
<!-- ``` -->
<!-- <\!-- real	13m4.432s -\-> -->
<!-- <\!-- user	13m28.788s -\-> -->
<!-- <\!-- sys	0m7.607s -\-> -->


<!-- # Experiment 22 - all conjectured breakpoints -->


<!-- First download trilonet and extract to the `scripts/` directory. Then, in shell, -->
<!-- navigate to `virus-project/scripts/trilonet3/TriLoNet/TriLoNet` -->

<!-- ``` -->
<!-- N="22" # Experiment number -->
<!-- suffix="with-all-conjectured-breakpoints" -->
<!-- indir="../../../../analysis/trilonet" -->
<!-- outdir="../../../../analysis/trilonet/experiment-${N}-${suffix}" -->
<!-- mkdir ${outdir} -->

<!-- for taxonset in "set1c" "set1b" -->
<!-- do -->
<!--     for K in {0.5,1,2,4,5,6.5,8,10,20} -->
<!--     do -->
<!--         output_name="experiment-${N}--${taxonset}-k${K}-${suffix}" -->
<!--         java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --b60800,61000,81000,82700,109900,110300  --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt -->
<!--     done -->
<!-- done -->
<!-- ``` -->

<!-- (Note that 127986 is the "cleaned sequences length" for set1c, and 127992 is -->
<!-- for set1b. This is why we were getting the out of bounds error in previous experiments.) -->

<!-- # Experiment 23 - no breakpoints supplied -->

<!-- ``` -->
<!-- N="23" # Experiment number -->
<!-- suffix="with-no-breakpoints" -->
<!-- indir="../../../../analysis/trilonet" -->
<!-- outdir="../../../../analysis/trilonet/experiment-${N}-${suffix}" -->
<!-- mkdir ${outdir} -->

<!-- for taxonset in "set1c" "set1b" -->
<!-- do -->
<!--     for K in {0.5,1,2,4,5,6.5,8,10,20} -->
<!--     do -->
<!--         output_name="experiment-${N}--${taxonset}-k${K}-${suffix}" -->
<!--         java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt -->
<!--     done -->
<!-- done -->
<!-- ``` -->

<!-- # Experiment 24 - uniform 1500 bp mesh -->

<!-- ``` -->
<!-- N="24" # Experiment number -->
<!-- suffix="with-uniform-1500-bp-mesh" -->
<!-- indir="../../../../analysis/trilonet" -->
<!-- outdir="../../../../analysis/trilonet/experiment-${N}-${suffix}" -->
<!-- mkdir ${outdir} -->

<!-- for taxonset in "set1c" "set1b" -->
<!-- do -->
<!--     for K in {0.5,1,2,4,5,6.5,8,10,20} -->
<!--     do -->
<!--         output_name="experiment-${N}--${taxonset}-k${K}-${suffix}" -->
<!--         java -jar TriLoNet.jar ${indir}/${taxonset}.fasta --b1500,3000,4500,6000,7500,9000,10500,12000,13500,15000,16500,18000,19500,21000,22500,24000,25500,27000,28500,30000,31500,33000,34500,36000,37500,39000,40500,42000,43500,45000,46500,48000,49500,51000,52500,54000,55500,57000,58500,60000,61500,63000,64500,66000,67500,69000,70500,72000,73500,75000,76500,78000,79500,81000,82500,84000,85500,87000,88500,90000,91500,93000,94500,96000,97500,99000,100500,102000,103500,105000,106500,108000,109500,111000,112500,114000,115500,117000,118500,120000,121500,123000,124500,126000,127500 --k${K} ${outdir}/${output_name}.dot ${outdir}/${output_name}.txt -->
<!--     done -->
<!-- done -->
<!-- ``` -->



# Experiments 21 - 24 (combined run)
Here we run several experiments at once, comparing sets 1b and 1c using (1)
different breakpoint choices and (2) different values of kappa. 

To replicate these experiments, first download trilonet and extract to the
`scripts/` directory. Then, in shell, navigate to
`virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

Then run the following code:

```
indir="../../../../analysis/trilonet"

for N in 21 22 23 24
do
    case "$N" in
        21)
            suffix="with-major-breakpoint-only"
            breakpoints="--b81000,82700"
            ;;
        22)
            suffix="with-all-conjectured-breakpoints"
            breakpoints="--b60800,61000,81000,82700,109900,110300"
            ;;
        23)
            suffix="with-no-breakpoints"
            breakpoints=""
            ;;
        24)
            suffix="with-uniform-1500-bp-mesh"
            breakpoints="--b1500,3000,4500,6000,7500,9000,10500,12000,13500,15000,16500,18000,19500,21000,22500,24000,25500,27000,28500,30000,31500,33000,34500,36000,37500,39000,40500,42000,43500,45000,46500,48000,49500,51000,52500,54000,55500,57000,58500,60000,61500,63000,64500,66000,67500,69000,70500,72000,73500,75000,76500,78000,79500,81000,82500,84000,85500,87000,88500,90000,91500,93000,94500,96000,97500,99000,100500,102000,103500,105000,106500,108000,109500,111000,112500,114000,115500,117000,118500,120000,121500,123000,124500,126000,127500"
            ;;
    esac
    outdir="../../../../analysis/trilonet/experiment-${N}-${suffix}"
    mkdir -p "$outdir"
    for taxonset in set1c set1b
    do
        for K in 0.5 1 2 4 5 6.5 8 10 20
        do
            output_name="experiment-${N}--${taxonset}-k${K}-${suffix}"
            java -jar TriLoNet.jar \
                "${indir}/${taxonset}.fasta" \
                $breakpoints \
                "--k${K}" \
                "${outdir}/${output_name}.dot" \
                "${outdir}/${output_name}.txt"
        done
    done
done
```

The above code saves the TriLoNet output files to four separate folders in the
directory `analysis/trilonet/`. Total runtime is about 45 minutes (my
desktop). The code for plotting the results of these experiments is found in
`analysis/trilonet/plots-for-experiments-21-through-24.jl`

Comment: The uniform breakpoints case (experiment 24) only goes up to 127500,
which is less than the full genome length. This is because 127986 is the
"cleaned sequences length" for set1c, and 127992 is for set1b. (This is why we
were getting the out of bounds error in previous experiments.)


# Dependence of TriLoNet on row order of MSA
This experiment uses a toy example to show that the output network of TriLoNet
depends on the order that species are listed in the MSA.

First, run the following from `virus-project/scripts/trilonet3/TriLoNet/TriLoNet`

```
dir="../../../../analysis/trilonet/trinet-row-test"
mkdir -p $dir
cat > $dir/row34.fasta <<'EOF'
>1
AAAAC
>2
CAAAC
>3
GAACA
>4
TAACA
EOF

cat > $dir/row43.fasta <<'EOF'
>1
AAAAC
>2
CAAAC
>4
TAACA
>3
GAACA
EOF

for name in row34 row43; do
    java -jar TriLoNet.jar "${dir}/${name}.fasta" --k6 "${dir}/${name}.dot" "${dir}/${name}.txt"
done

```

This gives two different networks, for row 34 and row43 respectively:

```
((((4)#H1,1),3),(#H1,2))root;
((((3)#H1,1),4),(#H1,2))root;
```

To visualize them, run the following Julia code from `$dir`:

```
using PhyloPlots, RCall, PhyloNetworks
net34 = PhyloPlots.readTopology("((((4)#H1,1),3),(#H1,2))root;")
net43 = PhyloPlots.readTopology("((((3)#H1,1),4),(#H1,2))root;")

PhyloNetworks.rotate!(net34,-3)
PhyloNetworks.rotate!(net34,-4)
PhyloNetworks.rotate!(net43,-3)
PhyloNetworks.rotate!(net43,-4)

R"png(filename='row34vsrow43.png', width=36, height=18, units='cm', res=300)"
R"par(mfrow=c(1,2), mar=c(1,1,3,1))"
PhyloPlots.plot(net34, tipcex=1.5, useedgelength=false, showedgenumber=false,
    shownodenumber=false, showgamma=true, arrowlen=0.1, style=:fulltree)
R"title('net34')"
PhyloPlots.plot(net43, tipcex=1.5, useedgelength=false, showedgenumber=false,
    shownodenumber=false, showgamma=true, arrowlen=0.1, style=:fulltree)
R"title('net43')"
R"dev.off()"

```

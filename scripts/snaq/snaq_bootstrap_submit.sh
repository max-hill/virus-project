# does 4 bootstrap reps, for set1b and set1c (so 8 total)
# they are done sequentially, so will use different seeds.
#
# to get 100 bootstrap replicates, execute this 25 times
# like this below, from ~/private/concordance/herpesvirus :
# ./gitrepo/scripts/snaq/snaq_bootstrap_submit.sh 0 &
# ./gitrepo/scripts/snaq/snaq_bootstrap_submit.sh 4 &
# ...
# ./gitrepo/scripts/snaq/snaq_bootstrap_submit.sh 96 &

startrep=$1
endrep=$(( $startrep + 3 ))
for (( i=$startrep; i<=$endrep; i++ ))
do
  for j in {b,c}
  do
    datainput="genetrees" # "SNPs"
    echo "set1$j, bootstrap rep $i, from $datainput"
    julia --project --color=yes gitrepo/scripts/snaq/snaq_bootstrap.jl $datainput $j $i
  done
done

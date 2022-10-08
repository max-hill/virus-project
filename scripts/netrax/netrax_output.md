# Interpreting *netrax* output

NetRAX supports 4 (rooted) network topology rearrangement moves: *rNNI*, *rSPR*, *arc insertion*, *arc removal*, and their respective undo-operations.

NetRAX iterates over these *move types* (in the following order) as it searches the space of networks.

| previous | next |
| :-: | :-: |
| arc removal | rSPR |
| rSPR | rNNI |
| rNNI | arc insertion |
| arc insertion | arc removal |

There may be several *move candidates* for each move type, and
the candidate that reduces the BIC the most is accepted. The topology is then updated, and the other parameters (e.g. branch lengths, reticulation probabilities) are optimized.

If no candidate reduces the BIC, then the next move type in the order is considered.

If the previous move accepted was **not** an arc insertion, then that same move type is (greedily) explored again instead of moving ahead in the order.

---

Here are some commands for extracting information from stdout when NetRAX is run with a set of starting networks. In our example, stdout has been piped to `netrax-shell-log.txt`.

To begin with, we might want to examine each run (which corresponds to one of the start networks) by itself but in more detail.

## Demarcate different runs
We used 15 start networks, and so we have output for 15 runs (run-0, run-1, ..., run-14) in a single file.

The output for each run tracks the sequence of move types explored and the set of candidates (and their respective scores) evaluated for each move. At the end of each run, a breakdown of how many moves of each type was taken is reported (keyword: "statistics"). We use this to extract the line ranges associated with each run.
```shell
$ # Breakdown of accepted moves
$ grep -i "statistics" netrax-shell-log.txt -n -A 4 |
> head -5
1896:Statistics on which moves were taken:
1897-RSPRMove: 26
1898-RNNIMove: 14
1899-ArcRemovalMove: 0
1900-ArcInsertionMove: 8

$ # Line no.s for demarcation
$ grep -i "statistics" netrax-shell-log.txt -n |
> grep -oE "[0-9]+" |
> head -2
1896:Statistics on which moves were taken:
3043:Statistics on which moves were taken:

$ # Accepted moves + scores for each run
$ sed -ne "1,1896p" netrax-shell-log.txt |
> grep -i "took" -n -A 1 |
> head -6
64: Took ArcInsertionMove
65-IMPROVED GLOBAL BEST SCORE FOUND SO FAR (1 reticulations): 519098.3824       
--
118: Took RNNIMove
119-IMPROVED GLOBAL BEST SCORE FOUND SO FAR (1 reticulations): 519055.9691      
--

$ # Each move taken is recorded with "took ..."
$ # This tallies with the breakdown above: 48 = 26 + 14 + 8
$ sed -ne "1,1896p" netrax-shell-log.txt |
> grep -i "took" |
> wc -l
48
```

## View no. of reticulations, log-likelihood and BIC for the best network for all runs.
```shell
$ grep -i "Best inferred network has" netrax-shell-log.txt |
> head -2
Best inferred network has 8 reticulations, logl = -257743.2999, bic = 518012.4033
Best inferred network has 6 reticulations, logl = -257753.9598, bic = 517922.0854

$ # Lines 1-15 are for the 15 runs
$ # Line 16 is for the "best" run (run-12 in our example) 
$ grep -i "Best inferred network has" netrax-shell-log.txt |
> wc -l
16
```

## Extract log-likelihoods and BICs for best networks (of each run)
```shell
$ # log-likelihoods
$ # sed to extract odd lines
$ grep -i "Best inferred network has" netrax-shell-log.txt |
> grep -oE "[- ][0-9]+\.[0-9]+" |
> sed -n "p;n" |
> head -2
-257743.2999
-257753.9598

$ # BICs
$ # sed to extract even lines
$ grep -i "Best inferred network has" netrax-shell-log.txt |
> grep -oE "[- ][0-9]+\.[0-9]+" |
> sed -n "n;p" |
> head -2
 518012.4033
 517922.0854
```

## Extract BIC trajectory for a single run
```shell
$ sed -ne "1,1896p" netrax-shell-log.txt |
> grep -i "improved global" -n |
> grep -oE "[0-9]+\.[0-9]+"
> | head -2
520260.7902
519098.3824

$ sed -ne "1,1896p" netrax-shell-log.txt |
> grep -i "improved global" -n |
> grep -oE "[0-9]+\.[0-9]+"
> | wc -l
50
```
Note that we have 2 more BICs (50) than no. of moves taken (48). This is because the BIC is computed for the best tree (0 reticulations) before any rearrangement move is made, and also an extra time for the best network (8 reticulations) using a more intensive optimization procedure ("slow mode").

### Example:

Suppose we want to compare the BICs after the nth and (n+1)th arc insertion moves for the first run. We extract the BIC values calculated after arc insertion moves and pipe them to `bic.txt`.

```shell
$ sed -ne "1,1896p" netrax-shell-log.txt |
> grep -i "took arcinsertionmove" -A 1 |
> grep -oE "[0-9]+\.[0-9]+" > bic.txt
```

These values can be read, say in *Julia*, and processed using array operations. As we can see, increasing the no. of reticulations yields somewhat diminishing marginal improvement in BIC.

```shell
julia> bic = map(x->parse(Float64,x),bic.txt)
8-element Array{Float64,1}:
 519098.3824
 518734.6868
 518633.4823
 518339.0693
 518157.1892
 518112.4854
 518073.0673
 518012.5455

julia> bic[2:8]-bic[1:7]
7-element Array{Float64,1}:
 -363.6955999999773
 -101.20450000005076
 -294.41300000000047
 -181.8800999999512
  -44.70380000001751
  -39.41810000000987
  -60.52179999998771
```
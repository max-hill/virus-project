library(ape)
library(stringr)

# Extract embedding cost profile from RF-Net output
embeddingcostProfile <- function(file){
  # read rfnet output: newick strings of "best" networks as the no. of reticulations increases
  rfnets <- readLines(file)
  r <- 0:(length(rfnets)-1) # no. of reticulations
  dist <- vector("numeric",length(rfnets)) # total embedding cost
  # extract metadata (no. of reticulations, total embedding cost) from each newick string
  for (i in 1:length(rfnets)) {
    net <- rfnets[i]
    metadata <- unlist(str_match_all(net,"[[:lower:]]+=[[[:digit:]]\\.]+"))
    r[i] <- strsplit(metadata[1],"=")[[1]][2]
    dist[i] <- strsplit(metadata[2],"=")[[1]][2]
  }
  # no. of reticulations, r vs total embedding cost, dist
  profile <- data.frame(r=r,dist=dist)
  return(profile)
}

# Save embedding cost profiles
# L2R, regular mode
write.table(
  embeddingcostProfile("results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick"),
  "results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted",quote=F)

# L2R, regular mode, fast heuristic
write.table(
  embeddingcostProfile("results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted_f.newick"),
  "results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_f",quote=F)

# L2R, tree-child mode
write.table(
  embeddingcostProfile("results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted_tc.newick"),
  "results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_tc",quote=F)

# L2R, tree-child mode, fast heuristic
write.table(
  embeddingcostProfile("results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted_tc_f.newick"),
  "results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_tc_f",quote=F)

# R2L, regular mode
write.table(
  embeddingcostProfile("results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick"),
  "results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted",quote=F)

# R2L, regular mode, fast heuristic
write.table(
  embeddingcostProfile("results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_f.newick"),
  "results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_f",quote=F)

# R2L, tree-child mode
write.table(
  embeddingcostProfile("results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_tc.newick"),
  "results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_tc",quote=F)

# R2L, tree-child mode, fast heuristic
write.table(
  embeddingcostProfile("results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_tc_f.newick"),
  "results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_tc_f",quote=F)
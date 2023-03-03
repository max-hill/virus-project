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

partition_sizes <- c("15genes_blksize10000","97genes_blksize1500")
partition_size <- partition_sizes[2]

# Save embedding cost profiles
# L2R, regular mode
write.table(
  embeddingcostProfile(sprintf("results/bhv/rfnet/%s_set1c_ogrooted.newick",partition_size)),
  sprintf("results/bhv/rfnet/embedding cost/%s_set1c_ogrooted",partition_size),quote=F)

# L2R, regular mode, fast heuristic
write.table(
  embeddingcostProfile(sprintf("results/bhv/rfnet/%s_set1c_ogrooted_f.newick",partition_size)),
  sprintf("results/bhv/rfnet/embedding cost/%s_set1c_ogrooted_f",partition_size),quote=F)

# L2R, tree-child mode
write.table(
  embeddingcostProfile(sprintf("results/bhv/rfnet/%s_set1c_ogrooted_tc.newick",partition_size)),
  sprintf("results/bhv/rfnet/embedding cost/%s_set1c_ogrooted_tc",partition_size),quote=F)

# L2R, tree-child mode, fast heuristic
write.table(
  embeddingcostProfile(sprintf("results/bhv/rfnet/%s_set1c_ogrooted_tc_f.newick",partition_size)),
  sprintf("results/bhv/rfnet/embedding cost/%s_set1c_ogrooted_tc_f",partition_size),quote=F)

# R2L, regular mode
write.table(
  embeddingcostProfile(sprintf("results/bhv/rfnet/%s_rev_set1c_ogrooted.newick",partition_size)),
  sprintf("results/bhv/rfnet/embedding cost/%s_rev_set1c_ogrooted",partition_size),quote=F)

# R2L, regular mode, fast heuristic
write.table(
  embeddingcostProfile(sprintf("results/bhv/rfnet/%s_rev_set1c_ogrooted_f.newick",partition_size)),
  sprintf("results/bhv/rfnet/embedding cost/%s_rev_set1c_ogrooted_f",partition_size),quote=F)

# R2L, tree-child mode
write.table(
  embeddingcostProfile(sprintf("results/bhv/rfnet/%s_rev_set1c_ogrooted_tc.newick",partition_size)),
  sprintf("results/bhv/rfnet/embedding cost/%s_rev_set1c_ogrooted_tc",partition_size),quote=F)

# R2L, tree-child mode, fast heuristic
write.table(
  embeddingcostProfile(sprintf("results/bhv/rfnet/%s_rev_set1c_ogrooted_tc_f.newick",partition_size)),
  sprintf("results/bhv/rfnet/embedding cost/%s_rev_set1c_ogrooted_tc_f",partition_size),quote=F)
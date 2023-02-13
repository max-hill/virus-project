library(ape)
library(stringr)

# rfnet output: "best" networks as the no. of reticulations increases
rfnets <- readLines("results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted.newick")

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

# save
write.table(profile,
            "results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted",
            quote=F)

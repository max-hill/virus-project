library(seqinr)
library(stringr)

# make gene MSAs to input to Bacter

# full msa for taxa of interest
full_fasta <- read.fasta(file="data/bhv/msa/set1c.fasta")
taxa <- names(full_fasta)
ntaxa <- length(taxa)

# partition file demarcating genes (may define genes L2R or R2L)
part <- read.table("data/bhv/part/15genes_blksize10000_rev.txt",sep=",",strip.white=T)
ngenes <- nrow(part)

for (i in 1:ngenes) {
  gene_msa <- vector("list",ntaxa) # gene msa for taxa of interest
  lo_up <- as.integer(
    str_split_1(str_split_1(part$V2[i],"=")[2],"-")) # lower/upper site positions
  for (j in 1:ntaxa) {
    gene_msa[[j]] <- full_fasta[[j]][lo_up[1]:lo_up[2]]
  }
  write.fasta(gene_msa,taxa,sprintf("data/bhv/msa/R2L/%s.fasta",
                                    str_split_1(part$V2[i],"=")[1]))
}

library(ape)
library(phangorn)

# outgroup rooted genetrees
ogrooted_gts <- read.tree("data/bhv/genetrees/ogrooted/15genes_blksize10000_rev_set1c.treefile")

# genetree embeddings
gt_embeddings <- read.tree("results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted.newick.embeddings.tre")

# genetree embeddings, tree-child restriction
gt_embeddings_tc <- read.tree("results/bhv/rfnet/15genes_blksize10000_rev_set1c_ogrooted_tc.newick.embeddings.tre")

# partition file demarcating genes
part <- read.table("data/bhv/part/15genes_blksize10000_rev.txt",sep=",",strip.white=T)

# compute embedding cost
cost <- vector("numeric",length(ogrooted_gts))
cost_tc <- vector("numeric",length(ogrooted_gts))
for (i in 1:length(ogrooted_gts)) {
  cost[i] <- RF.dist(ogrooted_gts[[i]],gt_embeddings[[i]],rooted=F)
  cost_tc[i] <- RF.dist(ogrooted_gts[[i]],gt_embeddings_tc[[i]],rooted=F)
}

# compare genetrees to embeddings
par(mar=c(2,2,2,0))
layout(matrix(1:4,nrow=2,ncol=2,byrow=T))

for (i in 1:length(ogrooted_gts)) {
  gt_em <- gt_embeddings[[i]]
  gt_em$tip.label[which(gt_em$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(gt_em,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
  
  rgt <- ogrooted_gts[[i]]
  rgt$tip.label[which(rgt$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(rgt,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}
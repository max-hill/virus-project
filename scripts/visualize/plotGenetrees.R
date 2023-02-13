library(ape)

# genetree embeddings
gt_embeddings <- read.tree("results/bhv/rfnet/15genes_blksize10000_set1c_ogrooted_r20_any.newick.embeddings.tre")

# partition file demarcating genes
part <- read.table("data/bhv/part/15genes_blksize10000.txt",sep=",",strip.white=T)

# plot embeddings
par(mar=c(2,2,2,0))
layout(matrix(1:4,nrow=2,ncol=2,byrow=T))
for (i in c(4,6,10,15)) {
  gt_em <- gt_embeddings[[i]]
  gt_em$tip.label[which(gt_em$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(gt_em,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}
for (i in c(1,3,5,7)) {
  gt_em <- gt_embeddings[[i]]
  gt_em$tip.label[which(gt_em$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(gt_em,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}
for (i in c(8,9,11,12)) {
  gt_em <- gt_embeddings[[i]]
  gt_em$tip.label[which(gt_em$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(gt_em,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}
for (i in c(13,14,2)) {
  gt_em <- gt_embeddings[[i]]
  gt_em$tip.label[which(gt_em$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(gt_em,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}

# plot consensus tree of rooted genetrees
par(mar=c(2,2,3,0))
layout(matrix(1:2,nrow=1,ncol=2,byrow=T))
emct100 <- consensus(gt_embeddings,rooted=T) # strict consensus
emct100$tip.label[which(emct100$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
plot(emct100)
mtext("Strict consensus",line=.5,cex=0.8)
emct50 <- consensus(gt_embeddings,p=0.5,rooted=T) # majority-rule consensus
emct50$tip.label[which(emct50$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
plot(emct50)
mtext("50% consensus",line=.5,cex=0.8)
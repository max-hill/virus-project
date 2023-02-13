library(ape)

# unrooted genetrees
unrooted_gts <- read.tree("data/bhv/genetrees/unrooted/15genes_blksize10000_rev_set1c.treefile")
ngenes <- length(unrooted_gts)
ogrooted_gts <- vector("list",ngenes)

# root genetrees
for (i in 1:ngenes) {
  ugt <- unrooted_gts[[i]]
  ogrooted_gts[[i]] <- root(ugt,"BHV5",resolve.root=T)
}
class(ogrooted_gts) <- "multiPhylo"
write.tree(ogrooted_gts,"data/bhv/genetrees/ogrooted/15genes_blksize10000_rev_set1c.treefile")

# partition file demarcating genes
part <- read.table("data/bhv/part/15genes_blksize10000_rev.txt",sep=",",strip.white=T)

# plot rooted genetrees
par(mar=c(2,2,2,0))
layout(matrix(1:4,nrow=2,ncol=2,byrow=T))
for (i in c(4,6,10,15)) {
  rgt <- ogrooted_gts[[i]]
  rgt$tip.label[which(rgt$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(rgt,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}
for (i in c(1,3,5,7)) {
  rgt <- ogrooted_gts[[i]]
  rgt$tip.label[which(rgt$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(rgt,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}
for (i in c(8,9,11,12)) {
  rgt <- ogrooted_gts[[i]]
  rgt$tip.label[which(rgt$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(rgt,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}
for (i in c(13,14,2)) {
  rgt <- ogrooted_gts[[i]]
  rgt$tip.label[which(rgt$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
  plot(rgt,use.edge.length=F)
  mtext(sub("=",", ",sprintf("%s",part$V2[i])),line=.5,cex=0.8)
}

# plot consensus tree of rooted genetrees
par(mar=c(2,2,3,0))
layout(matrix(1:2,nrow=1,ncol=2,byrow=T))
gct100 <- consensus(ogrooted_gts,rooted=T) # strict consensus
gct100$tip.label[which(gct100$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
plot(gct100)
mtext("Strict consensus",line=.5,cex=0.8)
gct50 <- consensus(ogrooted_gts,p=0.5,rooted=T) # majority-rule consensus
gct50$tip.label[which(gct50$tip.label == "Titanium_IBR_MLV_vaccine",arr.ind=T)] <- "Ti"
plot(gct50)
mtext("50% consensus",line=.5,cex=0.8)
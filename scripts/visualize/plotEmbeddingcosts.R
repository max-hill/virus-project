# RF-net embedding cost profiles for "15 genes, 10Kbp blocks" partition
# genes defined left to right, RF-net in regular mode
l2r <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted",
                  header=T)
# ... , RF-net in regular mode with fast heuristic
l2r_f <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_f",
                    header=T)
# ..., RF-net in tree-child mode
l2r_tc <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_tc",
                     header=T)
# ..., RF-net in tree-child mode with fast heuristic
l2r_tc_f <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_tc_f",
                       header=T)

# genes defined right to left, RF-net in regular mode
r2l <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted",
                  header=T)
# ..., RF-net in regular mode with fast heuristic
r2l_f <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted_f",
                    header=T)
# ..., RF-net in tree-child mode
r2l_tc <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted_tc",
                     header=T)
# ..., RF-net in tree-child mode with fast heuristic
r2l_tc_f <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted_tc_f",
                       header=T)

# Calculate how total embedding cost decreases with the no. of reticulations
l2r_dec <- abs(l2r$dist[2:nrow(l2r)]-l2r$dist[1:(nrow(l2r)-1)])
l2r_f_dec <- abs(l2r_f$dist[2:nrow(l2r_f)]-l2r_f$dist[1:(nrow(l2r_f)-1)])
l2r_tc_dec <- abs(l2r_tc$dist[2:nrow(l2r_tc)]-l2r_tc$dist[1:(nrow(l2r_tc)-1)])
l2r_tc_f_dec <- abs(l2r_tc_f$dist[2:nrow(l2r_tc_f)]-l2r_tc_f$dist[1:(nrow(l2r_tc_f)-1)])

r2l_dec <- abs(r2l$dist[2:nrow(r2l)]-r2l$dist[1:(nrow(r2l)-1)])
r2l_f_dec <- abs(r2l_f$dist[2:nrow(r2l_f)]-r2l_f$dist[1:(nrow(r2l_f)-1)])
r2l_tc_dec <- abs(r2l_tc$dist[2:nrow(r2l_tc)]-r2l_tc$dist[1:(nrow(r2l_tc)-1)])
r2l_tc_f_dec <- abs(r2l_tc_f$dist[2:nrow(r2l_tc_f)]-r2l_tc_f$dist[1:(nrow(r2l_tc_f)-1)])

###############################################################################

pdf(file="figures/TEC.pdf",width=12,height=6)
# Plot embedding cost profiles for different RF-net modes
par(oma=c(2,2,0,0),mar=c(3,3,3,2)+0.1)
layout(matrix(c(1,2),1,2,byrow=T))
plot(x=1:10,y=l2r$dist,type="l",xlab="",ylab="",main="L2R",
     col="red",lwd=5,xaxt="n",cex.axis=1.5)
axis(1,at=1:10,cex.axis=1.5)
lines(x=1:9,y=l2r_f$dist,col="red",lty="dashed",lwd=5)
lines(x=1:9,y=l2r_tc$dist,col="blue")
lines(x=1:9,y=l2r_tc_f$dist,col="blue",lty="dashed")

plot(x=1:10,y=r2l$dist,type="l",xlab="",ylab="",main="R2L",
     col="red",lwd=5,xaxt="n",cex.axis=1.5)
axis(1,at=1:10,cex.axis=1.5)
lines(x=1:9,y=r2l_f$dist,col="red",lty="dashed",lwd=5)
lines(x=1:9,y=r2l_tc$dist,col="blue")
lines(x=1:9,y=r2l_tc_f$dist,col="blue",lty="dashed")
mtext(text="Total embedding cost",side=2,line=0,outer=T,cex=1.5)
mtext(text="No. of reticulations",side=1,line=0,outer=T,cex=1.5)

legend(x=5,y=23,legend=c("reg","reg + f","tc","tc + f"),
       col=c("red","red","blue","blue"),
       lty=c("solid","dashed","solid","dashed"),bty="n",cex=1.5,lwd=3)
dev.off()

###############################################################################

pdf(file="figures/TEC_decrease.pdf",width=12,height=6)
# Plot how total embedding cost decreases with the no. of reticulations
par(oma=c(2,2,0,0),mar=c(3,3,3,2)+0.1)
layout(matrix(c(1,2),1,2,byrow=T))
plot(x=1:9,y=l2r_dec,type="l",xlab="",ylab="",main="L2R",
     col="red",lwd=5,ylim=c(0,6),xaxt="n",cex.axis=1.5)
axis(1,at=1:9,cex.axis=1.5)
lines(x=1:8,y=l2r_f_dec,col="red",lty="dashed",lwd=5)
lines(x=1:8,y=l2r_tc_dec,col="blue")
lines(x=1:8,y=l2r_tc_f_dec,col="blue",lty="dashed")

plot(x=1:9,y=r2l_dec,type="l",xlab="",ylab="",main="R2L",
     col="red",lwd=5,ylim=c(0,6),xaxt="n",cex.axis=1.5)
axis(1,at=1:9,cex.axis=1.5)
lines(x=1:8,y=r2l_f_dec,col="red",lty="dashed",lwd=5)
lines(x=1:8,y=r2l_tc_dec,col="blue")
lines(x=1:8,y=r2l_tc_f_dec,col="blue",lty="dashed")
mtext(text="Decrease in total embedding cost",side=2,line=0,outer=T,cex=1.5)
mtext(text="No. of reticulations",side=1,line=0,outer=T,cex=1.5)

legend(x=5,y=6,legend=c("reg","reg + f","tc","tc + f"),
       col=c("red","red","blue","blue"),
       lty=c("solid","dashed","solid","dashed"),bty="n",cex=1.5,lwd=3)
dev.off()
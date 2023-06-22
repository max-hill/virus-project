library(stringr)

plot_detected_segments <- function(paths,ylab,main) {
  plot(x=0,xlim=c(0,144551),ylim=c(1,length(paths)),type="n",xlab="",ylab="",
       xaxt="n",yaxt="n")
  title(main=main,line=1,cex.main=1.5)
  axis(side=1,at=seq(from=10000,to=140000,by=10000),
       labels=format(seq(from=10000,to=140000,by=10000),scientific=T),
       cex.axis=1.5)
  for (i in 1:length(paths)) {
    # 34 is a magic number! It works only for the NEXUS files for set1c.
    nexusinfo_str <- readLines(paths[i])[34]
    start_end_pts <- data.frame(str_match_all(nexusinfo_str,"region=\\{(?<start>\\d+),(?<end>\\d+)\\}")[[1]])
    segments(x0=as.numeric(start_end_pts$start),y0=i,
             x1=as.numeric(start_end_pts$end),y1=i,lwd=4,
             col=rgb(0,0,0,0.3))
  }
}

# posterior support threshold 50%
main50 <- expression("Posterior support" >= "50%")
ylabs50 <- c("1","2","3","4","5","4a","5a","5b")
path50 <- "results/bhv/bacter/experimentMC3_rep%s/mc3.set1c.support50_burnin10.tree"
paths50 <- sapply(X=ylabs50,FUN=function(ylab) sprintf(path50,ylab))
paths50[1] <- "results/bhv/bacter/experimentMC3/mc3.set1c.support50_burnin10.tree"

# posterior support threshold 20%
main20 <- expression("Posterior support" >= "20%")
ylabs20 <- c("1","2","3","4","5","4a","5a","5b")
path20 <- "results/bhv/bacter/experimentMC3_rep%s/mc3.set1c.support20_burnin10.tree"
paths20 <- sapply(X=ylabs20,FUN=function(ylab) sprintf(path20,ylab))
paths20[1] <- "results/bhv/bacter/experimentMC3/mc3.set1c.support20_burnin10.tree"

pdf(file='figures/bacter/detectedsegments_support_50_20.pdf',width=12,height=6)
par(oma=c(2,2.2,0,0),mar=c(3,3,3,1)+0.1)
layout(matrix(c(1,2),1,2,byrow=T))
plot_detected_segments(paths50,ylab50,main50)
axis(side=2,at=1:length(paths50),labels=ylabs50,cex.axis=1.5,las=1)
plot_detected_segments(paths20,ylab20,main20)
mtext(text="Runs",side=2,line=1,outer=T,cex=1.5)
mtext(text="Site no.",side=1,line=0,outer=T,cex=1.5)
dev.off()